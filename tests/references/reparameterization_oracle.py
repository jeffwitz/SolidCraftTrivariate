#!/usr/bin/env python3
"""Independent scalar-QP and map-jet reference for the two-brick P6 model."""
from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys
import platform


BOUND = 0.75
ROOT = pathlib.Path(__file__).resolve().parents[2]


def wire_float(value):
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
    return value


def max_error(expected, actual):
    if expected is None or actual is None:
        return 0.0 if expected is None and actual is None else math.inf
    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        return abs(expected-actual)
    if isinstance(expected, list) and isinstance(actual, list) and len(expected)==len(actual):
        return max((max_error(x,y) for x,y in zip(expected,actual)),default=0.0)
    return 0.0 if expected==actual else math.inf


def call_rust_probe(payload: list[dict]) -> list[dict]:
    subprocess.run(["rtk","cargo","build","--release","--locked","--offline","--example","p6b_oracle_probe"],
                   cwd=ROOT,check=True,capture_output=True,text=True)
    binary=ROOT/"target/release/examples/p6b_oracle_probe"
    completed=subprocess.run([str(binary),"p6"],input=json.dumps(payload,allow_nan=False),
                             cwd=ROOT,check=True,capture_output=True,text=True)
    return json.loads(completed.stdout)


def solve(theta: float, weights: tuple[float, float], regularization: float) -> dict:
    if (not math.isfinite(theta) or any(not math.isfinite(w) or w < 0 for w in weights)
            or not math.isfinite(regularization) or regularization < 0):
        return {"accepted": False, "reason": "invalid input"}
    alpha = tuple(w / 3.0 for w in weights)
    if any(w > 0.0 and a == 0.0 for w, a in zip(weights, alpha)):
        return {"accepted": False, "reason": "data weight below resolvable objective range"}
    curvature = sum(alpha)
    if any(w > 0.0 for w in weights) and curvature == 0.0:
        return {"accepted": False, "reason": "data curvature below resolvable range"}
    denom = curvature + 2.0 * regularization
    if denom == 0:
        return {"accepted": True, "unique": False, "regime": "NonUnique",
                "a_range": [-BOUND, BOUND], "derivative": None, "kkt": 0.0}
    if not math.isfinite(denom):
        return {"accepted": False, "reason": "overflow"}
    response = curvature / denom
    if curvature > 0.0 and response == 0.0:
        return {"accepted": False, "reason": "data response below resolvable range"}
    raw = response * theta
    if not math.isfinite(raw):
        return {"accepted": False, "reason": "overflow"}
    if curvature > 0.0 and not math.isfinite(BOUND / response):
        return {"accepted": False, "reason": "activation threshold exceeds finite range"}
    tol = 32.0 * sys.float_info.epsilon * max(BOUND, abs(raw), 1.0)
    dist = abs(abs(raw) - BOUND)
    if abs(raw) == BOUND:
        regime = "ActivationThreshold"
        a = raw
        derivative = None
    elif raw > BOUND:
        regime = "NumericalUncertainty" if dist <= tol else "UpperBound"
        a = BOUND
        derivative = None if regime == "NumericalUncertainty" else 0.0
    elif raw < -BOUND:
        regime = "NumericalUncertainty" if dist <= tol else "LowerBound"
        a = -BOUND
        derivative = None if regime == "NumericalUncertainty" else 0.0
    else:
        regime = "NumericalUncertainty" if dist <= tol else "Interior"
        a = raw
        derivative = response if regime == "Interior" else None

    grad = 2.0 * curvature * (a - theta) + 4.0 * regularization * a
    grad1 = 2.0 * alpha[1] * (-a + theta) - 2.0 * regularization * a
    equality_multiplier = -grad1
    upper = max(-grad, 0.0) if a == BOUND and raw > BOUND else 0.0
    lower = max(grad, 0.0) if a == -BOUND and raw < -BOUND else 0.0
    stationarity = grad + upper - lower
    try:
        obj = alpha[0] * (a-theta)**2 + alpha[1] * (-a+theta)**2 + 2*regularization*a*a
    except OverflowError:
        return {"accepted": False, "reason": "objective overflow"}
    if not math.isfinite(obj):
        return {"accepted": False, "reason": "objective overflow"}
    return {"accepted": True, "unique": True, "regime": regime,
            "coefficients": [a, -a], "raw": raw, "response": response,
            "activation_threshold": (BOUND/response if curvature > 0.0 else None),
            "objective": obj, "equality_multiplier":equality_multiplier,
            "upper_multiplier": upper,
            "lower_multiplier": lower, "tangent_curvature": 2*denom,
            "kkt": max(abs(stationarity), abs(upper*(BOUND-a)),
                       abs(lower*(BOUND+a)), max(abs(a)-BOUND, 0.0)),
            "derivative": None if derivative is None else [derivative, -derivative]}


def map_value(result: dict, brick: int, u: float) -> float:
    a = result["coefficients"][brick]
    return brick + u + a*u*(1-u)


def gauss_legendre(order: int):
    pairs=[]
    for i in range(1,(order+1)//2+1):
        x=math.cos(math.pi*(i-0.25)/(order+0.5))
        for _ in range(32):
            p0,p1=1.0,x
            for k in range(2,order+1):
                p0,p1=p1,((2*k-1)*x*p1-(k-1)*p0)/k
            derivative=order*(x*p1-p0)/(x*x-1.0)
            step=p1/derivative
            x-=step
            if abs(step)<2e-16: break
        weight=2.0/((1.0-x*x)*derivative*derivative)
        pairs.extend(((-x,weight),(x,weight)) if 2*i<=order else ((x,weight),))
    return [(0.5*(x+1.0),0.5*w) for x,w in sorted(pairs)]


def nonaffine_physical_integral(theta: float, order: int) -> float:
    """Pull back a physical moment from one fixed curved, sheared solid."""
    solution=solve(theta,(1.0,2.0),0.1)
    points=gauss_legendre(order)
    total=0.0
    alpha,beta,gamma=0.07,0.11,0.05
    for brick,a in enumerate(solution["coefficients"]):
        for u,wu in points:
            source=brick+u+a*u*(1.0-u)
            speed=1.0+a*(1.0-2.0*u)
            jac=(1.0+beta*source)*speed
            for v,wv in points:
                y=v+alpha*source*(2.0-source)+beta*source*v
                for w,ww in points:
                    z=w+gamma*source*v
                    total+=wu*wv*ww*jac*(source*source+y*z)
    return total


def nonaffine_physical_sensitivity() -> dict:
    result={}
    for order in (8,12):
        values={f"{theta:.6f}":nonaffine_physical_integral(theta,order)
                for theta in (0.2-1e-3,0.2,0.2+1e-3)}
        derivative=(values["0.201000"]-values["0.199000"])/0.002
        result[str(order)]={"values":values,"centered_derivative":derivative,
                            "max_variation":max(values.values())-min(values.values())}
    return {"mapping":"(s, v+0.07*s*(2-s)+0.11*s*v, w+0.05*s*v)",
            "observable":"integral (X^2+Y*Z) dV over fixed image",
            "expected_sensitivity":0.0,"quadrature":result}


def run() -> dict:
    cases = [
        ("asymmetric_interior", 0.2, (1.0, 2.0), 0.1),
        ("symmetric_interior", 0.4, (1.0, 1.0), 0.0),
        ("upper_strict", 1.0, (1.0, 1.0), 0.0),
        ("lower_strict", -1.0, (1.0, 1.0), 0.0),
        ("upper_threshold", 0.75, (1.0, 1.0), 0.0),
        ("lower_threshold", -0.75, (1.0, 1.0), 0.0),
        ("below_threshold_ulp", math.nextafter(0.75, 0.0), (1.0, 1.0), 0.0),
        ("above_threshold_ulp", math.nextafter(0.75, math.inf), (1.0, 1.0), 0.0),
        ("regularization_only", 100.0, (0.0, 0.0), 1.0),
        ("flat", 0.4, (0.0, 0.0), 0.0),
        ("single_data_side", 0.3, (0.0, 2.0), 0.2),
    ]
    rows = []
    all_inputs = cases + [
        ("nan_theta", math.nan, (1.0,1.0), 0.0),
        ("negative_weight",0.0,(-1.0,1.0),0.0),
        ("infinite_regularizer",0.0,(1.0,1.0),math.inf),
        ("subnormal_weight",0.0,(math.ulp(0.0),0.0),0.0),
        ("finite_product_overflow",1e308,(1.0,1.0),0.0),
    ]
    payload=[{"case":name,"theta":wire_float(theta),"weights":[wire_float(w) for w in weights],
              "regularization":wire_float(lam)} for name,theta,weights,lam in all_inputs]
    native=call_rust_probe(payload)
    if len(native)!=len(all_inputs):
        raise AssertionError("Rust P6 probe returned a different case count")
    native_by_name={x["case"]:x for x in native}
    max_kkt = 0.0
    max_jet_error = 0.0
    native_max_errors={}
    for name, theta, weights, lam in cases:
        res = solve(theta, weights, lam)
        row = {"case": name, **res}
        rust=native_by_name[name]
        row["rust"]={k:rust[k] for k in ("accepted","error","unique","regime","coefficients","derivative",
                                               "objective","kkt_residual","tangent_curvature","equality_rank",
                                               "activation_threshold","minimum_map_speed","map_values","map_jets") if k in rust}
        if rust["accepted"] != res["accepted"]:
            raise AssertionError(f"P6 acceptance mismatch: {name}")
        if res["accepted"]:
            if rust["unique"] != res["unique"]:
                raise AssertionError(f"P6 uniqueness mismatch: {name}")
            if "regime" in res and rust["regime"] != res["regime"]:
                raise AssertionError(f"P6 regime mismatch: {name}: {rust['regime']} vs {res['regime']}")
            for ref_key,native_key in (("coefficients","coefficients"),("derivative","derivative"),
                                      ("objective","objective"),("kkt","kkt_residual"),
                                      ("tangent_curvature","tangent_curvature"),
                                      ("activation_threshold","activation_threshold"),
                                      ("minimum_map_speed","minimum_map_speed")):
                if ref_key in res:
                    err=max_error(res[ref_key],rust.get(native_key))
                    native_max_errors[ref_key]=max(native_max_errors.get(ref_key,0.0),err)
                    if err>1e-12:
                        raise AssertionError(f"P6 native/reference mismatch {name} {ref_key}: {err}")
            if res.get("coefficients") is not None:
                expected_maps=[map_value(res,b,u) for b in (0,1) for u in (0.125,0.37,0.5,0.83)]
                map_error=max_error(expected_maps,rust.get("map_values"))
                native_max_errors["map_values"]=max(native_max_errors.get("map_values",0.0),map_error)
                if map_error>1e-12:
                    raise AssertionError(f"P6 map mismatch: {name}")
                expected_jets=(None if res.get("derivative") is None else
                    [res["derivative"][b]*u*(1-u) for b in (0,1) for u in (0.125,0.37,0.5,0.83)])
                jet_error=max_error(expected_jets,rust.get("map_jets"))
                native_max_errors["map_jets"]=max(native_max_errors.get("map_jets",0.0),jet_error)
                if jet_error>1e-12:
                    raise AssertionError(f"P6 map jet mismatch: {name}")
        max_kkt = max(max_kkt, res.get("kkt", 0.0))
        if res.get("derivative") is not None:
            errors = []
            for h in (1e-4, 5e-5, 2.5e-5):
                plus, minus = solve(theta+h, weights, lam), solve(theta-h, weights, lam)
                for brick in (0, 1):
                    for u in (0.125, 0.37, 0.5, 0.83):
                        fd = (map_value(plus, brick, u)-map_value(minus, brick, u))/(2*h)
                        exact = res["derivative"][brick]*u*(1-u)
                        errors.append(abs(fd-exact))
            row["centered_map_jet_error_max"] = max(errors)
            max_jet_error = max(max_jet_error, max(errors))
        if res.get("coefficients") is not None:
            row["minimum_map_speed"] = 1.0-max(abs(v) for v in res["coefficients"])
            row["endpoint_error"] = max(abs(map_value(res,b,0.0)-float(b))
                                         for b in (0,1))
        rows.append(row)
    invalid_rows = []
    for name,theta,weights,lam in all_inputs[len(cases):]:
        ref=solve(theta,weights,lam); rust=native_by_name[name]
        if rust["accepted"] != ref["accepted"]:
            raise AssertionError(f"P6 invalid-case admission mismatch: {name}")
        invalid_rows.append({"case":name,**ref,"rust":rust})
    return {"schema":"p6_reparameterization_oracle_v1",
            "provenance":{"source_baseline":"ff2c36f","checkout_revision":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
                          "python":platform.python_version(),"protocol":"docs/lot-p6b-oracle.md",
                          "commands":["rtk python3 tests/references/reparameterization_oracle.py",
                                      "rtk cargo build --release --locked --offline --example p6b_oracle_probe",
                                      "rtk cargo test --locked --offline reparameterization"]},
            "case_count": len(rows), "invalid_count": len(invalid_rows),
            "native_probe_case_count":len(native),"native_reference_comparison":"passed",
            "max_native_reference_absolute_errors":native_max_errors,
            "max_kkt_residual": max_kkt, "max_centered_map_jet_error": max_jet_error,
            "acceptance":{"kkt_residual_max":1e-12,"centered_map_jet_error_max":2e-8},
            "nonaffine_physical_sensitivity":nonaffine_physical_sensitivity(),
            "cases": rows, "invalid_cases": invalid_rows}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
