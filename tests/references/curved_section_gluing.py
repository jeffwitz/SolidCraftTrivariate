"""Independent symbolic gluing classification for the actual polynomial pocket sections."""

from pathlib import Path
import argparse
import json

import sympy as s

U, V, T = s.symbols("u v t")
ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "docs/site/data/plate_blind_pocket/plate_blind_pocket.model.json"


def chart(block):
    du, dv, dw = block["degrees"]
    assert [du, dv, dw] == [2, 2, 1]
    net = block["controls"]
    assert len(net) == 18 and len({x[-1] for x in net}) == 1 and net[0][-1] > 0
    for i in range(0, 18, 2):
        assert net[i][:2] == net[i+1][:2] and net[i+1][2] > net[i][2]
    B = lambda x, i: s.binomial(2,i) * x**i * (1-x)**(2-i)
    return s.Matrix([sum(B(U,i)*B(V,j)*s.Rational(net[6*i+2*j][k], net[6*i+2*j][3])
                         for i in range(3) for j in range(3)).expand() for k in range(2)])


def det(a,b):
    return s.expand(a[0]*b[1] - a[1]*b[0])


def jets(x, side, reverse=False):
    normal, tang = (U,V) if side[0] == "u" else (V,U)
    end = int(side[1])
    edge = {normal:end, tang:1-T if reverse else T}
    xy = x.subs(edge).applyfunc(s.expand)
    inward = ((1-2*end)*x.diff(normal)).subs(edge).applyfunc(s.expand)
    return xy, xy.diff(T), inward


def witness(left, right, first, second, rev):
    xy, tangent, a = jets(left, first)
    yz, tb, b = jets(right, second, rev)
    assert xy == yz and tangent == tb
    # Inward transverse directions have opposite signed area relative to shared tangent.
    al, ar = det(a,tangent), det(b,tangent)
    g = s.gcd(al,ar)
    assert g != 0
    l, r = [s.cancel(z/g) for z in (al,ar)]
    # r*a - l*b + gamma*tangent = 0.
    k = 0 if tangent[0] != 0 else 1
    gamma = s.cancel((l*b[k]-r*a[k])/tangent[k])
    assert all(s.cancel(z) == 0 for z in r*a-l*b+gamma*tangent)
    def degree(z):
        return 0 if z == 0 else int(s.degree(z,T))
    def affine(x):
        return s.denom(s.cancel(x)).is_number and degree(s.cancel(x)) <= 1
    reason = "coprime transverse factors exceed degree one"
    accepted = False
    if affine(l) and affine(r):
        assert l.subs(T,0)*l.subs(T,1)>0 and r.subs(T,0)*r.subs(T,1)>0
        assert (l*r).subs(T,0)<0
        q = s.symbols("a:4")
        delta = l*(q[2]+q[3]*T)-r*(q[0]+q[1]*T)-gamma
        co = s.together(delta).as_numer_denom()[0].expand().as_poly(T).all_coeffs()
        sols = s.linsolve(co,q)
        accepted = bool(sols)
        reason = "linear gluing factors verified" if accepted else "no affine beta split"
    return {"as_g1": accepted, "reason":reason, "alpha_left":str(l), "alpha_right":str(r),
            "gamma":str(gamma), "normalized_affine_either_direction": bool(
            (affine(l/r) and affine(gamma/r)) or (affine(r/l) and affine(gamma/l)))}


def classify(model, native=None):
    charts = {b["id"]: chart(b) for b in model["blocks"]}
    rows = []
    records = {} if native is None else {tuple(tuple(s) for s in r["faces"]): r for r in native}
    for face in model["topology"]["interfaces"]:
        (a,first),(b,second) = face["first"], face["second"]
        if first[0] == "w":
            continue
        orient = face["orientation"]
        assert not orient["swap"] and not orient["flip_second"]
        rows.append({"first":face["first"],"second":face["second"],
                     **witness(charts[a], charts[b], first, second, orient["flip_first"])})
        if native is not None:
            record = records.pop(((a,first),(b,second)))
            g = record["as_g1"]
            assert g["status"] == ("as_g1" if rows[-1]["as_g1"] else "not_as_g1")
            if rows[-1]["as_g1"]:
                def polynomial(name):
                    return sum(s.Rational(c)*T**i for i,c in enumerate(g["gluing_factors"][name]))
                l,r,beta = map(polynomial, ("alpha_first","alpha_second","beta"))
                bl,br = map(polynomial, ("beta_first","beta_second"))
                assert all(s.degree(x,T) <= 1 for x in (l,r,bl,br))
                assert s.expand(l*br+r*bl-beta) == 0
                _, tangent, dl = jets(charts[a],first)
                _, _, dr = jets(charts[b],second,orient["flip_first"])
                assert all(s.expand(v)==0 for v in l*dl+r*dr+beta*tangent), (a,b,l,r,beta)
    if native is not None:
        assert all(face[0][1][0]=="w" for face in records)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native", type=Path, help="verify all native witnesses on the analytic jets")
    args = parser.parse_args()
    data = json.dumps(classify(json.loads(MODEL.read_text()),
                               json.loads(args.native.read_text()) if args.native else None), indent=2) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")
