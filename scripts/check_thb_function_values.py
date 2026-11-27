#!/usr/bin/env python3
"""Independently check a serialized ``thb_strip_exact_v1`` field basis.

The oracle evaluates tensor B-splines with Cox--de Boor recursion, derives
uniform knot-insertion coefficients, and reconstructs the truncated
hierarchical functions from the active support rule. It does not import or
call the Rust extraction implementation.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import itertools
import math
from fractions import Fraction
from pathlib import Path


def frac(value):
    if isinstance(value, list) and len(value) == 2:
        return Fraction(int(value[0]), int(value[1]))
    return Fraction(value)


def tensor_shape(knots, degree=3):
    return [len(axis) - degree - 1 for axis in knots]


def unravel(index, shape):
    k = index % shape[2]
    j = (index // shape[2]) % shape[1]
    i = index // (shape[1] * shape[2])
    return i, j, k


def flatten(ijk, shape):
    return (ijk[0] * shape[1] + ijk[1]) * shape[2] + ijk[2]


def basis_values(knots, degree, x):
    """Return all degree-p basis values and derivatives at x."""
    count = len(knots) - degree - 1
    end = float(knots[-degree - 1])
    # Evaluate the left limit at the clamped right endpoint so the derivative
    # follows the same closed-span convention as Bernstein extraction.
    if x == end:
        x = math.nextafter(end, -math.inf)
    previous = [
        1.0 if float(knots[i]) <= x < float(knots[i + 1]) else 0.0
        for i in range(len(knots) - 1)
    ]
    for p in range(1, degree + 1):
        current = [0.0] * (len(knots) - p - 1)
        for i in range(len(current)):
            left_den = float(knots[i + p] - knots[i])
            right_den = float(knots[i + p + 1] - knots[i + 1])
            if left_den:
                current[i] += (x - float(knots[i])) / left_den * previous[i]
            if right_den:
                current[i] += (float(knots[i + p + 1]) - x) / right_den * previous[i + 1]
        previous = current
    values = previous[:count]

    lower = [
        1.0 if float(knots[i]) <= x < float(knots[i + 1]) else 0.0
        for i in range(len(knots) - 1)
    ]
    for p in range(1, degree):
        current = [0.0] * (len(knots) - p - 1)
        for i in range(len(current)):
            left_den = float(knots[i + p] - knots[i])
            right_den = float(knots[i + p + 1] - knots[i + 1])
            if left_den:
                current[i] += (x - float(knots[i])) / left_den * lower[i]
            if right_den:
                current[i] += (float(knots[i + p + 1]) - x) / right_den * lower[i + 1]
        lower = current
    derivatives = [0.0] * count
    for i in range(count):
        left_den = float(knots[i + degree] - knots[i])
        right_den = float(knots[i + degree + 1] - knots[i + 1])
        if left_den:
            derivatives[i] += degree / left_den * lower[i]
        if right_den:
            derivatives[i] -= degree / right_den * lower[i + 1]
    return values, derivatives


def insert_once(knots, degree, matrix, value):
    """Apply the standard exact B-spline knot-insertion control map."""
    p = degree
    u = value
    n = len(matrix) - 1
    multiplicity = knots.count(u)
    k = next(i for i in range(len(knots) - 1) if knots[i] <= u < knots[i + 1])
    result = [dict() for _ in range(len(matrix) + 1)]

    def put(row, col, coefficient):
        if coefficient:
            result[row][col] = result[row].get(col, Fraction(0)) + coefficient

    for i in range(0, k - p + 1):
        result[i] = dict(matrix[i])
    for i in range(k - multiplicity + 1, n + 2):
        result[i] = dict(matrix[i - 1])
    for i in range(k - p + 1, k - multiplicity + 1):
        denominator = knots[i + p] - knots[i]
        alpha = (u - knots[i]) / denominator if denominator else Fraction(0)
        for col, coefficient in matrix[i].items():
            put(i, col, alpha * coefficient)
        for col, coefficient in matrix[i - 1].items():
            put(i, col, (1 - alpha) * coefficient)
    return knots[: k + 1] + [u] + knots[k + 1 :], result


def refinement_matrix(old_knots, new_knots, degree=3):
    knots = list(old_knots)
    count = len(knots) - degree - 1
    matrix = [{i: Fraction(1)} for i in range(count)]
    for value in sorted(set(new_knots) - set(old_knots)):
        target_count = new_knots.count(value) - old_knots.count(value)
        for _ in range(target_count):
            knots, matrix = insert_once(knots, degree, matrix, value)
    if knots != list(new_knots):
        raise ValueError("the serialized spline levels are not a pure knot refinement")
    return matrix


def tensor_refinement(axes, old_shape):
    new_shape = [len(axis) for axis in axes]
    result = []
    for i in range(new_shape[0]):
        for j in range(new_shape[1]):
            for k in range(new_shape[2]):
                row = {}
                for a, ca in axes[0][i].items():
                    for b, cb in axes[1][j].items():
                        for c, cc in axes[2][k].items():
                            col = flatten((a, b, c), old_shape)
                            row[col] = ca * cb * cc
                result.append(row)
    return result


def support_inside(index, shape, knots, region, degree=3):
    ijk = unravel(index, shape)
    return all(
        knots[axis][basis] >= region[axis][0]
        and knots[axis][basis + degree + 1] <= region[axis][1]
        for axis, basis in enumerate(ijk)
    )


def tensor_eval(knots, shape, point):
    axes = [basis_values(axis, 3, x) for axis, x in zip(knots, point)]
    values = [0.0] * math.prod(shape)
    grads = [[0.0, 0.0, 0.0] for _ in values]
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                index = flatten((i, j, k), shape)
                values[index] = axes[0][0][i] * axes[1][0][j] * axes[2][0][k]
                grads[index] = [
                    axes[0][1][i] * axes[1][0][j] * axes[2][0][k],
                    axes[0][0][i] * axes[1][1][j] * axes[2][0][k],
                    axes[0][0][i] * axes[1][0][j] * axes[2][1][k],
                ]
    return values, grads


def direct_thb(doc, point, fine_to_active_coarse):
    pknots = [[frac(x) for x in axis] for axis in doc["prepared_knots"]]
    fknots = [[frac(x) for x in axis] for axis in doc["fine_knots"]]
    pshape = [int(x) for x in doc["prepared_shape"]]
    fshape = tensor_shape(fknots)
    pvalues, pgrads = tensor_eval(pknots, pshape, point)
    fvalues, fgrads = tensor_eval(fknots, fshape, point)
    coarse = [int(x) for x in doc["coarse_active_ids"]]
    fine = [int(x) for x in doc["fine_active_ids"]]
    polynomial = [0.0] * len(doc["field_weights"])
    derivative = [[0.0, 0.0, 0.0] for _ in polynomial]
    for dof, basis in enumerate(coarse):
        polynomial[dof] = pvalues[basis]
        derivative[dof] = list(pgrads[basis])
    for pos, basis in enumerate(fine):
        if fvalues[basis] == 0.0 and not any(fgrads[basis]):
            continue
        for coarse_dof, coefficient in fine_to_active_coarse[pos]:
            factor = float(coefficient)
            polynomial[coarse_dof] -= factor * fvalues[basis]
            for axis in range(3):
                derivative[coarse_dof][axis] -= factor * fgrads[basis][axis]
    offset = len(coarse)
    for pos, basis in enumerate(fine):
        polynomial[offset + pos] = fvalues[basis]
        derivative[offset + pos] = list(fgrads[basis])
    weights = [float(frac(x)) for x in doc["field_weights"]]
    weighted = [w * n for w, n in zip(weights, polynomial)]
    denom = sum(weighted)
    ddenom = [sum(weights[i] * derivative[i][axis] for i in range(len(weights))) for axis in range(3)]
    rational = [value / denom for value in weighted]
    rational_grad = [[
        (weights[i] * derivative[i][axis] - rational[i] * ddenom[axis]) / denom
        for axis in range(3)
    ] for i in range(len(weights))]
    return rational, rational_grad


def validate_active_sets(doc):
    pknots = [[frac(x) for x in axis] for axis in doc["prepared_knots"]]
    fknots = [[frac(x) for x in axis] for axis in doc["fine_knots"]]
    pshape = [int(x) for x in doc["prepared_shape"]]
    fshape = tensor_shape(fknots)
    region = [[frac(x) for x in row] for row in doc["interior_region"]]
    coarse = [i for i in range(math.prod(pshape)) if not support_inside(i, pshape, pknots, region)]
    fine = [i for i in range(math.prod(fshape)) if support_inside(i, fshape, fknots, region)]
    if coarse != [int(x) for x in doc["coarse_active_ids"]]:
        raise ValueError("coarse active set disagrees with independent support test")
    if fine != [int(x) for x in doc["fine_active_ids"]]:
        raise ValueError("fine active set disagrees with independent support test")
    return coarse, fine


def validate_leaf_rows(doc):
    global_weights = [frac(x) for x in doc["field_weights"]]
    for leaf in doc["leaves"]:
        dofs = [int(x) for x in leaf["dofs"]]
        weights = [frac(x) for x in leaf["field_weights"]]
        if len(dofs) != len(weights) or any(dof < 0 or dof >= len(global_weights) for dof in dofs):
            raise ValueError("leaf local basis numbering or weight count is invalid")
        if any(weights[i] != global_weights[dof] for i, dof in enumerate(dofs)):
            raise ValueError("leaf field weights disagree with global THB coefficients")
        rows = leaf["extraction"]
        if len(rows) != 64:
            raise ValueError("tricubic leaf extraction must contain 64 Bernstein rows")
        for row in rows:
            total = Fraction(0)
            for column, coefficient in row:
                if int(column) < 0 or int(column) >= len(dofs):
                    raise ValueError("leaf extraction column is outside its local basis")
                value = frac(coefficient)
                if value < 0:
                    raise ValueError("leaf extraction contains a negative coefficient")
                total += value
            if total != 1:
                raise ValueError("leaf extraction row does not preserve constants")


def verify_refusal_controls(doc):
    """Ensure the independent audit rejects corrupted IDs and extraction rows."""
    bad_ids = copy.deepcopy(doc)
    bad_ids["coarse_active_ids"][0] = int(bad_ids["coarse_active_ids"][0]) + 1
    try:
        validate_active_sets(bad_ids)
    except (ValueError, KeyError, IndexError):
        pass
    else:
        raise AssertionError("corrupt support ID was not refused")

    bad_row = copy.deepcopy(doc)
    entry = next(
        row[0]
        for leaf in bad_row["leaves"]
        for row in leaf["extraction"]
        if row
    )
    entry[1] = [str(frac(entry[1]).numerator + frac(entry[1]).denominator), str(frac(entry[1]).denominator)]
    try:
        validate_leaf_rows(bad_row)
    except (ValueError, KeyError, IndexError):
        pass
    else:
        raise AssertionError("corrupt extraction row was not refused")


def extracted_thb(leaf, uvw):
    intervals = [[frac(x) for x in row] for row in leaf["intervals"]]
    extraction = leaf["extraction"]
    weights = [float(frac(x)) for x in leaf["field_weights"]]
    basis = [bernstein(3, t) for t in uvw]
    polynomial = [0.0] * len(weights)
    derivative_local = [[0.0, 0.0, 0.0] for _ in weights]
    for i in range(4):
        for j in range(4):
            for k in range(4):
                row_id = (i * 4 + j) * 4 + k
                b = basis[0][0][i] * basis[1][0][j] * basis[2][0][k]
                db = [
                    basis[0][1][i] * basis[1][0][j] * basis[2][0][k],
                    basis[0][0][i] * basis[1][1][j] * basis[2][0][k],
                    basis[0][0][i] * basis[1][0][j] * basis[2][1][k],
                ]
                for column, exact in extraction[row_id]:
                    column = int(column)
                    coefficient = float(frac(exact))
                    polynomial[column] += coefficient * b
                    for axis in range(3):
                        derivative_local[column][axis] += coefficient * db[axis]
    derivatives = [
        [derivative_local[i][axis] / float(intervals[axis][1] - intervals[axis][0]) for axis in range(3)]
        for i in range(len(weights))
    ]
    weighted = [weights[i] * polynomial[i] for i in range(len(weights))]
    denom = sum(weighted)
    ddenom = [sum(weights[i] * derivatives[i][axis] for i in range(len(weights))) for axis in range(3)]
    rational = [value / denom for value in weighted]
    rational_grad = [[
        (weights[i] * derivatives[i][axis] - rational[i] * ddenom[axis]) / denom
        for axis in range(3)
    ] for i in range(len(weights))]
    return rational, rational_grad


def direct_parent_geometry(model, knots, point):
    lower = model["lower"]
    upper = model["upper"]
    height = float(model["height_mm"])
    axes = [basis_values(knots[0], 3, point[0]), basis_values([0, 0, 1, 1], 1, point[1]), basis_values([0, 0, 1, 1], 1, point[2])]
    shape = [len(lower["controls"]), 2, 2]
    homogeneous = [0.0] * 4
    derivatives = [[0.0] * 4 for _ in range(3)]
    for i in range(shape[0]):
        for j in range(2):
            curve = lower["controls"] if j == 0 else upper["controls"]
            xw, yw, weight = map(float, curve[i])
            for k in range(2):
                basis = axes[0][0][i] * axes[1][0][j] * axes[2][0][k]
                dbasis = [
                    axes[0][1][i] * axes[1][0][j] * axes[2][0][k],
                    axes[0][0][i] * axes[1][1][j] * axes[2][0][k],
                    axes[0][0][i] * axes[1][0][j] * axes[2][1][k],
                ]
                control = [xw, yw, (height * weight if k else 0.0), weight]
                for c in range(4):
                    homogeneous[c] += basis * control[c]
                    for axis in range(3):
                        derivatives[axis][c] += dbasis[axis] * control[c]
    weight = homogeneous[3]
    position = [homogeneous[c] / weight for c in range(3)]
    jacobian = [[
        (derivatives[axis][c] * weight - homogeneous[c] * derivatives[axis][3]) / weight**2
        for c in range(3)
    ] for axis in range(3)]
    return position, jacobian


def local_rational_geometry(leaf, uvw):
    controls = [[float(frac(x)) for x in control] for control in leaf["geometry_controls_homogeneous"]]
    degrees = [3, 1, 1]
    shape = [4, 2, 2]
    axes = [bernstein(degree, uvw[axis]) for axis, degree in enumerate(degrees)]
    h = [0.0] * 4
    dh = [[0.0] * 4 for _ in range(3)]
    for i in range(4):
        for j in range(2):
            for k in range(2):
                control = controls[flatten((i, j, k), shape)]
                basis = axes[0][0][i] * axes[1][0][j] * axes[2][0][k]
                dbasis = [
                    axes[0][1][i] * axes[1][0][j] * axes[2][0][k],
                    axes[0][0][i] * axes[1][1][j] * axes[2][0][k],
                    axes[0][0][i] * axes[1][0][j] * axes[2][1][k],
                ]
                for c in range(4):
                    h[c] += basis * control[c]
                    for axis in range(3):
                        dh[axis][c] += dbasis[axis] * control[c]
    position = [h[c] / h[3] for c in range(3)]
    local_jacobian = [[(dh[a][c] * h[3] - h[c] * dh[a][3]) / h[3]**2 for c in range(3)] for a in range(3)]
    intervals = [[frac(x) for x in row] for row in leaf["intervals"]]
    jacobian = [[local_jacobian[a][c] / float(intervals[a][1] - intervals[a][0]) for c in range(3)] for a in range(3)]
    return position, jacobian


def bernstein(degree, t):
    if degree == 1:
        return [1-t, t], [-1.0, 1.0]
    if degree == 3:
        return (
            [(1-t)**3, 3*t*(1-t)**2, 3*t*t*(1-t), t**3],
            [-3*(1-t)**2, 3*(1-t)*(1-3*t), 3*t*(2-3*t), 3*t*t],
        )
    raise ValueError("unsupported Bernstein degree")


def gauss2():
    # Two-point tensor Gauss is enough for the parametric mass of the unity
    # sum; basis/derivative agreement is checked at a denser targeted grid.
    nodes = [-0.5773502691896257, 0.5773502691896257]
    weights = [1.0, 1.0]
    return [(0.5*(x+1), 0.5*w) for x, w in zip(nodes, weights)]


def check(path, tolerance=1e-10, model_path=None, output_path=None):
    doc = json.loads(Path(path).read_text())
    extraction_sha256 = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if doc.get("schema") != "thb_strip_exact_v1":
        raise ValueError("expected thb_strip_exact_v1 document")
    if doc.get("degrees") != [3, 3, 3] or doc.get("storage_order") != "w-fast, then v, then u":
        raise ValueError("unsupported degree or storage-order convention")
    validate_leaf_rows(doc)
    verify_refusal_controls(doc)
    pknots = [[frac(x) for x in axis] for axis in doc["prepared_knots"]]
    fknots = [[frac(x) for x in axis] for axis in doc["fine_knots"]]
    iknots = [[frac(x) for x in axis] for axis in doc["initial_knots"]]
    region = [[frac(x) for x in row] for row in doc["interior_region"]]
    pshape = [int(x) for x in doc["prepared_shape"]]
    initial_shape = tensor_shape(iknots)
    fshape = tensor_shape(fknots)
    prepared_range = [[int(x) for x in row] for row in doc["refined_cell_range"]]
    refined_spans = [2 * (hi - lo) for lo, hi in prepared_range]
    if math.prod(initial_shape) != 80 or math.prod(pshape) != 539:
        raise ValueError("pilot does not have the preregistered 80-to-539 tensor hierarchy")
    if refined_spans != [8, 4, 4]:
        raise ValueError("pilot refinement region is not the preregistered 8x4x4 interior")
    source_scalar_dofs = None
    source_model = None
    if model_path is not None:
        model_bytes = Path(model_path).read_bytes()
        source_model = json.loads(model_bytes)
        model_sha256 = hashlib.sha256(model_bytes).hexdigest()
        if source_model.get("schema") != "nurbs_sketch_strip_v1":
            raise ValueError("source model is not the admitted NURBS strip fixture")
        lower = source_model.get("lower", {}).get("controls", [])
        upper = source_model.get("upper", {}).get("controls", [])
        if len(lower) != 5 or len(upper) != 5:
            raise ValueError("source strip parent does not have 20 scalar controls")
        source_scalar_dofs = len(lower) * 2 * 2
    expected_coarse, expected_fine = validate_active_sets(doc)
    axes = [refinement_matrix(pknots[a], fknots[a]) for a in range(3)]
    transfer = tensor_refinement(axes, pshape)
    coarse_position = {basis: pos for pos, basis in enumerate(expected_coarse)}
    fine_to_active_coarse = [
        [(coarse_position[old_id], coefficient)
         for old_id, coefficient in transfer[fine_id].items()
         if old_id in coarse_position]
        for fine_id in expected_fine
    ]
    allweights = [float(frac(x)) for x in doc["field_weights"]]
    if len(allweights) != len(expected_coarse) + len(expected_fine):
        raise ValueError("global THB weight count does not match active functions")

    eps = 1e-8
    sample_t = [eps, 0.5, 1-eps]
    max_value_error = 0.0
    max_derivative_error = 0.0
    max_partition_error = 0.0
    max_gradient_sum = 0.0
    max_extracted_sum = 0.0
    max_geometry_position_error = 0.0
    max_geometry_jacobian_error = 0.0
    sample_count = 0
    for leaf in doc["leaves"]:
        intervals = [[frac(x) for x in row] for row in leaf["intervals"]]
        dofs = [int(x) for x in leaf["dofs"]]
        for uvw in itertools.product(sample_t, repeat=3):
            point = [float(intervals[a][0]) + uvw[a] * float(intervals[a][1] - intervals[a][0]) for a in range(3)]
            direct_v, direct_g = direct_thb(doc, point, fine_to_active_coarse)
            extracted_v, extracted_g = extracted_thb(leaf, uvw)
            local_direct_v = [direct_v[dof] for dof in dofs]
            local_direct_g = [direct_g[dof] for dof in dofs]
            max_value_error = max(max_value_error, *(abs(a-b) for a,b in zip(local_direct_v, extracted_v)))
            max_derivative_error = max(max_derivative_error, *(
                abs(local_direct_g[i][axis] - extracted_g[i][axis])
                for i in range(len(dofs)) for axis in range(3)
            ))
            max_partition_error = max(max_partition_error, abs(sum(direct_v)-1.0))
            max_gradient_sum = max(max_gradient_sum, *(abs(sum(g[a] for g in direct_g)) for a in range(3)))
            max_extracted_sum = max(max_extracted_sum, abs(sum(extracted_v)-1.0))
            if source_model is not None:
                parent_position, parent_jacobian = direct_parent_geometry(source_model, iknots, point)
                leaf_position, leaf_jacobian = local_rational_geometry(leaf, uvw)
                max_geometry_position_error = max(
                    max_geometry_position_error,
                    *(abs(a-b) for a,b in zip(parent_position, leaf_position)),
                )
                max_geometry_jacobian_error = max(
                    max_geometry_jacobian_error,
                    *(abs(parent_jacobian[a][c]-leaf_jacobian[a][c]) for a in range(3) for c in range(3)),
                )
            sample_count += 1

    # Integrate the partition of unity over the parametric support. This is a
    # numerical mass check; pointwise equality is checked independently above.
    gauss = gauss2()
    parametric_mass = 0.0
    for leaf in doc["leaves"]:
        intervals = [[frac(x) for x in row] for row in leaf["intervals"]]
        jacobian = math.prod(float(hi-lo) for lo,hi in intervals)
        for x, wx in gauss:
            for y, wy in gauss:
                for z, wz in gauss:
                    point = [float(intervals[0][0])+x*float(intervals[0][1]-intervals[0][0]),
                             float(intervals[1][0])+y*float(intervals[1][1]-intervals[1][0]),
                             float(intervals[2][0])+z*float(intervals[2][1]-intervals[2][0])]
                    values, _ = direct_thb(doc, point, fine_to_active_coarse)
                    parametric_mass += sum(values) * jacobian * wx * wy * wz
    domain_measure = math.prod(float(knots[-4] - knots[3]) for knots in iknots)
    report = {
        "schema": "thb_function_values_oracle_v1",
        "input": Path(path).name,
        "source_model": Path(model_path).name if model_path is not None else None,
        "input_sha256": extraction_sha256,
        "oracle_script": "scripts/check_thb_function_values.py",
        "source_model_sha256": model_sha256 if model_path is not None else None,
        "independent_method": "Cox-de Boor values/derivatives plus exact standard knot insertion; extracted Bernstein rows evaluated separately",
        "active_counts": {"coarse": len(expected_coarse), "fine": len(expected_fine), "total": len(allweights)},
        "tensor_space_counts": {"initial": math.prod(initial_shape), "prepared": math.prod(pshape)},
        "source_parent_scalar_dofs": source_scalar_dofs,
        "refined_fine_span_shape": refined_spans,
        "leaves": len(doc["leaves"]),
        "sample_count": sample_count,
        "sample_local_coordinates": sample_t,
        "corrupt_support_id_refused": True,
        "corrupt_extraction_row_refused": True,
        "max_basis_value_error": max_value_error,
        "max_parametric_derivative_error": max_derivative_error,
        "max_partition_of_unity_error": max_partition_error,
        "max_derivative_sum_error": max_gradient_sum,
        "max_extracted_partition_error": max_extracted_sum,
        "max_source_geometry_position_error_mm": max_geometry_position_error if source_model is not None else None,
        "max_source_geometry_parametric_jacobian_error": max_geometry_jacobian_error if source_model is not None else None,
        "integrated_parametric_partition_mass": parametric_mass,
        "mass_quadrature": "tensor Gauss-2",
        "parametric_domain_measure": domain_measure,
        "parametric_mass_error": abs(parametric_mass-domain_measure),
        "tolerance": tolerance,
        "status": "pass" if max(max_value_error, max_derivative_error, max_partition_error, max_gradient_sum, max_extracted_sum, max_geometry_position_error, max_geometry_jacobian_error, abs(parametric_mass-domain_measure)) <= tolerance else "fail",
    }
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(serialized, end="")
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        Path(output_path).write_text(serialized)
    if report["status"] != "pass":
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("document", type=Path)
    parser.add_argument("--model", type=Path, help="source recipe-built model for the 20-control parent check")
    parser.add_argument("--output", type=Path, help="optional machine-readable report path")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    check(args.document, args.tolerance, args.model, args.output)


if __name__ == "__main__":
    main()
