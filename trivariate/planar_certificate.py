"""Exact Bernstein certificates for local regularity of rational planar maps."""

from fractions import Fraction
from math import comb, isfinite


def _as_fraction(value):
    """Convert an explicitly supported scalar to its exact rational value."""
    if type(value) is int:
        return Fraction(value)
    if type(value) is float and isfinite(value):
        return Fraction.from_float(value)
    if type(value) is Fraction:
        return value
    raise ValueError("controls must contain finite int, float, or Fraction values")


def _shape(degrees, controls):
    if (not isinstance(degrees, (list, tuple)) or len(degrees) != 2
            or any(type(d) is not int or d < 0 for d in degrees)):
        raise ValueError("degrees must be two nonnegative integers")
    du, dv = degrees
    if not isinstance(controls, (list, tuple)) or len(controls) != (du + 1) * (dv + 1):
        raise ValueError("control count does not match the rectangular degree grid")
    grid = []
    for i in range(du + 1):
        row = []
        for j in range(dv + 1):
            point = controls[i * (dv + 1) + j]
            if not isinstance(point, (list, tuple)) or len(point) != 3:
                raise ValueError("each control must be a homogeneous triple")
            values = tuple(_as_fraction(value) for value in point)
            if values[2] <= 0:
                raise ValueError("all Bernstein weight coefficients must be positive")
            row.append(values)
        grid.append(row)
    return (du, dv), grid


def _constant(value, degrees):
    du, dv = degrees
    return [[value for _ in range(dv + 1)] for _ in range(du + 1)]


def _derivative(values, degrees, axis):
    du, dv = degrees
    if axis == 0:
        if du == 0:
            return (0, dv), _constant(Fraction(0), (0, dv))
        return (du - 1, dv), [
            [du * (values[i + 1][j] - values[i][j]) for j in range(dv + 1)]
            for i in range(du)
        ]
    if dv == 0:
        return (du, 0), _constant(Fraction(0), (du, 0))
    return (du, dv - 1), [
        [dv * (values[i][j + 1] - values[i][j]) for j in range(dv)]
        for i in range(du + 1)
    ]


def _multiply(left, left_degrees, right, right_degrees):
    """Multiply tensor-product Bernstein polynomials exactly."""
    lm, ln = left_degrees
    rm, rn = right_degrees
    result = _constant(Fraction(0), (lm + rm, ln + rn))
    for i in range(lm + 1):
        for j in range(ln + 1):
            for k in range(rm + 1):
                for ell in range(rn + 1):
                    factor = Fraction(comb(lm, i) * comb(rm, k), comb(lm + rm, i + k))
                    factor *= Fraction(comb(ln, j) * comb(rn, ell), comb(ln + rn, j + ell))
                    result[i + k][j + ell] += left[i][j] * right[k][ell] * factor
    return (lm + rm, ln + rn), result


def _add(left, right):
    degrees, a = left
    other_degrees, b = right
    if degrees != other_degrees:
        raise AssertionError("determinant terms must have matching degrees")
    return degrees, [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def _product3(a, da, b, db, c, dc):
    ab_degrees, ab = _multiply(a, da, b, db)
    return _multiply(ab, ab_degrees, c, dc)


def _determinant_numerator(degrees, grid):
    """Return Bernstein coefficients of det(H_u,H_v,H)."""
    du, dv = degrees
    components = [
        [[grid[i][j][component] for j in range(dv + 1)] for i in range(du + 1)]
        for component in range(3)
    ]
    if du == 0 or dv == 0:
        return (0, 0), [[Fraction(0)]]
    hu = [_derivative(component, degrees, 0) for component in components]
    hv = [_derivative(component, degrees, 1) for component in components]
    # det(Hu,Hv,H) = Hu_x(Hv_y W-Hv_w Y) - ... ; each triple
    # product has bidegree (3du-1,3dv-1).
    terms = []
    for positive, negative in ((0, 1), (1, 2), (2, 0)):
        first = _product3(hu[positive][1], hu[positive][0],
                          hv[negative][1], hv[negative][0],
                          components[(positive + 2) % 3], degrees)
        second = _product3(hu[positive][1], hu[positive][0],
                           hv[(negative + 1) % 3][1], hv[(negative + 1) % 3][0],
                           components[(positive + 1) % 3], degrees)
        # Cyclic expansion: x(y' w - w' y) + y(w' x - x' w) + w(x' y - y' x).
        terms.append(_add(first, (second[0], [[-x for x in row] for row in second[1]])))
    total = terms[0]
    for term in terms[1:]:
        total = _add(total, term)
    return total


def certify_planar(degrees, controls):
    """Certify strict local orientation of a rational Bézier surface patch.

    Controls are homogeneous ``(wx, wy, weight)`` triples, flattened with the
    second parameter index varying fastest. A valid result proves only local
    regularity and positive orientation throughout the closed parameter square.
    """
    try:
        degrees, grid = _shape(degrees, controls)
    except (ValueError, TypeError, OverflowError) as error:
        return {"status": "invalid", "reason": str(error)}

    numerator_degrees, coefficients = _determinant_numerator(degrees, grid)
    du, dv = degrees
    if du == 0 or dv == 0:
        return {"status": "invalid", "reason": "zero_parametric_derivative"}

    # Corner values are the corresponding Bernstein corner coefficients.
    corners = (coefficients[0][0], coefficients[0][-1],
               coefficients[-1][0], coefficients[-1][-1])
    if any(value <= 0 for value in corners):
        return {"status": "invalid", "reason": "nonpositive_corner_determinant"}

    numerator_lower = min(value for row in coefficients for value in row)
    if numerator_lower <= 0:
        return {"status": "indeterminate", "reason": "bernstein_envelope_not_strictly_positive"}

    weight_upper = max(grid[i][j][2] for i in range(du + 1) for j in range(dv + 1))
    return {
        "status": "valid",
        "determinant_lower_bound": numerator_lower / (weight_upper ** 3),
        "method": "exact_tensor_bernstein_convex_hull",
        "domain": "closed unit square",
        "claim": "strict local regularity and positive orientation only",
        "numerator_degrees": list(numerator_degrees),
    }
