"""Global certificate for a regular planar patch with three straight edges."""

from fractions import Fraction
from math import comb

from .planar_certificate import certify_planar


def _cross(a, b):
    return a[0]*b[1] - a[1]*b[0]


def _difference(a, b):
    return (a[0]-b[0], a[1]-b[1])


def certify_arc_strip(controls, center, radius, tolerance=Fraction(1, 10**12)):
    """Bound arc fidelity and prove injectivity for a degree-(2,1) strip.

    The v=0 side is straight, v=1 is the arc. The other sides are rational
    linear segments. Angles and endpoints are reported in the caller's frame.
    """
    local = certify_planar((2, 1), controls)
    if local["status"] != "valid":
        return local
    center = tuple(map(Fraction, center))
    radius, tolerance = Fraction(radius), Fraction(tolerance)
    if len(center) != 2 or radius <= 0 or tolerance < 0:
        raise ValueError("expected a planar center, positive radius and nonnegative tolerance")

    triples = [tuple(map(Fraction,p)) for p in controls]
    base, arc = triples[::2], triples[1::2]
    def physical(points):
        return [(x/q, y/q) for x,y,q in points]
    b0, bm, b1 = physical(base)
    a0, am, a1 = physical(arc)
    direction = _difference(b1,b0)
    offset = _difference(bm,b0)
    if (_cross(direction,offset) or not 0 < sum(x*y for x,y in zip(direction,offset))
            < sum(x*x for x in direction)):
        return {"status":"indeterminate", "reason":"base_not_proved_monotone_and_straight"}

    # Three sides of the boundary follow the sides of a strictly convex
    # quadrilateral. Strict interior control inequalities exclude crossings.
    corners = (b0,b1,a1,a0)
    edges = [_difference(corners[(i+1)%4], corners[i]) for i in range(4)]
    if any(_cross(edges[i], edges[(i+1)%4]) <= 0 for i in range(4)):
        return {"status":"indeterminate", "reason":"boundary_corner_order"}
    for i in (0,1,3):
        if _cross(edges[i], _difference(am,corners[i])) <= 0:
            return {"status":"indeterminate", "reason":"arc_side_separation"}

    centered = [(x-center[0]*q, y-center[1]*q, q) for x,y,q in arc]
    turns = [_cross(centered[i][:2], centered[j][:2]) for i,j in ((0,1),(0,2),(1,2))]
    if not (all(t>0 for t in turns) or all(t<0 for t in turns)):
        return {"status":"indeterminate", "reason":"arc_polar_monotonicity"}

    # Bernstein product of the centered quadratic homogeneous arc with itself.
    residual = [Fraction(0)]*5
    for i,(x,y,q) in enumerate(centered):
        for j,(xx,yy,qq) in enumerate(centered):
            weight = Fraction(comb(2,i)*comb(2,j), comb(4,i+j))
            residual[i+j] += weight*(x*xx + y*yy - radius**2*q*qq)
    error = max(map(abs,residual))/(min(p[2] for p in centered)**2*radius)
    if error > tolerance:
        return {"status":"indeterminate", "reason":"radial_fidelity_bound_above_budget",
                "radius_error_mm":error}

    return {"status":"valid", "determinant_lower_bound":local["determinant_lower_bound"],
            "radius_error_mm":error, "arc_endpoints":(a0,a1),
            "injectivity":"positive_jacobian_and_simple_oriented_boundary"}
