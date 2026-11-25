"""Exact synchronized constant-distance offset obstruction on a real plate span."""

import json
from pathlib import Path
import sympy as sp

root = Path(__file__).resolve().parents[2]
recipe = json.loads(
    (root / "tests/references/g2_cross_sketches.recipe.json").read_text()
)
curve = recipe["sketches"][0]["loops"][0]["curve"]
assert curve["degree"] == 3 and curve["knots"] == list(range(-3, 12))
t = sp.symbols("t")
# Six times the four cardinal cubic basis functions on the first active span.
basis = [(1 - t) ** 3, 3 * t**3 - 6 * t**2 + 4, -3 * t**3 + 3 * t**2 + 3 * t + 1, t**3]
X, Y, W = [
    sp.expand(sum(b * p[d] for b, p in zip(basis, curve["controls"][:4])))
    for d in range(3)
]
dx, dy = [sp.diff(h, t) * W - h * sp.diff(W, t) for h in (X, Y)]
speed_numerator = sp.Poly(dx**2 + dy**2, t)
content, primitive = speed_numerator.primitive()
gcd = sp.gcd(primitive, primitive.diff())
assert primitive.degree() == 8 and gcd.degree() == 0
print(
    json.dumps(
        dict(
            span=[0, 1],
            homogeneous=[str(h) for h in (X, Y, W)],
            squared_speed_denominator="W^4",
            squared_speed_numerator_content=int(content),
            primitive_coefficients_descending=[int(c) for c in primitive.all_coeffs()],
            gcd_with_derivative_degree=int(gcd.degree()),
            squarefree_degree=int(primitive.degree()),
            consequence="The constant-distance planar rail has no rational coordinate functions in the source parameter. Knot multiplicity changes do not remove the square root.",
        ),
        indent=2,
    )
)

# Prove that even forgetting the source correspondence does not rationalize
# this rail. Recover the foot parameter from the two offset point coordinates.
# A specialized gcd of degree one bounds the generic degree from above,
# provided both original polynomial degrees are preserved by specialization.
for side in (-1, 1):
    dx0, dy0 = (h.subs(t, 0) for h in (dx, dy))
    length0 = sp.sqrt(dx0**2 + dy0**2)
    u0 = (X / W).subs(t, 0) + side * dy0 / (4 * length0)
    v0 = (Y / W).subs(t, 0) - side * dx0 / (4 * length0)
    a, b = X - u0*W, Y - v0*W
    distance = sp.Poly(a**2 + b**2 - W**2 / 16, t, extension=True)
    normal = sp.Poly(a*dx + b*dy, t, extension=True)
    assert (distance.degree(), normal.degree()) == (6, 7)
    assert sp.gcd(distance, normal).monic().as_expr() == t
