"""Independent exact pointwise checks of the continuous certificate bounds."""

from fractions import Fraction as F
from itertools import product
from math import comb, prod
import unittest

from trivariate.audit import audit
from trivariate.build import build


def bernstein(d, i, t):
    if not 0 <= i <= d:
        return 0
    return comb(d, i) * t**i * (1-t)**(d-i)


def homogeneous_position_and_derivatives(block, point):
    """Evaluate the polynomial and its partials using exact binary fractions."""
    ds = block["degrees"]
    result = [[F(0)]*4 for _ in range(4)]
    for idx, ctrl in zip(product(*(range(d+1) for d in ds)), block["controls"], strict=True):
        basis = [bernstein(d, i, t) for d, i, t in zip(ds, idx, point)]
        terms = [prod(basis)]
        for axis, (d, i, t) in enumerate(zip(ds, idx, point)):
            derivative = d*(bernstein(d-1, i-1, t) - bernstein(d-1, i, t))
            terms.append(derivative * basis[(axis+1)%3] * basis[(axis+2)%3])
        for row, term in zip(result, terms):
            for c, val in enumerate(ctrl):
                row[c] += F(val)*term
    value, *partials = result
    xyz = [c/value[3] for c in value[:3]]
    columns = [[(p[c]-xyz[c]*p[3])/value[3] for c in range(3)] for p in partials]
    return xyz, columns, value[3]


def determinant(columns):
    a, b, c = columns
    return sum(a[i]*(b[(i+1)%3]*c[(i+2)%3] - b[(i+2)%3]*c[(i+1)%3]) for i in range(3))


class AnnularCertificateTests(unittest.TestCase):
    def test_report_bounds_cover_extreme_heights_and_rounded_products(self):
        for heights in ([1, 1], [1, 2], [1000, 1], [1000, 546], [1000, 1000]):
            model = build({"family":"annular_pair_v1", "heights_mm":heights})
            report = audit(model)["certificate"]
            self.assertEqual(report["status"], "valid")
            def bound(name):
                return F(**report["bounds"][name])
            radial = bound("radius_error_mm")
            axial = bound("height_error_mm")
            minimum_det = bound("determinant_lower_bound")
            self.assertLessEqual(max(radial, axial), F(1,10**12))
            self.assertGreater(minimum_det, 0)
            start = 0
            for block, height in zip(model["blocks"], heights):
                for u,v,w in ((F(0),F(0),F(0)), (F(1),F(1),F(1)), (F(1,3),F(2,5),F(7,9))):
                    xyz, dx, weight = homogeneous_position_and_derivatives(block, (u,v,w))
                    x,y,z = xyz
                    self.assertGreaterEqual(weight, bound("minimum_weight"))
                    self.assertGreaterEqual(determinant(dx), minimum_det)
                    self.assertLessEqual((2-v-radial)**2, x*x+y*y)
                    self.assertLessEqual(x*x+y*y, (2-v+radial)**2)
                    self.assertLessEqual(abs(z-start-w*height), axial)
                start += height
