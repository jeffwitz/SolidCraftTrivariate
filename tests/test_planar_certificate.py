"""Independent oracles for the Bernstein determinant certificate."""

from fractions import Fraction as F
from itertools import product
import unittest

from trivariate.pierced import pierced_plate_blocks
from trivariate.planar_certificate import certify_planar
from tests.test_annular_certificate import homogeneous_position_and_derivatives


class PlanarCertificateTests(unittest.TestCase):
    def test_exact_affine_area_scale_and_homogeneous_rescaling(self):
        points = [[2*u+v, -u+3*v, 1] for u,v in product(range(2),repeat=2)]
        for factor in (F(1), F(17,13)):
            result = certify_planar((1,1), [[factor*c for c in pt] for pt in points])
            self.assertEqual(result["status"], "valid")
            self.assertEqual(result["determinant_lower_bound"], 7)

    def test_negative_corner_and_inconclusive_envelope(self):
        folded = [[0,0,1], [0,1,1], [1,0,1], [1,-1,1]]
        self.assertEqual(certify_planar((1,1), folded)["status"], "invalid")
        # y'(v)=(v-1/2)**2+1/10000 is positive but its control polygon crosses zero.
        eps = F(1,10000)
        y = [F(0), (F(1,4)+eps)/3, 2*eps/3, F(1,12)+eps]
        points = [[x,t,1] for x,t in product(range(2),y)]
        self.assertEqual(certify_planar((1,3),points)["status"], "indeterminate")

    def test_all_sectors_against_exact_quotient_rule_at_closed_grid(self):
        for block in pierced_plate_blocks():
            controls = [[x,y,q] for x,y,_,q in block["controls"][::2]]
            result = certify_planar(block["degrees"][:2], controls)
            self.assertEqual(result["status"], "valid")
            lower = result["determinant_lower_bound"]
            self.assertGreater(lower, 0)
            for u,v in product((F(0),F(1,3),F(1)),repeat=2):
                _, (du,dv,_), _ = homogeneous_position_and_derivatives(block,(u,v,F(0)))
                self.assertGreaterEqual(du[0]*dv[1]-du[1]*dv[0], lower)

    def test_malformed_inputs(self):
        for degrees,points in (([1,True], [[1,1,1]]*4), ([1,1], [[1,1,0]]*4),
                               ([1,1], [[1,1,1]]*3), ([1,1], [[1,float("nan"),1]]*4)):
            self.assertEqual(certify_planar(degrees, points)["status"], "invalid")
