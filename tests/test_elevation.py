"""Degree elevation checked against independent tensor polynomial evaluation."""

from fractions import Fraction
from itertools import product
from math import comb
import unittest

from trivariate.pierced import pierced_plate_blocks
from trivariate.trace import elevate_trace, extract_trace, traces_match


def evaluate(trace, u, v):
    p, q = trace["degrees"]
    result = [Fraction(0)]*4
    for (i, j), ctrl in zip(product(range(p+1), range(q+1)), trace["controls"], strict=True):
        factor = comb(p,i)*u**i*(1-u)**(p-i) * comb(q,j)*v**j*(1-v)**(q-j)
        for k, value in enumerate(ctrl):
            result[k] += factor*Fraction(value)
    return result


class ElevationTests(unittest.TestCase):
    def test_rational_patch_preserved_on_unisolvent_grid(self):
        original = {"degrees":[1, 2], "controls":[
            [i+j, j*j-i, 3*i+j, Fraction(2+i+2*j, 7)]
            for i,j in product(range(2),range(3))
        ]}
        elevated = elevate_trace(original, (3,4))
        for i,j in product(range(4),range(5)):
            uv = Fraction(i,3),Fraction(j,4)
            self.assertEqual(evaluate(elevated,*uv), evaluate(original,*uv))
        with self.assertRaises(ValueError):
            elevate_trace(original,(1,1))

    def test_east_plate_outer_face_matches_elevated_affine_rectangle(self):
        face = extract_trace(pierced_plate_blocks()[0], "v0")
        affine = {"degrees":[1,1],"controls":[[10,y,z,1] for y in (-10,10) for z in (0,8)]}
        self.assertTrue(traces_match(face, elevate_trace(affine, face["degrees"])))
