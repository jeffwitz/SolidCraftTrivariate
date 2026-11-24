"""Global section certificates checked on canonical and perturbed arcs."""

from fractions import Fraction as F
import unittest

from trivariate.arc_strip import certify_arc_strip
from trivariate.elbow import elbow_blocks
from trivariate.pierced import pierced_plate_blocks


class ArcStripTests(unittest.TestCase):
    def test_six_native_sections_and_certified_radial_bound(self):
        for blocks,center,r in ((elbow_blocks(),(14,14),6),(pierced_plate_blocks(),(0,0),4)):
            for block in blocks:
                ctrl = [[x,y,q] for x,y,_,q in block["controls"][::2]]
                result = certify_arc_strip(ctrl,center,r)
                self.assertEqual(result["status"],"valid")
                err = result["radius_error_mm"]
                self.assertLessEqual(err,F(1,10**12))
                for u in (F(0),F(2,7),F(1,2),F(1)):
                    basis = ((1-u)**2,2*u*(1-u),u*u)
                    x,y,q = [sum(b*F(p[k]) for b,p in zip(basis,ctrl[1::2])) for k in range(3)]
                    squared = (x/q-center[0])**2+(y/q-center[1])**2
                    self.assertLessEqual((r-err)**2,squared)
                    self.assertLessEqual(squared,(r+err)**2)

    def test_deformed_misplaced_and_reversed_arc(self):
        block = pierced_plate_blocks()[0]
        ctrl = [[x,y,q] for x,y,_,q in block["controls"][::2]]
        self.assertNotEqual(certify_arc_strip(ctrl,(1,0),4)["status"],"valid")
        self.assertNotEqual(certify_arc_strip(ctrl[::-1],(0,0),4)["status"],"valid")
        ctrl[3][0] += 1e-5
        self.assertNotEqual(certify_arc_strip(ctrl,(0,0),4)["status"],"valid")
