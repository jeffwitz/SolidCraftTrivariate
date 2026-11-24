"""Check all denominator derivatives on an asymmetric quadratic patch."""

from fractions import Fraction
from itertools import product
import unittest

from trivariate.spline import RationalBezierMap


class AsymmetricWeightTest(unittest.TestCase):
    def test_center_derivatives(self):
        # At the center the quadratic Bernstein basis is (1,2,1)/4
        # and its derivative is (-1,0,1). Both linear bases are (1,1)/2.
        controls = []
        homogeneous = [0] * 4
        derivatives = [[0] * 4 for _ in range(3)]
        for i, j, k in product(range(3), range(2), range(2)):
            weight = 1 + i + 2 * j + 3 * k
            point = [weight * i, weight * (j + k), weight * (i + k), weight]
            controls.append(point)
            middle = (1, 2, 1)[i]
            # Integer numerators with a common denominator of 16.
            factors = (4 * (i - 1), 2 * middle * (2*j - 1), 2 * middle * (2*k - 1))
            for component, value in enumerate(point):
                homogeneous[component] += middle * value
                for axis, factor in enumerate(factors):
                    derivatives[axis][component] += factor * value

        position, columns = RationalBezierMap((2, 1, 1), controls).evaluate(0.5, 0.5, 0.5)
        for component in range(3):
            expected = Fraction(homogeneous[component], homogeneous[3])
            self.assertAlmostEqual(position[component], float(expected), delta=1e-12)
            for axis, derivative in enumerate(derivatives):
                slope = (derivative[component] - expected * derivative[3]) / homogeneous[3]
                self.assertAlmostEqual(columns[axis][component], float(slope), delta=1e-12)
