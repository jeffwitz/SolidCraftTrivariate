"""Independent geometric checks for the four-sector pierced plate template."""

import itertools
import math
import unittest

from trivariate.pierced import pierced_plate_blocks


TOL = 1e-12
_SECTORS = ("east", "north", "west", "south")


def _basis2(t):
    return ((1.0 - t) ** 2, 2.0 * t * (1.0 - t), t**2)


def _basis2_derivative(t):
    return (-2.0 * (1.0 - t), 2.0 - 4.0 * t, 2.0 * t)


def _evaluate(block, u, v, w):
    """Evaluate homogeneous Bernstein coefficients and quotient derivatives."""
    bu, bv, bw = _basis2(u), (1.0 - v, v), (1.0 - w, w)
    dbu = _basis2_derivative(u)
    dbv, dbw = (-1.0, 1.0), (-1.0, 1.0)
    homogeneous = [0.0] * 4
    derivatives = [[0.0] * 4 for _ in range(3)]
    for i, j, k in itertools.product(range(3), range(2), range(2)):
        index = (i * 2 + j) * 2 + k
        point = block["controls"][index]
        factors = (bu[i] * bv[j] * bw[k], dbu[i] * bv[j] * bw[k],
                   bu[i] * dbv[j] * bw[k], bu[i] * bv[j] * dbw[k])
        for component in range(4):
            homogeneous[component] += factors[0] * point[component]
            for axis in range(3):
                derivatives[axis][component] += factors[axis + 1] * point[component]
    weight = homogeneous[3]
    position = tuple(value / weight for value in homogeneous[:3])
    columns = tuple(
        tuple((d[component] - position[component] * d[3]) / weight
              for component in range(3))
        for d in derivatives
    )
    return position, columns


def _determinant(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _quarter_turn(point, turns):
    x, y = point
    for _ in range(turns):
        x, y = -y, x
    return x, y


def _expected_inner(u, turns):
    """Closed-form rational circular arc, independent of control evaluation."""
    q = math.sqrt(0.5)
    b0, b1, b2 = _basis2(u)
    denominator = b0 + q * b1 + b2
    x = 4.0 * (q * b0 + b1 + q * b2) / denominator
    y = 4.0 * (-q * b0 + q * b2) / denominator
    return _quarter_turn((x, y), turns)


def _gauss_legendre(order):
    """Compute Gauss-Legendre nodes and weights using Newton iteration."""
    nodes, weights = [], []
    for index in range(1, order + 1):
        root = math.cos(math.pi * (index - 0.25) / (order + 0.5))
        for _ in range(32):
            previous, current = 1.0, root
            for degree in range(2, order + 1):
                previous, current = current, (
                    (2 * degree - 1) * root * current - (degree - 1) * previous
                ) / degree
            if order == 1:
                previous, current = 1.0, root
            derivative = order * (root * current - previous) / (root * root - 1.0)
            update = current / derivative
            root -= update
            if abs(update) < 2e-16:
                break
        previous, current = 1.0, root
        for degree in range(2, order + 1):
            previous, current = current, (
                (2 * degree - 1) * root * current - (degree - 1) * previous
            ) / degree
        if order == 1:
            previous, current = 1.0, root
        derivative = order * (root * current - previous) / (root * root - 1.0)
        nodes.append((root + 1.0) / 2.0)
        weights.append(1.0 / ((1.0 - root * root) * derivative * derivative))
    return nodes, weights


def _integrated_volume(blocks, order=16):
    nodes, weights = _gauss_legendre(order)
    total = 0.0
    for block in blocks:
        subtotal = 0.0
        for u, u_weight in zip(nodes, weights):
            for v, v_weight in zip(nodes, weights):
                _, columns = _evaluate(block, u, v, 0.0)
                subtotal += u_weight * v_weight * _determinant(columns)
        total += subtotal
    return total


class PiercedPlateTests(unittest.TestCase):
    def setUp(self):
        self.blocks = pierced_plate_blocks()

    def test_four_ordered_blocks_and_exact_shared_radial_edges(self):
        self.assertEqual([block["id"] for block in self.blocks],
                         [f"sector_{name}" for name in _SECTORS])
        for block in self.blocks:
            self.assertEqual(block["degrees"], [2, 1, 1])
            self.assertEqual(len(block["controls"]), 12)
        for index in range(4):
            current, following = self.blocks[index], self.blocks[(index + 1) % 4]
            current_edge = current["controls"][8:12]
            following_edge = following["controls"][0:4]
            self.assertEqual(current_edge, following_edge)

    def test_radius_altitude_square_boundary_and_positive_orientation(self):
        grid = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)
        min_determinant = math.inf
        max_radius_error = 0.0
        for turns, block in enumerate(self.blocks):
            for u, v, w in itertools.product(grid, repeat=3):
                point, columns = _evaluate(block, u, v, w)
                determinant = _determinant(columns)
                min_determinant = min(min_determinant, determinant)
                self.assertGreater(determinant, 0.0)
                self.assertAlmostEqual(point[2], 8.0 * w, delta=TOL)

                if v == 1.0:
                    expected = _expected_inner(u, turns)
                    max_radius_error = max(
                        max_radius_error,
                        abs(math.hypot(point[0], point[1]) - 4.0),
                        *(abs(actual - target)
                          for actual, target in zip(point[:2], expected)),
                    )
                if v == 0.0:
                    x, y = point[:2]
                    self.assertLessEqual(max(abs(x), abs(y)), 10.0 + TOL)
                    self.assertAlmostEqual(max(abs(x), abs(y)), 10.0, delta=TOL)
                    if turns % 2 == 0:
                        self.assertAlmostEqual(x, 10.0 if turns == 0 else -10.0, delta=TOL)
                    else:
                        self.assertAlmostEqual(y, 10.0 if turns == 1 else -10.0, delta=TOL)
        self.assertLessEqual(max_radius_error, TOL)
        self.min_determinant = min_determinant
        self.max_radius_error = max_radius_error

    def test_volume_by_independent_gauss_legendre_quadrature(self):
        measured_12 = _integrated_volume(self.blocks, order=12)
        measured_16 = _integrated_volume(self.blocks, order=16)
        expected = 8.0 * (400.0 - 16.0 * math.pi)
        error = abs(measured_16 - expected)
        self.measured_volume = measured_16
        self.volume_error = error
        self.quadrature_delta = abs(measured_16 - measured_12)
        self.assertLessEqual(self.quadrature_delta, 1e-8)
        self.assertLessEqual(error, 1e-8)


if __name__ == "__main__":
    unittest.main()
