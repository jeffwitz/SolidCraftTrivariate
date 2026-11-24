"""Independent geometry and area oracles for the two-block profile elbow."""

import itertools
import math
import unittest

from trivariate.elbow import elbow_blocks
from trivariate.planar_certificate import certify_planar
from trivariate.trace import extract_trace, traces_match


TOL = 1e-12


def _bernstein(degree, t):
    if degree == 1:
        return (1.0 - t, t)
    return ((1.0 - t) ** 2, 2.0 * t * (1.0 - t), t**2)


def _bernstein_derivative(degree, t):
    if degree == 1:
        return (-1.0, 1.0)
    return (-2.0 * (1.0 - t), 2.0 - 4.0 * t, 2.0 * t)


def _map(block, u, v, w):
    """Evaluate the homogeneous tensor map and analytic quotient derivatives."""
    basis = (_bernstein(2, u), _bernstein(1, v), _bernstein(1, w))
    deriv = (_bernstein_derivative(2, u), _bernstein_derivative(1, v),
             _bernstein_derivative(1, w))
    homogeneous = [0.0] * 4
    dh = [[0.0] * 4 for _ in range(3)]
    for i, j, k in itertools.product(range(3), range(2), range(2)):
        control = block["controls"][(i * 2 + j) * 2 + k]
        terms = (basis[0][i] * basis[1][j] * basis[2][k],
                 deriv[0][i] * basis[1][j] * basis[2][k],
                 basis[0][i] * deriv[1][j] * basis[2][k],
                 basis[0][i] * basis[1][j] * deriv[2][k])
        for c in range(4):
            homogeneous[c] += terms[0] * control[c]
            for axis in range(3):
                dh[axis][c] += terms[axis + 1] * control[c]
    weight = homogeneous[3]
    point = tuple(value / weight for value in homogeneous[:3])
    columns = tuple(
        tuple((d[c] - point[c] * d[3]) / weight for c in range(3)) for d in dh
    )
    return point, columns


def _gauss(order):
    """Return Gauss-Legendre nodes and weights on [0,1]."""
    nodes, weights = [], []
    for index in range(1, order + 1):
        x = math.cos(math.pi * (index - 0.25) / (order + 0.5))
        for _ in range(32):
            p0, p1 = 1.0, x
            for degree in range(2, order + 1):
                p0, p1 = p1, ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree
            if order == 1:
                p0, p1 = 1.0, x
            dp = order * (x * p1 - p0) / (x * x - 1.0)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 2e-16:
                break
        p0, p1 = 1.0, x
        for degree in range(2, order + 1):
            p0, p1 = p1, ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree
        if order == 1:
            p0, p1 = 1.0, x
        dp = order * (x * p1 - p0) / (x * x - 1.0)
        nodes.append((x + 1.0) / 2.0)
        weights.append(1.0 / ((1.0 - x * x) * dp * dp))
    return nodes, weights


def _area(block, order):
    nodes, weights = _gauss(order)
    result = 0.0
    for u, wu in zip(nodes, weights):
        for v, wv in zip(nodes, weights):
            _, columns = _map(block, u, v, 0.5)
            result += wu * wv * (columns[0][0] * columns[1][1]
                                 - columns[0][1] * columns[1][0])
    return result


def _expected_arc(u):
    """Closed-form rational arc from -135 to -90 degrees, center (14,14)."""
    q = math.cos(math.pi / 8.0)
    half = math.sqrt(0.5)
    b0, b1, b2 = _bernstein(2, u)
    middle = (14.0 + 6.0 * math.cos(-5.0 * math.pi / 8.0) / q,
              14.0 + 6.0 * math.sin(-5.0 * math.pi / 8.0) / q)
    xh = b0 * (14.0 - 6.0 * half) + b1 * q * middle[0] + b2 * 14.0
    yh = b0 * (14.0 - 6.0 * half) + b1 * q * middle[1] + b2 * 8.0
    wh = b0 + b1 * q + b2
    return xh / wh, yh / wh


class ElbowTests(unittest.TestCase):
    def setUp(self):
        self.blocks = elbow_blocks()

    def test_two_sections_share_exact_diagonal_trace(self):
        self.assertEqual([block["id"] for block in self.blocks], ["east", "north"])
        for block in self.blocks:
            self.assertEqual(block["degrees"], [2, 1, 1])
            self.assertEqual(len(block["controls"]), 12)
        east = extract_trace(self.blocks[0], "u0")
        north = extract_trace(self.blocks[1], "u1")
        self.assertTrue(traces_match(east, north))

    def test_radius_height_orientation_and_planar_certificate(self):
        grid = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)
        max_radius_error = 0.0
        min_det = math.inf
        for block in self.blocks:
            plane_controls = [
                [control[0], control[1], control[3]]
                for control in block["controls"][::2]
            ]
            certificate = certify_planar(block["degrees"][:2], plane_controls)
            self.assertEqual(certificate["status"], "valid", certificate)
            self.assertGreater(certificate["determinant_lower_bound"], 0.0)
            for u, v, w in itertools.product(grid, repeat=3):
                point, columns = _map(block, u, v, w)
                determinant = (
                    columns[0][0] * (columns[1][1] * columns[2][2]
                                     - columns[1][2] * columns[2][1])
                    - columns[0][1] * (columns[1][0] * columns[2][2]
                                       - columns[1][2] * columns[2][0])
                    + columns[0][2] * (columns[1][0] * columns[2][1]
                                       - columns[1][1] * columns[2][0])
                )
                min_det = min(min_det, determinant)
                self.assertGreater(determinant, 0.0)
                self.assertAlmostEqual(point[2], 8.0 * w, delta=TOL)
                if v == 1.0:
                    max_radius_error = max(
                        max_radius_error,
                        abs(math.hypot(point[0] - 14.0, point[1] - 14.0) - 6.0),
                    )
                    if block["id"] == "east":
                        expected = _expected_arc(u)
                        max_radius_error = max(
                            max_radius_error,
                            *(abs(actual - target)
                              for actual, target in zip(point[:2], expected)),
                        )
        self.assertLessEqual(max_radius_error, TOL)
        self.min_determinant = min_det
        self.max_radius_error = max_radius_error

    def test_gauss_area_matches_independent_profile_formula(self):
        area12 = sum(_area(block, 12) for block in self.blocks)
        area16 = sum(_area(block, 16) for block in self.blocks)
        expected = 14.0**2 - 9.0 * math.pi
        self.area = area16
        self.area_error = abs(area16 - expected)
        self.quadrature_delta = abs(area16 - area12)
        self.assertLessEqual(self.quadrature_delta, 1e-10)
        self.assertLessEqual(self.area_error, 1e-10)


if __name__ == "__main__":
    unittest.main()
