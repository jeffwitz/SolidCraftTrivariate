"""Independent geometric and volume checks for the fixed perforated bracket."""

import itertools
import json
import math
from pathlib import Path
import unittest

from trivariate.bracket import bracket_blocks


ROOT = Path(__file__).resolve().parents[1]
FIELDCRAFT_REVISION = "35e834652d8548384c2da8ef1b28d588c9af2ffd"
_LAYERS = ((0, 10), (10, 30), (30, 40))


def _basis(degree, t):
    return tuple(
        math.comb(degree, index) * t**index * (1.0 - t) ** (degree - index)
        for index in range(degree + 1)
    )


def _basis_derivative(degree, t):
    if degree == 0:
        return (0.0,)
    lower = _basis(degree - 1, t)
    return tuple(
        degree * ((lower[index - 1] if index else 0.0)
                  - (lower[index] if index < degree else 0.0))
        for index in range(degree + 1)
    )


def _map(block, parameters):
    """Independently evaluate homogeneous Bernstein sums and quotient rule."""
    degrees = block["degrees"]
    basis = [_basis(degree, value) for degree, value in zip(degrees, parameters)]
    dbasis = [_basis_derivative(degree, value)
              for degree, value in zip(degrees, parameters)]
    homogeneous = [0.0] * 4
    derivative = [[0.0] * 4 for _ in range(3)]
    index = 0
    for i, j, k in itertools.product(
        range(degrees[0] + 1), range(degrees[1] + 1), range(degrees[2] + 1)
    ):
        point = block["controls"][index]
        index += 1
        multi = (i, j, k)
        factor = basis[0][i] * basis[1][j] * basis[2][k]
        for component in range(4):
            homogeneous[component] += factor * point[component]
        for axis in range(3):
            term = dbasis[axis][multi[axis]]
            for other in range(3):
                if other != axis:
                    term *= basis[other][multi[other]]
            for component in range(4):
                derivative[axis][component] += term * point[component]
    weight = homogeneous[3]
    position = tuple(value / weight for value in homogeneous[:3])
    columns = tuple(
        tuple((column[c] - position[c] * column[3]) / weight for c in range(3))
        for column in derivative
    )
    return position, columns


def _determinant(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _gauss(order):
    """Compute Gauss-Legendre nodes and weights on [0,1] by Newton steps."""
    nodes, weights = [], []
    for index in range(1, order + 1):
        x = math.cos(math.pi * (index - 0.25) / (order + 0.5))
        for _ in range(32):
            p0, p1 = 1.0, x
            for degree in range(2, order + 1):
                p0, p1 = p1, ((2 * degree - 1) * x * p1
                              - (degree - 1) * p0) / degree
            if order == 1:
                p0, p1 = 1.0, x
            derivative = order * (x * p1 - p0) / (x * x - 1.0)
            step = p1 / derivative
            x -= step
            if abs(step) < 2e-16:
                break
        p0, p1 = 1.0, x
        for degree in range(2, order + 1):
            p0, p1 = p1, ((2 * degree - 1) * x * p1
                          - (degree - 1) * p0) / degree
        if order == 1:
            p0, p1 = 1.0, x
        derivative = order * (x * p1 - p0) / (x * x - 1.0)
        nodes.append((x + 1.0) / 2.0)
        weights.append(1.0 / ((1.0 - x * x) * derivative * derivative))
    return nodes, weights


def _block_volume(block, order, nodes_weights=None):
    nodes, weights = nodes_weights or _gauss(order)
    total = 0.0
    for i, u in enumerate(nodes):
        for j, v in enumerate(nodes):
            for k, w in enumerate(nodes):
                _, columns = _map(block, (u, v, w))
                total += weights[i] * weights[j] * weights[k] * _determinant(columns)
    return total


def _inside_profile_and_holes(point):
    """Closed-form membership oracle for the rounded L profile and two bores."""
    x, y, z = point
    tol = 1e-10
    if not (-tol <= x <= 60.0 + tol and -tol <= y <= 50.0 + tol
            and -tol <= z <= 40.0 + tol):
        return False
    if x <= 8.0 or y <= 8.0:
        in_profile = True
    elif 8.0 <= x <= 14.0:
        boundary = 14.0 - math.sqrt(max(0.0, 36.0 - (x - 14.0) ** 2))
        in_profile = y <= boundary + tol
    else:
        in_profile = False
    if not in_profile:
        return False
    if -tol <= y <= 8.0 + tol and (x - 40.0) ** 2 + (z - 20.0) ** 2 < 16.0 - tol:
        return False
    if -tol <= x <= 8.0 + tol and (y - 32.0) ** 2 + (z - 20.0) ** 2 < 16.0 - tol:
        return False
    return True


class BracketTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blocks = bracket_blocks()

    def test_fixed_36_block_inventory_and_integerized_controls(self):
        self.assertEqual(len(self.blocks), 36)
        identifiers = [block["id"] for block in self.blocks]
        self.assertEqual(len(set(identifiers)), 36)
        self.assertTrue(all(type(value) is int
                            for block in self.blocks
                            for control in block["controls"]
                            for value in control))
        for z0, z1 in _LAYERS:
            self.assertIn(f"elbow_east_z{z0}_{z1}", identifiers)
            self.assertIn(f"arm_x_16_30_z{z0}_{z1}", identifiers)
        for leg in ("x", "y"):
            for sector in ("east", "north", "west", "south"):
                self.assertIn(f"hole_{leg}_sector_{sector}_z10_30", identifiers)

    def test_volume_matches_analytic_formula(self):
        nodes_weights_12 = _gauss(12)
        nodes_weights_16 = _gauss(16)
        volume12 = sum(_block_volume(block, 12, nodes_weights_12)
                       for block in self.blocks)
        volume16 = sum(_block_volume(block, 16, nodes_weights_16)
                       for block in self.blocks)
        expected = 34080.0 - 616.0 * math.pi
        self.assertLessEqual(abs(volume16 - expected), 1e-7)
        self.assertLessEqual(abs(volume16 - volume12), 1e-7)
        self.volume_order12 = volume12
        self.volume_order16 = volume16
        self.volume_error = abs(volume16 - expected)

    def test_volume_matches_locked_fieldcraft_fixture_when_available(self):
        checkout_path = ROOT / "upstream/fieldcraft"
        fixture_path = ROOT / "upstream/fieldcraft/tests/corpus/filleted-bracket-blocks/fixture.json"
        if not checkout_path.is_dir():
            self.skipTest(
                "optional FieldCraft checkout is absent; run `make bootstrap` to install references"
            )
        if not fixture_path.is_file():
            self.fail(f"FieldCraft checkout is present but its locked fixture is missing: {fixture_path}")
        fixture = json.loads(fixture_path.read_text())
        fixture_volume = fixture["values"]["surface.solidcraft_mm3"]["value"]
        lock = (ROOT / "references.lock").read_text()
        self.assertIn(FIELDCRAFT_REVISION, lock)
        nodes_weights_16 = _gauss(16)
        volume16 = sum(_block_volume(block, 16, nodes_weights_16)
                       for block in self.blocks)
        self.assertLessEqual(abs(volume16 - fixture_volume), 1e-7)
        self.fixture_volume = fixture_volume

    def test_hole_and_fillet_radii_and_transverse_hole_centers(self):
        max_fillet_error = 0.0
        max_hole_error = 0.0
        for block in self.blocks:
            identifier = block["id"]
            if identifier.startswith("elbow_"):
                for u in (0.0, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0):
                    point, _ = _map(block, (u, 1.0, 0.5))
                    max_fillet_error = max(
                        max_fillet_error,
                        abs(math.hypot(point[0] - 14.0, point[1] - 14.0) - 6.0),
                    )
            elif identifier.startswith("hole_x_"):
                for u in (0.0, 0.2, 0.5, 0.8, 1.0):
                    point, _ = _map(block, (u, 1.0, 0.5))
                    max_hole_error = max(
                        max_hole_error,
                        abs(math.hypot(point[0] - 40.0, point[2] - 20.0) - 4.0),
                    )
            elif identifier.startswith("hole_y_"):
                for u in (0.0, 0.2, 0.5, 0.8, 1.0):
                    point, _ = _map(block, (u, 1.0, 0.5))
                    max_hole_error = max(
                        max_hole_error,
                        abs(math.hypot(point[1] - 32.0, point[2] - 20.0) - 4.0),
                    )
        self.assertLessEqual(max_fillet_error, 1e-12)
        self.assertLessEqual(max_hole_error, 1e-12)
        self.max_fillet_error = max_fillet_error
        self.max_hole_error = max_hole_error

    def test_discriminating_inside_and_outside_points_and_sampled_blocks(self):
        inside = (
            (30, 4, 5), (4, 30, 35), (10, 9, 20),
            (35, 4, 20), (40, 4, 15), (4, 27, 20),
        )
        outside = (
            (10, 10, 5), (61, 4, 20), (4, 51, 20),
            (40, 4, 20), (4, 32, 20), (43, 4, 20), (4, 35, 20),
        )
        for point in inside:
            self.assertTrue(_inside_profile_and_holes(point), point)
        for point in outside:
            self.assertFalse(_inside_profile_and_holes(point), point)

        parameters = (0.0, 0.5, 1.0)
        for block in self.blocks:
            for u, v, w in itertools.product(parameters, repeat=3):
                point, columns = _map(block, (u, v, w))
                self.assertGreater(_determinant(columns), 0.0,
                                   (block["id"], (u, v, w)))
                self.assertTrue(
                    _inside_profile_and_holes(point),
                    (block["id"], (u, v, w), point),
                )


if __name__ == "__main__":
    unittest.main()
