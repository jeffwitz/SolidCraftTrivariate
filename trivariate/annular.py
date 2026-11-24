"""Construction helpers for a pair of vertically adjacent annular patches."""

from itertools import product
from math import sqrt


def _heights(recipe):
    """Validate the versioned height parameters for an annular pair recipe."""
    if not isinstance(recipe, dict) or recipe.keys() != {"family", "heights_mm"}:
        raise ValueError("Expected family and heights_mm")
    heights = recipe["heights_mm"]
    if (
        not isinstance(heights, list)
        or len(heights) != 2
        or any(type(height) is not int or not 0 < height <= 1000 for height in heights)
    ):
        raise ValueError("Expected two positive integer heights in mm, each at most 1000")
    return heights


def _controls(z0, height):
    """Return controls ordered angle, radius, z, with z varying fastest."""
    weights = (1.0, sqrt(0.5), 1.0)
    points = []
    for angle_index, radial_index, z_index in product(range(3), range(2), range(2)):
        radius = (2.0, 1.0)[radial_index]
        weight = weights[angle_index]
        if angle_index == 0:
            x, y = radius, 0.0
        elif angle_index == 1:
            x, y = radius, radius
        else:
            x, y = 0.0, radius
        z = z0 + height * z_index
        points.append([weight * x, weight * y, weight * z, weight])
    return points


def build_annular_pair(recipe):
    """Build two degree-(2,1,1) quarter-annulus patches sharing a z face."""
    heights = _heights(recipe)
    lower, upper = heights
    return {
        "schema": "annular_pair_v1",
        "heights_mm": heights.copy(),
        "blocks": [
            {"id": "lower", "degrees": [2, 1, 1], "controls": _controls(0, lower)},
            {
                "id": "upper",
                "degrees": [2, 1, 1],
                "controls": _controls(lower, upper),
            },
        ],
        "interface": ["lower.w1", "upper.w0", "identity_uv"],
    }
