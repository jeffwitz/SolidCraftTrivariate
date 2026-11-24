"""Construct the two-prism seed in homogeneous Bernstein form."""

from itertools import product


def build(recipe):
    if not isinstance(recipe, dict) or recipe.keys() != {"family", "widths_mm"}:
        raise ValueError("Expected family and widths_mm")
    if recipe["family"] != "affine_pair_v1":
        raise ValueError("Unsupported family")
    widths = recipe["widths_mm"]
    if not isinstance(widths, list) or len(widths) != 2 or any(
        type(x) is not int or x <= 0 for x in widths
    ):
        raise ValueError("Expected two positive integer widths in mm")
    blocks = []
    start = 0
    for name, width in zip(("left", "right"), widths):
        points = [[start + i * width, j, k, 1] for i, j, k in product(range(2), repeat=3)]
        blocks.append({"id": name, "controls": points})
        start += width
    return {"schema": "affine_pair_v1", "widths_mm": widths.copy(), "blocks": blocks,
            "interface": ["left.u1", "right.u0", "identity_vw"]}
