"""Construct a candidate model from a versioned recipe."""

from itertools import product

from .annular import build_annular_pair
from .bracket import bracket_blocks
from .bracket_audit import design_faces
from .incidence import enumerate_faces


def build(recipe):
    if not isinstance(recipe, dict) or "family" not in recipe:
        raise ValueError("Expected a recipe with family")
    if recipe["family"] == "annular_pair_v1":
        return build_annular_pair(recipe)
    if recipe["family"] == "filleted_bracket_v1":
        if set(recipe) != {"family"}:
            raise ValueError("The fixed filleted bracket recipe only accepts family")
        blocks = bracket_blocks()
        topology = enumerate_faces(blocks)
        labels = design_faces({b["id"]: b for b in blocks}, topology["exterior"])
        return {"schema": recipe["family"], "blocks": blocks,
                "topology": topology, "design_faces": labels}
    if recipe.keys() != {"family", "widths_mm"}:
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
