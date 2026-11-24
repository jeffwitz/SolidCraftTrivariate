"""Independent coefficient audit for the ``annular_pair_v1`` family."""

from fractions import Fraction
from math import isfinite, sqrt

from .trace import extract_trace, traces_match


_HEIGHT_LIMIT_MM = 1000


def audit(model):
    """Audit two vertically adjacent rational quarter-annulus blocks.

    The oracle reconstructs every homogeneous Bernstein coefficient directly
    from the declared integer heights. It does not call the model builder or
    the spline evaluator. ``valid`` certifies this fixed coefficient family;
    it is not a general-purpose continuous injectivity certificate.
    """

    def invalid(reason):
        return {"status": "invalid", "reason": reason}

    if not isinstance(model, dict) or model.get("schema") != "annular_pair_v1":
        return {"status": "indeterminate", "reason": "unknown_schema"}
    if model.keys() != {"schema", "heights_mm", "blocks", "interface"}:
        return invalid("unexpected_model_fields")
    heights = model["heights_mm"]
    if not isinstance(heights, list) or len(heights) != 2 or any(
        type(height) is not int or not 1 <= height <= _HEIGHT_LIMIT_MM
        for height in heights
    ):
        return invalid("invalid_heights")

    blocks = model["blocks"]
    if not isinstance(blocks, list) or len(blocks) != 2:
        return invalid("block_count")
    z0 = 0
    lower_block = None
    upper_block = None
    for index, (block, block_id, height) in enumerate(
        zip(blocks, ("lower", "upper"), heights)
    ):
        if not isinstance(block, dict) or block.get("id") != block_id:
            return invalid("block_identity")
        if block.keys() != {"id", "degrees", "controls"}:
            return invalid("unexpected_block_fields")
        if (
            not isinstance(block["degrees"], list)
            or len(block["degrees"]) != 3
            or any(type(degree) is not int for degree in block["degrees"])
            or block["degrees"] != [2, 1, 1]
        ):
            return invalid("degree_mismatch")
        controls = block["controls"]
        if not isinstance(controls, list) or len(controls) != 12:
            return invalid("control_count")

        expected = []
        arc_weights = (1.0, sqrt(0.5), 1.0)
        # Index order is angle, radius, z, with z varying fastest.
        for i, (cx, cy) in enumerate(((1, 0), (1, 1), (0, 1))):
            weight = arc_weights[i]
            for j, radius in enumerate((2, 1)):
                for k, local_z in enumerate((z0, z0 + height)):
                    q = Fraction.from_float(weight)
                    x_integer = radius * cx
                    y_integer = radius * cy
                    expected.append([
                        float(q * x_integer), float(q * y_integer),
                        float(q * local_z), weight,
                    ])
        for point, reference in zip(controls, expected):
            if not isinstance(point, list) or len(point) != 4:
                return invalid("control_shape")
            if any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or (isinstance(value, float) and not isfinite(value))
                for value in point
            ):
                return invalid("nonfinite_control")
            if point != reference:
                return invalid("coefficient_mismatch")

        if index == 0:
            lower_block = block
        else:
            upper_block = block
        z0 += height

    try:
        lower_trace = extract_trace(lower_block, "w1")
        upper_trace = extract_trace(upper_block, "w0")
        if not traces_match(lower_trace, upper_trace):
            return invalid("shared_trace_mismatch")
    except (TypeError, ValueError):
        return invalid("shared_trace_invalid")

    if model["interface"] != ["lower.w1", "upper.w0", "identity_uv"]:
        return invalid("interface_incidence")
    return {
        "status": "valid",
        "method": "independent_homogeneous_coefficient_reconstruction",
        "height_total_mm": sum(heights),
        "shared_faces": 1,
        "exterior_sides": 10,
        "coefficient_comparison": "strict_binary64",
    }
