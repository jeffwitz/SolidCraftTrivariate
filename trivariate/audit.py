"""Independent exact oracle for the restricted integer affine family."""


def audit(model):
    def invalid(reason):
        return {"status": "invalid", "reason": reason}

    if not isinstance(model, dict) or model.get("schema") != "affine_pair_v1":
        return {"status": "indeterminate", "reason": "unknown_schema"}
    if model.keys() != {"schema", "widths_mm", "blocks", "interface"}:
        return invalid("unexpected_model_fields")
    widths = model.get("widths_mm")
    if not isinstance(widths, list) or len(widths) != 2 or any(
        type(value) is not int or value <= 0 for value in widths
    ):
        return invalid("invalid_dimensions")
    blocks = model.get("blocks")
    if not isinstance(blocks, list) or len(blocks) != 2:
        return invalid("block_count")
    for index, block in enumerate(blocks):
        if not isinstance(block, dict) or block.get("id") != ("left", "right")[index]:
            return invalid("block_identity")
        if block.keys() != {"id", "controls"}:
            return invalid("unexpected_block_fields")
        x0, x1 = sum(widths[:index]), sum(widths[:index + 1])
        controls = block.get("controls")
        if not isinstance(controls, list) or len(controls) != 8:
            return invalid("control_count")
        for n, point in enumerate(controls):
            expected = [x0 if n < 4 else x1, (n // 2) % 2, n % 2, 1]
            if not isinstance(point, list) or any(type(c) is not int for c in point) or point != expected:
                return invalid("affine_coefficient_mismatch")
    if model.get("interface") != ["left.u1", "right.u0", "identity_vw"]:
        return invalid("interface_incidence")
    return {"status": "valid", "method": "exact_integer_affine_coefficients",
            "min_det_mm3": min(widths), "volume_mm3": sum(widths),
            "shared_faces": 1, "exterior_sides": 10}
