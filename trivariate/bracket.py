"""Fixed 36-block trivariate construction of a perforated angle bracket."""

from fractions import Fraction
from math import lcm

from .elbow import elbow_blocks
from .pierced import pierced_plate_blocks


_LAYERS = ((0, 10), (10, 30), (30, 40))
_SECTORS = ("east", "north", "west", "south")


def _fraction(value):
    """Convert stored controls to exact rational binary/decimal values."""
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction)):
        raise ValueError("control coefficients must be finite real numbers")
    if isinstance(value, float):
        from math import isfinite

        if not isfinite(value):
            raise ValueError("control coefficients must be finite real numbers")
    return Fraction(value)


def _integerized(block_id, degrees, controls):
    scale = 1
    for point in controls:
        for value in point:
            scale = lcm(scale, value.denominator)
    return {
        "id": block_id,
        "degrees": list(degrees),
        "controls": [[int(value * scale) for value in point] for point in controls],
    }


def _box(block_id, bounds):
    """Create an affine box with u=x, v=y, w=z and k varying fastest."""
    x0, x1, y0, y1, z0, z1 = bounds
    controls = [
        (Fraction(x), Fraction(y), Fraction(z), Fraction(1))
        for x in (x0, x1)
        for y in (y0, y1)
        for z in (z0, z1)
    ]
    return _integerized(block_id, (1, 1, 1), controls)


def _extrude_elbow(block, z0, z1):
    height_scale = Fraction(z1 - z0, 8)
    controls = []
    for point in block["controls"]:
        x, y, z, weight = map(_fraction, point)
        controls.append((x, y, z0 * weight + height_scale * z, weight))
    return _integerized(
        f"elbow_{block['id']}_z{z0}_{z1}", block["degrees"], controls
    )


def _transform_hole(block, leg):
    controls = []
    for point in block["controls"]:
        x, y, z, weight = map(_fraction, point)
        if leg == "x":
            transformed = (40 * weight + x, 8 * weight - z, 20 * weight + y, weight)
        else:
            transformed = (z, 32 * weight + x, 20 * weight + y, weight)
        controls.append(transformed)
    return _integerized(
        f"hole_{leg}_{block['id']}_z10_30",
        block["degrees"],
        controls,
    )


def _layer_boxes(z0, z1, *, middle=False):
    suffix = f"z{z0}_{z1}"
    boxes = [
        _box(f"connector_x_14_16_{suffix}", (14, 16, 0, 8, z0, z1)),
        _box(f"connector_y_14_16_{suffix}", (0, 8, 14, 16, z0, z1)),
    ]
    x_intervals = ((16, 30), (50, 60)) if middle else ((16, 30), (30, 50), (50, 60))
    y_intervals = ((16, 22), (42, 50)) if middle else ((16, 22), (22, 42), (42, 50))
    boxes.extend(
        _box(f"arm_x_{x0}_{x1}_{suffix}", (x0, x1, 0, 8, z0, z1))
        for x0, x1 in x_intervals
    )
    boxes.extend(
        _box(f"arm_y_{y0}_{y1}_{suffix}", (0, 8, y0, y1, z0, z1))
        for y0, y1 in y_intervals
    )
    return boxes


def bracket_blocks():
    """Return 36 blocks for the fixed 60×50×40 mm perforated bracket.

    The layers occupy z ranges ``[0,10]``, ``[10,30]``, and ``[30,40]``.
    Each block's homogeneous controls are stored as integers after applying
    one positive common scale factor to the entire block.
    """
    elbows = elbow_blocks()
    blocks = []
    for z0, z1 in _LAYERS:
        blocks.extend(_extrude_elbow(block, z0, z1) for block in elbows)
        blocks.extend(_layer_boxes(z0, z1, middle=(z0, z1) == (10, 30)))

    plate_sectors = pierced_plate_blocks()
    for sector, plate in zip(_SECTORS, plate_sectors, strict=True):
        if plate["id"] != f"sector_{sector}":
            raise ValueError("plate sector ordering changed")
        blocks.append(_transform_hole(plate, "x"))
    for sector, plate in zip(_SECTORS, plate_sectors, strict=True):
        if plate["id"] != f"sector_{sector}":
            raise ValueError("plate sector ordering changed")
        blocks.append(_transform_hole(plate, "y"))

    if len(blocks) != 36 or len({block["id"] for block in blocks}) != 36:
        raise ValueError("bracket block inventory is inconsistent")
    return blocks
