"""Extract and compare oriented tensor-product Bezier face traces."""

from fractions import Fraction
from math import isfinite


def _number(value):
    """Return a rational representation of a finite JSON-style number."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("homogeneous coefficients must be finite real numbers")
    try:
        if not isfinite(value):
            raise ValueError("homogeneous coefficients must be finite real numbers")
        return Fraction(value)
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValueError("homogeneous coefficients must be finite real numbers") from error


def _validated_tensor(degrees, controls, rank):
    """Validate degrees and homogeneous controls for a tensor of given rank."""
    if (
        not isinstance(degrees, (list, tuple))
        or len(degrees) != rank
        or any(type(degree) is not int or degree < 0 for degree in degrees)
    ):
        raise ValueError(f"degrees must be {rank} nonnegative integers")
    try:
        controls = list(controls)
    except TypeError as error:
        raise ValueError("controls must be a flat coefficient sequence") from error
    expected = 1
    for degree in degrees:
        expected *= degree + 1
    if len(controls) != expected:
        raise ValueError("control count does not match degrees")
    rational_controls = []
    for point in controls:
        if not isinstance(point, (list, tuple)) or len(point) != 4:
            raise ValueError("each control must be a homogeneous quadruple")
        values = tuple(_number(value) for value in point)
        if values[3] <= 0:
            raise ValueError("homogeneous weights must be strictly positive")
        rational_controls.append(tuple(point))
    return tuple(degrees), rational_controls


def _validated_block(block):
    if not isinstance(block, dict) or not {"degrees", "controls"} <= block.keys():
        raise ValueError("block must contain degrees and controls")
    return _validated_tensor(block["degrees"], block["controls"], 3)


def _validated_trace(trace):
    degrees = trace.get("degrees") if isinstance(trace, dict) else None
    controls = trace.get("controls") if isinstance(trace, dict) else None
    return _validated_tensor(degrees, controls, 2)


def extract_trace(block, side):
    """Extract a face trace; its two axes retain parent-axis order."""
    degrees, controls = _validated_block(block)
    if not isinstance(side, str) or len(side) != 2 or side[0] not in "uvw" or side[1] not in "01":
        raise ValueError("side must be one of u0, u1, v0, v1, w0, or w1")
    fixed_axis = "uvw".index(side[0])
    fixed_index = int(side[1]) * degrees[fixed_axis]
    free_axes = tuple(axis for axis in range(3) if axis != fixed_axis)
    trace_degrees = tuple(degrees[axis] for axis in free_axes)
    rows = []
    for i in range(degrees[0] + 1):
        for j in range(degrees[1] + 1):
            for k in range(degrees[2] + 1):
                index = (i * (degrees[1] + 1) + j) * (degrees[2] + 1) + k
                multi = (i, j, k)
                if multi[fixed_axis] == fixed_index:
                    rows.append((multi, controls[index]))
    return {"degrees": list(trace_degrees), "controls": [list(row) for _, row in rows]}


def reorient_trace(trace, *, swap=False, flip_first=False, flip_second=False):
    """Apply a transpose and parameter reversals to a degree-two face trace."""
    for flag in (swap, flip_first, flip_second):
        if type(flag) is not bool:
            raise ValueError("orientation flags must be booleans")
    degrees, controls = _validated_trace(trace)
    if swap:
        degrees = (degrees[1], degrees[0])
        controls = [
            controls[j * (degrees[0] + 1) + i]
            for i in range(degrees[0] + 1)
            for j in range(degrees[1] + 1)
        ]
    first, second = degrees
    grid = [
        [controls[i * (second + 1) + j] for j in range(second + 1)]
        for i in range(first + 1)
    ]
    if flip_first:
        grid.reverse()
    if flip_second:
        grid = [list(reversed(row)) for row in grid]
    return {
        "degrees": list(degrees),
        "controls": [list(point) for row in grid for point in row],
    }


def traces_match(first, second, **orientation):
    """Compare rational trace coefficients up to one positive global scale."""
    first_degrees, first_controls = _validated_trace(first)
    oriented = reorient_trace(second, **orientation)
    second_degrees, second_controls = _validated_trace(oriented)
    if first_degrees != second_degrees:
        return False
    left = [_number(value) for point in first_controls for value in point]
    right = [_number(value) for point in second_controls for value in point]
    scale = left[3] / right[3]
    if scale <= 0:
        return False
    return all(a == scale * b for a, b in zip(left, right))
