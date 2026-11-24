"""Exact tensor-product Bernstein subdivision by de Casteljau."""

from fractions import Fraction
from math import isfinite


def _rational(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction)):
        raise ValueError(f"{label} must be a finite real number")
    if isinstance(value, float) and not isfinite(value):
        raise ValueError(f"{label} must be a finite real number")
    return Fraction(value)


def _axis_index(axis):
    if isinstance(axis, bool):
        raise ValueError("axis must be 0, 1, 2, or u, v, w")
    if isinstance(axis, int) and axis in (0, 1, 2):
        return axis
    if isinstance(axis, str) and axis in "uvw" and len(axis) == 1:
        return "uvw".index(axis)
    raise ValueError("axis must be 0, 1, 2, or u, v, w")


def _de_casteljau(line, cut):
    levels = [line]
    while len(levels[-1]) > 1:
        previous = levels[-1]
        levels.append([
            tuple((1 - cut) * a + cut * b for a, b in zip(previous[i], previous[i + 1]))
            for i in range(len(previous) - 1)
        ])
    degree = len(line) - 1
    left = [levels[index][0] for index in range(degree + 1)]
    right = [levels[degree - index][index] for index in range(degree + 1)]
    return left, right


def _subdivide_tensor(degrees, values, split_axis, cut):
    shape = tuple(degree + 1 for degree in degrees)
    parts = [[], []]
    other_axes = tuple(axis for axis in range(3) if axis != split_axis)
    for first in range(shape[other_axes[0]]):
        for second in range(shape[other_axes[1]]):
            line = []
            for index in range(shape[split_axis]):
                coordinate = [0, 0, 0]
                coordinate[split_axis] = index
                coordinate[other_axes[0]] = first
                coordinate[other_axes[1]] = second
                flat_index = (coordinate[0] * shape[1] + coordinate[1]) * shape[2] + coordinate[2]
                line.append(values[flat_index])
            left, right = _de_casteljau(line, cut)
            for target, split_line in zip(parts, (left, right)):
                for index, value in enumerate(split_line):
                    coordinate = [0, 0, 0]
                    coordinate[split_axis] = index
                    coordinate[other_axes[0]] = first
                    coordinate[other_axes[1]] = second
                    flat_index = (coordinate[0] * shape[1] + coordinate[1]) * shape[2] + coordinate[2]
                    target.append((flat_index, value))
    ordered = []
    for part in parts:
        result = [None] * len(values)
        for flat_index, value in part:
            result[flat_index] = value
        ordered.append(result)
    return ordered


def subdivide(block, axis, cut, *, field=None):
    """Split a homogeneous trivariate Bernstein block at an exact parameter.

    ``axis`` is ``u``, ``v``, ``w`` or its index ``0``, ``1``, ``2``.
    ``cut`` lies in the closed unit interval. If supplied, ``field`` contains
    physical scalar Bernstein coefficients in the same tensor ordering.
    Returned controls and field coefficients are :class:`Fraction` values.
    """
    if not isinstance(block, dict) or not {"degrees", "controls"} <= block.keys():
        raise ValueError("block must contain degrees and controls")
    degrees = block["degrees"]
    if (
        not isinstance(degrees, (list, tuple))
        or len(degrees) != 3
        or any(type(degree) is not int or degree < 0 for degree in degrees)
    ):
        raise ValueError("degrees must be three nonnegative integers")
    degrees = tuple(degrees)
    split_axis = _axis_index(axis)
    cut = _rational(cut, "cut")
    if not 0 < cut < 1:
        raise ValueError("cut must lie strictly inside the unit interval")

    controls = block["controls"]
    try:
        controls = list(controls)
    except TypeError as error:
        raise ValueError("controls must be a flat homogeneous sequence") from error
    expected = (degrees[0] + 1) * (degrees[1] + 1) * (degrees[2] + 1)
    if len(controls) != expected:
        raise ValueError("control count does not match degrees")
    homogeneous = []
    for point in controls:
        if not isinstance(point, (list, tuple)) or len(point) != 4:
            raise ValueError("each control must be a homogeneous quadruple")
        values = tuple(_rational(value, "control coefficient") for value in point)
        if values[3] <= 0:
            raise ValueError("homogeneous weights must be strictly positive")
        homogeneous.append(values)

    geometry_parts = _subdivide_tensor(degrees, homogeneous, split_axis, cut)
    result = [
        {"degrees": list(degrees), "controls": [list(point) for point in part]}
        for part in geometry_parts
    ]
    if field is not None:
        try:
            field_values = list(field)
        except TypeError as error:
            raise ValueError("field must be a flat scalar coefficient sequence") from error
        if len(field_values) != expected:
            raise ValueError("field coefficient count does not match degrees")
        field_values = [
            _rational(value, "field coefficient") * control[3]
            for value, control in zip(field_values, homogeneous, strict=True)
        ]
        field_parts = _subdivide_tensor(
            degrees, [(value,) for value in field_values], split_axis, cut
        )
        for output, numerator_part in zip(result, field_parts, strict=True):
            output["field"] = [
                numerator[0] / control[3]
                for numerator, control in zip(
                    numerator_part, output["controls"], strict=True
                )
            ]
    return tuple(result)
