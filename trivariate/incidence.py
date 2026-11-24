"""Whole-face topological incidence for exact trivariate Bezier blocks."""

from fractions import Fraction
from itertools import product
from math import isfinite

from .trace import elevate_trace, extract_trace, reorient_trace, traces_match


_SIDES = ("u0", "u1", "v0", "v1", "w0", "w1")
_ORIENTATIONS = tuple(
    {"swap": swap, "flip_first": flip_first, "flip_second": flip_second}
    for swap, flip_first, flip_second in product((False, True), repeat=3)
)


def _fraction(value):
    if type(value) is int:
        return Fraction(value)
    if type(value) is float and isfinite(value):
        return Fraction.from_float(value)
    if type(value) is Fraction:
        return value
    raise ValueError("homogeneous coefficients must be finite int, float, or Fraction values")


def _validate_blocks(blocks):
    if not isinstance(blocks, (list, tuple)):
        raise ValueError("blocks must be a list or tuple")
    identities = set()
    for block in blocks:
        if not isinstance(block, dict) or not isinstance(block.get("id"), str) or not block["id"]:
            raise ValueError("each block must have a nonempty string id")
        if block["id"] in identities:
            raise ValueError("duplicate block identity")
        identities.add(block["id"])
    return sorted(blocks, key=lambda block: block["id"])


def _face_vertices(trace):
    du, dv = trace["degrees"]
    if du < 1 or dv < 1:
        raise ValueError("degenerate_face_vertices")
    controls = trace["controls"]
    corner_indices = (0, dv, du * (dv + 1), (du + 1) * (dv + 1) - 1)
    vertices = []
    for index in corner_indices:
        wx, wy, wz, weight = (_fraction(value) for value in controls[index])
        if weight <= 0:
            raise ValueError("face weight must be positive")
        vertices.append((wx / weight, wy / weight, wz / weight))
    if len(set(vertices)) != 4:
        raise ValueError("degenerate_face_vertices")
    return frozenset(vertices)


def _oriented_trace_pair(first, second):
    first_degrees = tuple(first["degrees"])
    matches = []
    elevated_first = {}
    for orientation in _ORIENTATIONS:
        oriented = reorient_trace(second, **orientation)
        target_degrees = tuple(max(a, b) for a, b in zip(first_degrees, oriented["degrees"]))
        if target_degrees not in elevated_first:
            elevated_first[target_degrees] = elevate_trace(first, target_degrees)
        elevated_second = elevate_trace(oriented, target_degrees)
        if traces_match(elevated_first[target_degrees], elevated_second):
            matches.append(orientation)
    if not matches:
        raise ValueError("incompatible_face_traces")
    if len(matches) != 1:
        raise ValueError("ambiguous_face_orientation")
    return matches[0]


def enumerate_faces(blocks):
    """Enumerate whole-face interfaces and exterior faces for Bezier blocks.

    Face candidates are grouped by their four exact physical corner vertices.
    A shared corner set is only a candidate: full homogeneous traces must also
    agree after one unique square orientation and exact degree elevation.
    This establishes topology for matching whole faces only. It does not
    establish absence of intersections or partial-face contacts.
    """
    blocks = _validate_blocks(blocks)
    faces = {}
    for block in blocks:
        for side in _SIDES:
            trace = extract_trace(block, side)
            vertices = _face_vertices(trace)
            faces.setdefault(vertices, []).append((block, side, trace))

    interfaces = []
    exterior = []
    for occurrences in faces.values():
        if len(occurrences) > 2:
            raise ValueError("triple_face_incidence")
        if len(occurrences) == 1:
            block, side, _ = occurrences[0]
            exterior.append([block["id"], side])
            continue
        (first_block, first_side, first_trace), (second_block, second_side, second_trace) = occurrences
        orientation = _oriented_trace_pair(first_trace, second_trace)
        interfaces.append({
            "first": [first_block["id"], first_side],
            "second": [second_block["id"], second_side],
            "orientation": orientation,
        })
    return {"interfaces": interfaces, "exterior": exterior}
