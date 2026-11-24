"""Two rational Bézier blocks for the inside corner of a filleted profile."""

from math import cos, pi, sin, sqrt


_CENTER = 14.0
_RADIUS = 6.0
_HEIGHT = 8.0


def _east_controls():
    """Construct the east block from its base, circular arc, and extrusion."""
    quarter_mid = -112.5 * pi / 180.0
    arc_weight = cos(pi / 8.0)
    arc_middle = (
        _CENTER + _RADIUS * cos(quarter_mid) / arc_weight,
        _CENTER + _RADIUS * sin(quarter_mid) / arc_weight,
    )
    arc = (
        (_CENTER - _RADIUS / sqrt(2.0), _CENTER - _RADIUS / sqrt(2.0), 1.0),
        (arc_middle[0], arc_middle[1], arc_weight),
        (_CENTER, _CENTER - _RADIUS, 1.0),
    )
    base = ((0.0, 0.0, 1.0), (7.0, 0.0, 1.0), (14.0, 0.0, 1.0))
    controls = []
    for i in range(3):
        for x, y, weight in (base[i], arc[i]):
            for z in (0.0, _HEIGHT):
                controls.append([weight * x, weight * y, weight * z, weight])
    return controls


def elbow_blocks():
    """Return east and north extrusion blocks around a radius-6 profile bend.

    Both local blocks have degrees ``[2, 1, 1]`` and height 8 mm. The east
    block's ``u`` direction runs left to right, and ``v`` runs from its base
    to the circular arc. The north block is the coordinate-swapped east block
    with ``u`` reversed, giving a positive orientation and an identical shared
    diagonal trace.
    """
    east = _east_controls()
    north = []
    for i in range(3):
        for j in range(2):
            for k in range(2):
                source = east[((2 - i) * 2 + j) * 2 + k]
                x, y, z, weight = source
                north.append([y, x, z, weight])
    return [
        {"id": "east", "degrees": [2, 1, 1], "controls": east},
        {"id": "north", "degrees": [2, 1, 1], "controls": north},
    ]
