"""Reusable rational Bézier blocks for a square plate around a circular hole."""

from math import sqrt


_OUTER = ((10.0, -10.0, 1.0), (10.0, 0.0, 1.0), (10.0, 10.0, 1.0))
_INNER = (
    (2.0 * sqrt(2.0), -2.0 * sqrt(2.0), 1.0),
    (4.0 * sqrt(2.0), 0.0, sqrt(0.5)),
    (2.0 * sqrt(2.0), 2.0 * sqrt(2.0), 1.0),
)
_SECTORS = ("east", "north", "west", "south")


def _quarter_turn(x, y, turns):
    """Rotate a point by exact coordinate permutations in quarter turns."""
    for _ in range(turns):
        x, y = -y, x
    return x, y


def pierced_plate_blocks():
    """Return four degree-(2,1,1) blocks for a 20×20×8 mm pierced plate.

    Each block covers one 90-degree sector, with angle increasing along ``u``,
    the square boundary at ``v=0``, the rational circular boundary at
    ``v=1``, and altitude increasing along ``w``. Coordinates are local to
    the hole center and the plate's bottom plane. Controls are homogeneous
    ``(wx, wy, wz, weight)`` tuples flattened with ``w`` varying fastest.
    """
    blocks = []
    for turns, sector in enumerate(_SECTORS):
        controls = []
        for angle_index in range(3):
            for radial_index in range(2):
                x, y, weight = (
                    _OUTER[angle_index]
                    if radial_index == 0
                    else _INNER[angle_index]
                )
                x, y = _quarter_turn(x, y, turns)
                for z in (0.0, 8.0):
                    controls.append([weight * x, weight * y, weight * z, weight])
        blocks.append(
            {
                "id": f"sector_{sector}",
                "degrees": [2, 1, 1],
                "controls": controls,
            }
        )
    return blocks
