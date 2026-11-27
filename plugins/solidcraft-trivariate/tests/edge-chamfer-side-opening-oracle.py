"""Independent cross-section integral for the synthetic side-opening chamfer.

The source is a 22 x 12 x 12 mm box with a radius-2 mm through-bore centered
at (20.5, 4) in XY. Chamfering the vertical edge (22, 0) by equal setback d
removes the triangular prism x - y > 22 - d, with x <= 22 and y >= 0.
The integral subtracts the part of that prism already occupied by the bore.
"""

from __future__ import annotations

import json
import math
from scipy.integrate import quad


RADIUS = 2.0
BORE_X = 20.5
BORE_Y = 4.0
HEIGHT = 12.0
STOCK_X = 22.0
STOCK_Y = 12.0
STOCK_Z = 12.0


def stock_volume() -> float:
    # The bore is clipped by the x=22 stock face. Its circular cap is outside.
    cap_half_chord = STOCK_X - BORE_X
    cap = RADIUS**2 * math.acos(cap_half_chord / RADIUS)
    cap -= cap_half_chord * math.sqrt(RADIUS**2 - cap_half_chord**2)
    disk_in_stock = math.pi * RADIUS**2 - cap
    return STOCK_X * STOCK_Y * STOCK_Z - disk_in_stock * HEIGHT


def bore_area_in_chamfer(d: float) -> float:
    """Integrate the disk portion inside u>=0, y>=0, u+y<=d; u=22-x."""
    u_center = STOCK_X - BORE_X
    u_min, u_max = u_center - RADIUS, u_center + RADIUS

    def vertical_overlap(u: float) -> float:
        h = math.sqrt(max(0.0, RADIUS**2 - (u - u_center) ** 2))
        lower = max(0.0, BORE_Y - h)
        upper = min(d - u, BORE_Y + h, STOCK_Y)
        return max(0.0, upper - lower)

    lo, hi = max(0.0, u_min), min(d, u_max)
    if hi <= lo:
        return 0.0

    # Split at support endpoints and at exact line/circle intersections to
    # keep adaptive quadrature stable when the clipping topology changes.
    points = [lo, hi]
    # Circle intersecting u + y = d, where disk coordinates are (u, y).
    # The line's signed distance to the disk center is
    # (u_center + bore_y - d) / sqrt(2).
    delta = u_center + BORE_Y - d
    if abs(delta) <= RADIUS * math.sqrt(2.0):
        root = math.sqrt(max(0.0, 2.0 * RADIUS**2 - delta**2))
        for u in ((u_center + d - BORE_Y + root) / 2.0,
                  (u_center + d - BORE_Y - root) / 2.0):
            if lo < u < hi:
                points.append(u)
    points.sort()
    return sum(
        quad(vertical_overlap, a, b, epsabs=2e-13, epsrel=2e-13, limit=200)[0]
        for a, b in zip(points, points[1:])
    )


def expected_volume(d: float) -> float:
    wedge_area = d * d / 2.0
    removed = HEIGHT * (wedge_area - bore_area_in_chamfer(d))
    return stock_volume() - removed


if __name__ == "__main__":
    print(json.dumps({
        "method": "independent_xy_cross_section_integral",
        "stock_volume_mm3": stock_volume(),
        "chamfer_distances_mm": {
            str(d): expected_volume(d) for d in (6.0, 6.5, 7.0, 8.0)
        },
    }, indent=2, sort_keys=True))
