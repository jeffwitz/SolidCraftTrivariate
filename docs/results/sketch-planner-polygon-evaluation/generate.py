#!/usr/bin/env python3
"""Create the frozen, public synthetic polygonal evaluation corpus."""
import hashlib
import json
import math
import random
from pathlib import Path

OUT = Path(__file__).with_name("cases.json")
CASES = []


def add(name, family, height, outer, cuts, shift=(0, 0), notes=None):
    recipe = {
        "family": "sketch_layer_partition_v1",
        "height_mm": height,
        "outer": {"id": "stock", "vertices": outer},
        "cuts": [
            {"id": f"cut_{i+1}", "vertices": vertices, "z_interval_mm": interval}
            for i, (vertices, interval) in enumerate(cuts)
        ],
        "construction": [],
    }
    if shift != (0, 0):
        recipe["top_shift_mm"] = list(shift)
    CASES.append({"id": name, "family": family, "admissible": True,
                  "recipe": recipe, "design_notes": notes or ""})


# Convex polygonal through holes: vary hole vertex count, location, clearance,
# aspect ratio and count. Every cutter spans the entire extrusion.
RECT = [[-32, -22], [32, -22], [32, 22], [-32, 22]]
add("through_diamond_center", "convex_through_holes", 12, RECT,
    [([[0, -9], [11, 0], [0, 9], [-11, 0]], [0, 12])])
add("through_hex_offset", "convex_through_holes", 15, RECT,
    [([[7, -8], [15, -5], [18, 2], [11, 9], [3, 7], [0, 0]], [0, 15])],
    shift=(2, -1))
add("through_pentagon_tight_clearance", "convex_through_holes", 9,
    [[-18, -12], [18, -12], [18, 12], [-18, 12]],
    [([[0, -8], [9, -2], [6, 7], [-5, 9], [-10, 1]], [0, 9])])
add("through_two_rectangles", "convex_through_holes", 18,
    [[-40, -20], [40, -20], [40, 20], [-40, 20]],
    [([[-28, -6], [-19, -6], [-19, 6], [-28, 6]], [0, 18]),
     ([[13, -9], [27, -9], [27, 9], [13, 9]], [0, 18])])
add("through_octagon_high_aspect", "convex_through_holes", 10,
    [[-50, -14], [50, -14], [50, 14], [-50, 14]],
    [([[-5, -8], [5, -8], [9, -4], [9, 4], [5, 8], [-5, 8], [-9, 4], [-9, -4]],
      [0, 10])])

# Blind pockets: floor parity, profile shape, depth and multiple event layers.
add("blind_rotated_quad", "convex_blind_pockets", 20, RECT,
    [([[8, 5], [23, 9], [19, 21], [4, 17]], [13, 20])])
add("blind_triangle_shallow", "convex_blind_pockets", 16, RECT,
    [([[-7, -5], [9, -4], [1, 10]], [12, 16])])
add("blind_hex_deep", "convex_blind_pockets", 22,
    [[-36, -24], [36, -24], [36, 24], [-36, 24]],
    [([[3, -12], [13, -10], [19, -1], [14, 10], [2, 12], [-5, 3]], [4, 22])],
    shift=(3, 2))
add("blind_two_levels", "convex_blind_pockets", 24,
    [[-44, -24], [44, -24], [44, 24], [-44, 24]],
    [([[-31, -12], [-19, -12], [-19, 2], [-31, 2]], [5, 16]),
     ([[12, 4], [30, 4], [30, 16], [12, 16]], [11, 24])])
add("blind_close_to_boundary", "convex_blind_pockets", 14,
    [[-30, -18], [30, -18], [30, 18], [-30, 18]],
    [([[20, -12], [27, -9], [25, 8], [18, 5]], [6, 14])])

# Concavity is deliberately scored as its own extension, including concave
# stock sections and a concave internal cutter. No curved/3D Boolean cases.
add("concave_c_stock_two_holes", "concave_stock_or_cut", 13,
    [[-24, -20], [24, -20], [24, -13], [-11, -13], [-11, 13], [24, 13],
     [24, 20], [-24, 20]],
    [([[-21, -17], [-16, -17], [-16, -11], [-21, -11]], [0, 13]),
     ([[-21, 11], [-16, 11], [-16, 17], [-21, 17]], [0, 13])])
add("concave_l_stock_pocket", "concave_stock_or_cut", 17,
    [[-28, -22], [28, -22], [28, -12], [-8, -12], [-8, 22], [-28, 22]],
    [([[-23, 4], [-12, 4], [-12, 13], [-18, 13], [-18, 18], [-23, 18]], [9, 17])])
add("concave_cut_in_rectangle", "concave_stock_or_cut", 11, RECT,
    [([[-19, -12], [-5, -12], [-5, -7], [-11, -7], [-11, 1], [-19, 1]], [2, 9])])
add("concave_u_stock_offset_hole", "concave_stock_or_cut", 19,
    [[-32, -22], [32, -22], [32, 22], [18, 22], [18, -8], [-18, -8],
     [-18, 22], [-32, 22]],
    [([[21, -3], [27, -3], [27, 4], [21, 4]], [0, 19])], shift=(-2, 3))
add("concave_stock_blind_two_profiles", "concave_stock_or_cut", 21,
    [[-40, -24], [40, -24], [40, -16], [-14, -16], [-14, 16], [40, 16],
     [40, 24], [-40, 24]],
    [([[-34, -20], [-25, -20], [-25, -12], [-34, -12]], [0, 21]),
     ([[-33, 12], [-23, 12], [-23, 20], [-33, 20]], [8, 21])])

# Seeded randomized convex stress series. Each stock is a slanted integer
# affine image of a regular n-gon. Homothetic cutters are strictly contained;
# translated cutters are checked by the runner's exact Fraction predicates.
rng = random.Random(20261126)
for case_index in range(30):
    n = 5 + ((case_index * 7 + rng.randrange(0, 8)) % 8)
    radius = 100
    base = [[round(radius * math.cos(2 * math.pi * k / n)),
             round(radius * math.sin(2 * math.pi * k / n))] for k in range(n)]
    # A deterministic integer shear and quarter-turn create non-axis-aligned
    # stock while preserving exact convexity and integral coordinates.
    shear = rng.choice([-2, -1, 0, 1, 2])
    turn = rng.randrange(4)
    shift = [rng.randint(-250, 250), rng.randint(-250, 250)]
    def transform(p, scale):
        x, y = p[0] * scale, p[1] * scale
        x += shear * y
        for _ in range(turn):
            x, y = -y, x
        return [x + shift[0], y + shift[1]]
    outer = [transform(p, 10) for p in base]
    feature_count = 1 + (case_index % 3)
    offsets = {1: [(0, 0)], 2: [(-280, 0), (280, 0)],
               3: [(-300, -180), (300, -180), (0, 320)]}[feature_count]
    cuts = []
    for j, (ox, oy) in enumerate(offsets):
        pts = []
        for p in base:
            # Apply the same slant/rotation to cutter and stock. Offset is
            # transformed too so it remains in the stock's local frame.
            x, y = p[0] * 2 + ox, p[1] * 2 + oy
            x += shear * y
            for _ in range(turn):
                x, y = -y, x
            pts.append([x + shift[0], y + shift[1]])
        h = 9 + (case_index % 10)
        if case_index % 2 == 0:
            interval = [0, h]
            family = "seeded_convex_through_holes"
        else:
            floor = 2 + (case_index % 4)
            interval = [floor, h]
            family = "seeded_convex_blind_pockets"
        cuts.append((pts, interval))
    family = ("seeded_convex_through_holes" if case_index % 2 == 0
              else "seeded_convex_blind_pockets")
    add(f"seeded_convex_{case_index:02d}", family, h, outer, cuts,
        shift=(rng.randint(-2, 2), rng.randint(-2, 2)),
        notes=f"seed=20261126; n={n}; cuts={feature_count}; shear={shear}; quarter_turn={turn}")

data = json.dumps({"schema": "polygonal_paired_eval_v1", "generator_seed": 20261126,
                   "coordinate_contract": "integer millimetres, each coordinate in [-4096,4096]",
                   "cases": CASES}, indent=2, allow_nan=False) + "\n"
OUT.write_text(data)
print(f"cases={len(CASES)} sha256={hashlib.sha256(OUT.read_bytes()).hexdigest()}")
