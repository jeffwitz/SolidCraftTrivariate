"""Plot converging hexahedral results against fixed-partition IGA measurements."""

import json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
hexes = json.loads((root / "hexahedral-results.json").read_text())
iga = json.loads((root / "normal-refinement.json").read_text())
fig, axes = plt.subplots(2, 2, figsize=(9, 6.3), layout="constrained")
for col, rise in enumerate([0, 3]):
    for i, (hex_key, iga_key, label) in enumerate(
        [
            ("ux_mm", "loaded_end_mean_ux_mm", "Mean end displacement (mm)"),
            ("max_gap_N", "max_gap_N", "Largest interface force jump (N)"),
        ]
    ):
        ax = axes[i, col]
        for kind, color, marker in [
            ("HEX8", "#02675e", "o"),
            ("HEX20", "#984000", "s"),
        ]:
            rows = sorted(
                (r for r in hexes if r["rise_mm"] == rise and r["element"] == kind),
                key=lambda r: r["free_dofs"],
            )
            ax.plot(
                [r["free_dofs"] for r in rows],
                [r[hex_key] for r in rows],
                marker=marker,
                color=color,
                label=kind,
            )
        rows = [r for r in iga if r["rise_mm"] == rise]
        ax.plot(
            [r["free_vector_dofs"] for r in rows],
            [r[iga_key] for r in rows],
            "^--",
            color="#54317d",
            label="IGA, transverse refinement",
        )
        ax.set_xscale("log")
        ax.grid(alpha=0.2)
        if col == 0:
            ax.set_ylabel(label)
        if i == 1:
            ax.set_xlabel("Free displacement DOFs")
        if i == 0:
            ax.set_title(
                "Perpendicular tools"
                if rise == 0
                else "Inclined tool (+3 mm per 100 mm)"
            )
axes[0, 0].legend(fontsize=8)
fig.savefig(root / "hexahedral-comparison.png", dpi=170)
