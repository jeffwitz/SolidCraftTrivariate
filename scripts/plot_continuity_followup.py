"""Render the two native error-per-DOF sweeps from their published JSON."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/results/continuity-followup"
DEST = ROOT / "docs/rtd/reference/continuity-errors.svg"


def plot():
    pocket = json.loads((REPORT / "pocket-vector.json").read_text())
    corner = json.loads((REPORT / "reentrant-ritz.json").read_text())["states"]
    with plt.rc_context(
        {"font.size": 10, "svg.hashsalt": "continuity-followup", "svg.fonttype": "none"}
    ):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")
        for p, color in [(3, "#106b75"), (4, "#a34d14")]:
            for c, style in [(0, "--o"), (1, "-s")]:
                data = [r for r in pocket if r["p"] == p and r["continuity"] == c]
                axes[0].loglog(
                    [r["free_vector_dofs"] for r in data],
                    [r["relative_energy_norm"] for r in data],
                    style,
                    color=color,
                    label=f"p={p}, C{c} between charts",
                )
        axes[0].set(
            xlabel="Free displacement coefficients",
            ylabel="Relative energy-norm error",
            title="Blind pocket · vector elasticity",
        )
        for name, style, color in [
            ("uniform", "--o", "#606570"),
            ("thb", "-s", "#287844"),
        ]:
            data = [r["native"] for r in corner if r["curve"] == name]
            axes[1].loglog(
                [r["free_scalar_dofs"] for r in data],
                [r["relative_H1_seminorm"] for r in data],
                style,
                color=color,
                label=name.upper() if name == "thb" else "Uniform",
            )
        axes[1].set(
            xlabel="Free scalar coefficients",
            ylabel="Relative H¹ seminorm error",
            title="Reentrant edge · scalar harmonic field",
        )
        axes[1].yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=1))
        axes[1].yaxis.set_minor_formatter(PercentFormatter(xmax=1, decimals=1))
        for ax in axes:
            ax.grid(True, which="major", alpha=0.24)
            ax.legend(fontsize=8)
        fig.savefig(DEST, metadata={"Date": None})
        DEST.write_text("\n".join(line.rstrip() for line in DEST.read_text().splitlines()) + "\n")
        plt.close(fig)


if __name__ == "__main__":
    plot()
