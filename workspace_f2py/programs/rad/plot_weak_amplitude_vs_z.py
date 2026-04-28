#!/usr/bin/env python3

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


RAD_DIR = Path(__file__).resolve().parent
INPUT_FILE = RAD_DIR / "outputs" / "results" / "amplitudes" / "amplitude_calculation_results.txt"
OUTPUT_DIR = RAD_DIR / "outputs" / "analysis"
OUTPUT_FILE = OUTPUT_DIR / "weak_amplitude_vs_Z.png"

SUMMARY_RE = re.compile(
    r"DATA:\s+Z=(?P<z>[-+0-9.]+)\s+N=(?P<n>[-+0-9.]+)\s+R=(?P<r>[-+0-9.]+)\s+"
    r"(?P<initial>\S+)->(?P<final>\S+)\s+Inter:(?P<inter>\S+)\s+"
    r"A_S=\s*(?P<a_s>[-+0-9.Ee]+)\s+A_W=\s*(?P<a_w>[-+0-9.Ee]+)\s+eta=\s*(?P<eta>[-+0-9.Ee]+)"
)


def load_latest_entries(path: Path) -> list[dict[str, object]]:
    entries: dict[tuple[int, str, str, str], dict[str, object]] = {}

    for line in path.read_text().splitlines():
        match = SUMMARY_RE.search(line)
        if not match:
            continue

        z = int(round(float(match.group("z"))))
        initial = match.group("initial")
        final = match.group("final")
        inter = match.group("inter")

        key = (z, initial, final, inter)
        entries[key] = {
            "Z": z,
            "initial": initial,
            "final": final,
            "intermediate": inter,
            "A_W": float(match.group("a_w")),
            "eta": float(match.group("eta")),
        }

    return sorted(entries.values(), key=lambda item: int(item["Z"]))


def build_plot(entries: list[dict[str, object]], output_file: Path) -> None:
    if not entries:
        raise ValueError("No amplitude entries were found to plot.")

    z_values = [int(item["Z"]) for item in entries]
    a_w_values = [float(item["A_W"]) for item in entries]
    log_abs_a_w_values = [float("-inf") if value == 0.0 else math.log10(abs(value)) for value in a_w_values]
    labels = [f"{item['initial']}->{item['final']}" for item in entries]

    fit_x = [math.log10(z) for z in z_values]
    fit_y = log_abs_a_w_values
    mean_x = sum(fit_x) / len(fit_x)
    mean_y = sum(fit_y) / len(fit_y)
    fitted_power = sum((x - mean_x) * (y - mean_y) for x, y in zip(fit_x, fit_y)) / sum(
        (x - mean_x) ** 2 for x in fit_x
    )
    fitted_intercept = mean_y - fitted_power * mean_x

    z3_power = 3.0
    z3_intercept = mean_y - z3_power * mean_x
    z3_line_x = [z3_intercept + z3_power * math.log10(z) for z in z_values]

    output_file.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.scatter(log_abs_a_w_values, z_values, s=70, color="#0b6e4f", edgecolors="black", linewidths=0.8, zorder=3)
    ax.plot(log_abs_a_w_values, z_values, color="#1d4e89", linewidth=1.8, zorder=2, label="Data trend")
    ax.plot(z3_line_x, z_values, color="#b55d00", linewidth=1.8, linestyle=":", zorder=1, label=r"Reference: $|A_W|\propto Z^3$")

    for x_val, y_val, label, a_w_val in zip(log_abs_a_w_values, z_values, labels, a_w_values):
        ax.annotate(label, (x_val, y_val), xytext=(6, 4), textcoords="offset points", fontsize=9)
        ax.annotate(f"{a_w_val:+.1e}", (x_val, y_val), xytext=(6, -11), textcoords="offset points", fontsize=8, alpha=0.8)

    ax.set_xlabel(r"$\log_{10}(|A_W|)$")
    ax.set_ylabel("Atomic number Z")
    ax.set_title(r"Comparison of $\log_{10}(|A_W|)$ Across Available Systems")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="lower right", frameon=True)
    ax.text(
        0.03,
        0.97,
        fr"Fit from data: $|A_W|\propto Z^{{{fitted_power:.2f}}}$",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )

    fig.tight_layout()
    fig.savefig(output_file, dpi=200)
    plt.close(fig)


def main() -> None:
    entries = load_latest_entries(INPUT_FILE)
    build_plot(entries, OUTPUT_FILE)

    print(f"Loaded {len(entries)} unique amplitude entries from {INPUT_FILE}")
    print(f"Saved plot to {OUTPUT_FILE}")
    for item in entries:
        print(
            f"Z={item['Z']:>3}  "
            f"{item['initial']}->{item['final']:<3}  "
            f"A_W={item['A_W']:+.8e}  "
            f"eta={item['eta']:+.8e}"
        )


if __name__ == "__main__":
    main()
