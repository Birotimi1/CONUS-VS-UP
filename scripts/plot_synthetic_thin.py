#!/usr/bin/env python3
"""Plot saved forward predictions, using DATA/synthetics/LATEST_RUN.txt by default."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
BLUE = "#187a9a"
NAMES = ["hvsr", "ellip", "cph", "ugr", "rf"]
COORDS = ["f", "T", "T", "T", "t"]


def draw_vs(ax, model, bottom_km=75, depth_scale=1):
    depths = np.r_[0, np.cumsum(model[:-1, 0]), bottom_km]
    ax.plot(np.repeat(model[:, 2], 2),
            np.column_stack([depths[:-1], depths[1:]]).ravel()*depth_scale,
            color=BLUE, lw=2.2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", nargs="?", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir
    if run_dir is None:
        parent = ROOT / "DATA/synthetics"
        run_dir = parent / (parent / "LATEST_RUN.txt").read_text().strip()
    result = loadmat(run_dir / "synthetic_thin.mat", simplify_cells=True)
    tm, data, forward = (result[k] for k in ["true_model", "syn_data", "forward"])
    if not forward["all_forward_succeeded"]:
        raise ValueError("The forward run is incomplete; curves will not be plotted.")
    for name, coord in zip(NAMES, COORDS):
        x, y = np.asarray(data[name][coord]), np.asarray(data[name]["clean"])
        if x.shape != y.shape or not np.isfinite(y).all():
            raise ValueError(f"Invalid {name} output")
        np.testing.assert_array_equal(forward["pred"][name]["model_input"], forward["model"])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 13, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 3, figsize=(16.4, 10.0))
    axes = axes.ravel()
    model = tm["model_matrix"]
    draw_vs(axes[0], model)
    axes[0].set(title="(a) Prescribed shear-velocity structure", xlabel="$V_S$ (km/s)",
                ylabel="Depth (km)", xlim=(0, 5), ylim=(75, 0))
    axes[0].axhline(40, color=".5", lw=.9, ls="--")
    axes[0].text(.04, 40, " Moho: 40 km", transform=axes[0].get_yaxis_transform(),
                 va="bottom", color=".35", fontsize=9)
    inset = axes[0].inset_axes([.12, .08, .42, .32])
    draw_vs(inset, model, depth_scale=1000)
    inset.set(title="Upper 200 m", xlim=(0, 3), ylim=(200, 0), ylabel="Depth (m)")
    inset.set_xticks([0, 1.5, 3]); inset.set_yticks([0, 50, 130, 200])
    inset.tick_params(labelsize=8); inset.title.set_size(9); inset.yaxis.label.set_size(8)
    inset.grid(alpha=.18)
    titles = ["(b) Diffuse-field HVSR", "(c) Rayleigh ellipticity",
              "(d) Rayleigh phase velocity", "(e) Rayleigh group velocity",
              "(f) Synthetic P receiver function"]
    xlabels = ["Frequency (Hz)", "Period (s)", "Period (s)", "Period (s)", "Time after direct P (s)"]
    ylabels = ["H/V", "H/V", "Phase velocity (km/s)", "Group velocity (km/s)", "R/Z transfer amplitude"]
    for ax, name, coord, title, xlabel, ylabel in zip(axes[1:], NAMES, COORDS, titles, xlabels, ylabels):
        x, y = np.asarray(data[name][coord]), np.asarray(data[name]["clean"])
        ax.plot(x, y, color=BLUE, lw=1.8)
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
        if name != "rf":
            ax.set_xlim(x.min(), x.max())
    axes[1].set_xscale("log")
    axes[1].set_xticks([.2, .5, 1, 2, 5, 10, 25])
    axes[1].xaxis.set_major_formatter(ScalarFormatter())
    axes[1].tick_params(axis="x", which="minor", labelbottom=False)
    rf = forward["config"]["rf"]
    t, amp = data["rf"]["t"], data["rf"]["clean"]
    axes[5].fill_between(t, 0, amp, where=amp>=0, color=BLUE, alpha=.20)
    axes[5].fill_between(t, 0, amp, where=amp<0, color="#bf683a", alpha=.28)
    axes[5].axhline(0, color=".45", lw=.65)
    axes[5].set_xlim(*np.asarray(rf["plot_window"]).ravel())
    ps = tm["moho_primary_ps_delay_s"]
    axes[5].axvline(ps, color=".45", lw=1, ls=":")
    axes[5].text(ps+.25, .96, f"Moho Ps reference\n{ps:.2f} s", va="top", fontsize=9,
                 color=".35", transform=axes[5].get_xaxis_transform())
    axes[5].text(.98, .82, f"p = {rf['rayp_s_km']:.2f} s/km | alpha = {rf['gaussian_alpha']:.1f}",
                 transform=axes[5].transAxes, ha="right", fontsize=9)
    for ax in axes:
        ax.grid(alpha=.20)
    fig.suptitle("Known thin layered Earth: five forward observables", fontsize=20, weight="bold", y=.975)
    fig.text(.5, .935, "50 m soft cap | Moho at 40 km | Half-space starts at 70 km | One flat, isotropic, elastic model",
             ha="center", fontsize=11, color=".3")
    fig.text(.5, .02, "Clean synthetic predictions from HVf and CPS. No noise or inversion. Moho Ps marker is a ray-theory reference.",
             ha="center", fontsize=10, color=".35")
    fig.tight_layout(rect=(.015, .045, .99, .91), w_pad=2.5, h_pad=2.8)
    for extension in ["png", "pdf"]:
        fig.savefig(run_dir / f"five_observables.{extension}", dpi=240, facecolor="white")
    plt.close(fig)
    report = {"case": tm["case_name"], "model_input_identical_for_all_five": True,
              "source_layers": int(model.shape[0]), "numerical_layers": int(forward["model"].shape[0]),
              "moho_primary_ps_reference_s": float(ps), "noise_added": False, "observables": {}}
    for name in NAMES:
        y = np.asarray(data[name]["clean"])
        report["observables"][name] = {"samples": int(y.size), "minimum": float(y.min()), "maximum": float(y.max())}
    (run_dir / "validation.json").write_text(json.dumps(report, indent=2)+"\n")
    print(run_dir / "five_observables.png")


if __name__ == "__main__":
    main()
