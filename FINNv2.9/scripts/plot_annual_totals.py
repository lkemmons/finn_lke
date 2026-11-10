#!/usr/bin/env python
"""Publication figures of FINN annual emission totals (from compute_annual_totals.py).

Colour = FINN version, line style = fire-detection input (mod solid / modvrs dashed).
Years with daily coverage below MIN_COVERAGE (finn_datasets.py) are not plotted.
Writes <out>.pdf (vector, editable text) and <out>.png (300 dpi).

Examples:
    python plot_annual_totals.py --species CO                    # global, single column
    python plot_annual_totals.py --species CO --width double
    python plot_annual_totals.py --species OC --region SHAF
"""
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MultipleLocator

import finn_datasets as cfg

VERSION_COLORS = {"FINNv2.5": "#2a78d6", "FINNv2.9.2": "#eb6834", "FINNv2.9.1": "#1baf7a"}
SIM_STYLE = {  # linestyle, marker face ("color" = filled, "white" = open)
    "mod":    dict(ls="-",  mfc=None),
    "modvrs": dict(ls="--", mfc="white"),
}
SIM_LABELS = {"mod": "MODIS", "modvrs": "MODIS+VIIRS"}
REGION_LONG = {
    "Global": "Global", "BONA": "Boreal North America", "TENA": "Temperate North America",
    "CEAM": "Central America", "NHSA": "N. Hemisphere South America",
    "SHSA": "S. Hemisphere South America", "EURO": "Europe", "MIDE": "Middle East",
    "NHAF": "N. Hemisphere Africa", "SHAF": "S. Hemisphere Africa", "BOAS": "Boreal Asia",
    "CEAS": "Central Asia", "SEAS": "Southeast Asia", "EQAS": "Equatorial Asia",
    "AUST": "Australia and New Zealand",
}
WIDTHS = {"single": (3.5, 2.7, 8), "double": (7.2, 3.6, 9)}   # inches, inches, font pt


def set_style(fontsize):
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": fontsize, "axes.labelsize": fontsize, "legend.fontsize": fontsize - 1,
        "xtick.labelsize": fontsize - 1, "ytick.labelsize": fontsize - 1,
        "axes.linewidth": 0.6, "axes.edgecolor": "#333333",
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "xtick.major.size": 3, "xtick.minor.size": 1.5,
        "ytick.major.size": 3,
        "axes.grid": True, "axes.grid.axis": "y", "grid.color": "#dddddd", "grid.linewidth": 0.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
    })


def year_step(years, step):
    """Major x-tick spacing: `step` years, but every year for short ranges (<= 6 years)."""
    return 1 if years[1] - years[0] < 6 else step


def plot_region(ax, ann, region, years, ms=3, lw=1.2):
    """Draw the 4 dataset lines for one region onto ax; return legend handles."""
    handles = []
    xlim = (years[0] - 0.6, years[1] + 0.6)
    for ds in cfg.DATASETS:
        d = ann[(ann.dataset == ds["name"]) & (ann.region == region)
                & (ann.coverage >= cfg.MIN_COVERAGE)
                & ann.year.between(*years)].sort_values("year")
        if d.empty:
            continue
        c = VERSION_COLORS[ds["version"]]
        st = SIM_STYLE[ds["sim"]]
        (h,) = ax.plot(d.year, d.Tg, color=c, ls=st["ls"], lw=lw,
                       marker="o", ms=ms, mew=0.8, mec=c, mfc=st["mfc"] or c,
                       label=f"{ds['version']} {SIM_LABELS[ds['sim']]}",
                       zorder=3 if ds["sim"] == "modvrs" else 2)
        handles.append(h)
    ax.set_xlim(*xlim)
    ax.set_ylim(bottom=0)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    return handles


def plot_grid(ann, species, years):
    """4x4 panels: World upper left, then the 14 GFED regions in reading order,
    legend in the lower-right panel. Each panel has its own y scale."""
    set_style(8)
    mpl.rcParams.update({"xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, axes = plt.subplots(4, 4, figsize=(7.2, 7.0), sharex=True, layout="constrained")
    axes = axes.ravel()
    legend_idx = 15
    panels = dict(enumerate(REGION_LONG))          # 0 = Global, 1..14 = GFED regions

    handles = None
    for i, region in panels.items():
        ax = axes[i]
        h = plot_region(ax, ann, region, years, ms=2.2, lw=0.9)
        handles = handles or h
        ymax = max(line.get_ydata().max() for line in ax.get_lines())
        ax.set_ylim(0, ymax * 1.2)                 # headroom for the panel label
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4, min_n_ticks=3))
        ax.xaxis.set_major_locator(MultipleLocator(year_step(years, 5)))
        ax.text(0.03, 0.97, "World" if region == "Global" else region, transform=ax.transAxes,
                ha="left", va="top", fontweight="bold", fontsize=8,
                bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8), zorder=5)
    # x tick labels only on the lowest data panel of each column
    for col in range(4):
        rows = [r for r in range(4) if (r * 4 + col) in panels]
        for r in rows:
            axes[r * 4 + col].tick_params(labelbottom=(r == max(rows)))
    # legend panel
    lax = axes[legend_idx]
    lax.axis("off")
    lax.legend(handles=handles, loc="center", frameon=False, handlelength=2.6,
               labelspacing=0.9, fontsize=7.5)
    label = cfg.SPECIES[species]["label"]
    fig.supylabel(f"{label} emissions (Tg yr$^{{-1}}$)", fontsize=9)
    fig.supxlabel("Year", fontsize=9)
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, choices=list(cfg.SPECIES))
    p.add_argument("--region", default="Global", choices=list(REGION_LONG) + ["all"],
                   help="region code, or 'all' for the 4x4 World + 14 GFED regions figure")
    p.add_argument("--totals-dir", type=Path, default=Path("totals"))
    p.add_argument("--out-dir", type=Path, default=Path("figures"))
    p.add_argument("--width", choices=list(WIDTHS), default="single")
    p.add_argument("--years", nargs=2, type=int, default=[2002, 2025])
    args = p.parse_args()

    ann = pd.read_csv(args.totals_dir / f"annual_totals_{args.species}.csv")
    yrs = f"{args.years[0]}-{args.years[1]}"
    if args.region == "all":
        fig = plot_grid(ann, args.species, args.years)
        stem = args.out_dir / f"finn_annual_{args.species}_regions_{yrs}"
    else:
        w, h, fs = WIDTHS[args.width]
        set_style(fs)
        fig, ax = plt.subplots(figsize=(w, h), layout="constrained")
        handles = plot_region(ax, ann, args.region, args.years)
        ax.xaxis.set_major_locator(MultipleLocator(year_step(args.years, 4 if args.width == "single" else 2)))
        label = cfg.SPECIES[args.species]["label"]
        where = "Global" if args.region == "Global" else args.region
        ax.set_ylabel(f"{where} {label} emissions (Tg yr$^{{-1}}$)")
        ax.set_xlabel("Year")
        fig.legend(handles=handles, loc="outside upper center",
                   ncols=2 if len(handles) == 4 else min(len(handles), 3), frameon=False,
                   handlelength=2.6, columnspacing=1.2, handletextpad=0.5)
        stem = args.out_dir / f"finn_annual_{args.species}_{args.region}_{args.width}_{yrs}"

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}")
    print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
