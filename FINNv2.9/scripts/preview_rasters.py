#!/usr/bin/env python3
"""Render quick-look PNG previews of LCT and VCF GeoTIFFs.

Auto-detects the product from the filename (``modlct`` vs ``modvcf``)
and picks appropriate symbology:

  * LCT  → categorical IGBP colormap (17 classes + unclassified) with a
           legend.
  * VCF  → continuous 0..100 color ramp per band (3 bands: tree, herb,
           bare).  One PNG per band.

Reads a decimated version of the raster via rasterio's ``decimation``
read, so even a global 2-GB COG renders in a few seconds.

Example
-------
    python scripts/preview_rasters.py \\
        /glade/.../rasters/modlct_2020.tif \\
        /glade/.../rasters/modvcf_2020.tif \\
        --out-dir ./previews

Produces:
    previews/modlct_2020.png
    previews/modvcf_2020_tree.png
    previews/modvcf_2020_herb.png
    previews/modvcf_2020_bare.png
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize
from matplotlib.patches import Patch

log = logging.getLogger("preview_rasters")


# IGBP land-cover class names and NASA-standard-ish colors.  Index = value.
_LCT_CLASSES = [
    (0,  "Water",                     "#1f3f8f"),
    (1,  "Evergreen Needleleaf",      "#05450a"),
    (2,  "Evergreen Broadleaf",       "#086a10"),
    (3,  "Deciduous Needleleaf",      "#54a708"),
    (4,  "Deciduous Broadleaf",       "#78d203"),
    (5,  "Mixed Forest",              "#009900"),
    (6,  "Closed Shrublands",         "#c6b044"),
    (7,  "Open Shrublands",           "#dcd159"),
    (8,  "Woody Savannas",            "#dade48"),
    (9,  "Savannas",                  "#fbff13"),
    (10, "Grasslands",                "#b6ff05"),
    (11, "Permanent Wetlands",        "#27ff87"),
    (12, "Croplands",                 "#c24f44"),
    (13, "Urban & Built-up",          "#a5a5a5"),
    (14, "Cropland Mosaic",           "#ff6d4c"),
    (15, "Snow & Ice",                "#f0f0f0"),
    (16, "Barren",                    "#cdb380"),
    (17, "Unclassified",              "#000000"),
]
_LCT_VALUES = [v for v, _, _ in _LCT_CLASSES]
_LCT_LABELS = [n for _, n, _ in _LCT_CLASSES]
_LCT_COLORS = [c for _, _, c in _LCT_CLASSES]

# VCF band names (MODIS MOD44B order).  If the file has a different band
# count or BAND_NAMES tags, we fall back to generic "band_1", "band_2", …
_VCF_BAND_NAMES = ["tree", "herb", "bare"]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="preview_rasters.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    p.add_argument("inputs", type=Path, nargs="+",
                   help="one or more GeoTIFF files (modlct_* or modvcf_*)")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="output directory for PNGs (created if missing)")
    p.add_argument("--max-width", type=int, default=2000,
                   help="target pixel width of the decimated read "
                        "(default: %(default)s)")
    p.add_argument("--dpi", type=int, default=120,
                   help="PNG DPI (default: %(default)s)")
    p.add_argument("--force", choices=["lct", "vcf"], default=None,
                   help="override auto-detection from filename")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Decimated read
# ---------------------------------------------------------------------------

def _decimated_read(src: rasterio.DatasetReader, max_width: int,
                     band: int = 1) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Read a downsampled version of ``band`` so the longer axis fits in
    ``max_width`` pixels.  Returns (array, (west, east, south, north))."""
    w, h = src.width, src.height
    scale = max(1, max(w, h) // max_width)
    out_w, out_h = max(1, w // scale), max(1, h // scale)
    data = src.read(band, out_shape=(out_h, out_w),
                     resampling=rasterio.enums.Resampling.nearest)
    b = src.bounds
    return data, (b.left, b.right, b.bottom, b.top)


# ---------------------------------------------------------------------------
# Renderers
# ---------------------------------------------------------------------------

def render_lct(path: Path, out_png: Path, *, max_width: int, dpi: int) -> None:
    """Categorical IGBP map with a side legend."""
    with rasterio.open(path) as src:
        if src.count > 1:
            log.warning("%s has %d bands; previewing band 1 only",
                         path.name, src.count)
        data, (w, e, s, n) = _decimated_read(src, max_width, band=1)
        nodata = src.nodata

    # Mask nodata and out-of-range values
    if nodata is not None:
        data = np.where(data == nodata, 255, data)
    unique_vals = sorted(set(np.unique(data)) - {255})

    cmap = ListedColormap(_LCT_COLORS)
    bounds = np.array(_LCT_VALUES + [_LCT_VALUES[-1] + 1]) - 0.5
    norm = BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.imshow(data, extent=(w, e, s, n), cmap=cmap, norm=norm,
              interpolation="nearest", origin="upper")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title(f"{path.stem}  (LCT, IGBP classes)")
    ax.grid(True, alpha=0.25)

    # Legend — only classes that appear in the preview
    present = [i for i, v in enumerate(_LCT_VALUES) if v in unique_vals]
    handles = [Patch(facecolor=_LCT_COLORS[i], edgecolor="k",
                      linewidth=0.3, label=f"{_LCT_VALUES[i]:2d}  {_LCT_LABELS[i]}")
               for i in present]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=8, title="IGBP class", frameon=True)

    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info("  wrote %s  (%d classes present)", out_png.name, len(present))


def render_vcf(path: Path, out_dir: Path, *, max_width: int, dpi: int) -> None:
    """Continuous 0..100 ramp, one PNG per band."""
    with rasterio.open(path) as src:
        n_bands = src.count
        for i in range(1, n_bands + 1):
            band_name = (_VCF_BAND_NAMES[i - 1] if i - 1 < len(_VCF_BAND_NAMES)
                          else f"band_{i}")
            data, (w, e, s, n) = _decimated_read(src, max_width, band=i)
            nodata = src.nodata

            arr = data.astype("float32")
            if nodata is not None:
                arr = np.where(data == nodata, np.nan, arr)
            # Common MODIS MOD44B fill codes (>100) are masked
            arr = np.where(arr > 100, np.nan, arr)

            cmap_name = {
                "tree": "Greens", "herb": "YlGn", "bare": "YlOrBr",
            }.get(band_name, "viridis")
            cmap = plt.get_cmap(cmap_name).copy()
            cmap.set_bad(color="#dddddd")      # NaN / fill colour
            norm = Normalize(vmin=0, vmax=100)

            fig, ax = plt.subplots(figsize=(14, 7))
            im = ax.imshow(arr, extent=(w, e, s, n), cmap=cmap, norm=norm,
                            interpolation="nearest", origin="upper")
            ax.set_xlabel("longitude")
            ax.set_ylabel("latitude")
            ax.set_title(f"{path.stem}  (VCF {band_name}, %)")
            ax.grid(True, alpha=0.25)
            cbar = fig.colorbar(im, ax=ax, orientation="vertical",
                                 shrink=0.75, pad=0.02)
            cbar.set_label(f"{band_name} cover (%)")

            fig.tight_layout()
            out_png = out_dir / f"{path.stem}_{band_name}.png"
            fig.savefig(out_png, dpi=dpi, bbox_inches="tight",
                        facecolor="white")
            plt.close(fig)
            log.info("  wrote %s", out_png.name)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def _detect_kind(path: Path) -> str:
    """'lct' or 'vcf' from the filename; raises if neither matches."""
    stem = path.stem.lower()
    if "lct" in stem:
        return "lct"
    if "vcf" in stem:
        return "vcf"
    raise RuntimeError(
        f"cannot infer product type from {path.name!r}; pass --force lct|vcf")


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ok, fail = 0, 0
    for path in args.inputs:
        if not path.exists():
            log.error("not found: %s", path); fail += 1; continue
        kind = args.force or _detect_kind(path)
        log.info("rendering %s (%s)", path.name, kind.upper())
        try:
            if kind == "lct":
                render_lct(path, args.out_dir / f"{path.stem}.png",
                            max_width=args.max_width, dpi=args.dpi)
            else:
                render_vcf(path, args.out_dir,
                            max_width=args.max_width, dpi=args.dpi)
            ok += 1
        except Exception as e:
            log.error("failed on %s: %s", path.name, e, exc_info=args.verbose)
            fail += 1

    log.info("done — %d ok, %d failed", ok, fail)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

