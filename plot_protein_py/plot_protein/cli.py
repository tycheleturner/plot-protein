#!/usr/bin/env python3
# plot-protein Python code

"""
Protein plotting script (built using the initial plot-protein R framework)

Features:
  - Plot protein line, architecture (e.g., domains), post-translational modification sites.
  - Plot one or two mutation cohorts:
      * --mutations -> top (or single) cohort (always above the line)
      * --mutations_bottom -> bottom cohort (optional; plotted below the line)
  - Mutation files can have 5, 6, or 7 columns:
      5 cols: ProteinId, GeneName, Position, RefAA, AltAA
      6 cols: ProteinId, GeneName, Position, RefAA, AltAA, Annotation
      7 cols: ProteinId, GeneName, Position, RefAA, AltAA, Annotation, Score
  - Coloring modes (--color-by):
      * auto (default): use Annotation colors if present, else by cohort
      * annotation: force Annotation-based colors
      * score: color by numeric Score (continuous colormap + colorbar), this could be a conservation or severity score
      * cohort: ignore annotation/score, color by cohort (top vs bottom)
  - Recurrent positions (>= 2 occurrences within a cohort) use diamond markers.
  - Optional hiding of domains and PTMs.
  - Architecture (-a) and PTM (-p) files are optional; if omitted, those layers are not drawn.
  - Output format: PDF / PNG / SVG.
  - Optional faceting by domain regions: --facet-domains (this just means we show the domains with several plots at once)
  - Themes: --theme light/dark
  - Palettes: --palette default/colorblind
  - DPI control for PNGs: --dpi

Extra controls:
  - Jitter: --jitter auto/off, with --jitter-window and --jitter-amplitude (in case there are lots of variants this can make it easier to see them all)
  - Grid: --grid
  - Point size: --point-size
  - Title override: --title
  - Filtering: --include-annotations, --min-score
  - Export annotation colors: --annotation-colors-out
  - Export score colormap bins: --score-colors-out
"""

# --- imports ---
import argparse
import sys
from pathlib import Path
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# --- base visual style settings (light default) ---

MODERN_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "Helvetica", "sans-serif"],
    "axes.edgecolor": "#cccccc",
    "axes.linewidth": 0.8,
    "xtick.color": "#000000",
    "ytick.color": "#000000",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}
plt.rcParams.update(MODERN_STYLE)

# These will be overridden by configure_style() depending on theme/palette
COLOR_PROTEIN_LINE = "#4a4a4a"
COLOR_TOP_COHORT = "#1f77b4"
COLOR_BOTTOM_COHORT = "#d62728"
COLOR_PTM = "#E69F00"
COLOR_DOMAIN = "#009E73"
COLOR_TEXT = "#000000"
COLOR_AXIS_SPINE = "#000000"
COLOR_POINT_EDGE = "#ffffff"

# --- argument parsing ---
def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        prog="plot_protein",
        description="Plot protein mutations, domains, and post-translational modifications.",
    )

    ap.add_argument(
        "-m", "--mutations",
        default="mutationFile.txt",
        help=(
            "Primary (top) cohort mutation file. Whitespace-delimited (tabs and/or spaces), "
            "5/6/7 cols: ProteinId, GeneName, ProteinPositionOfMutation, "
            "ReferenceAminoAcid, AlternateAminoAcid[, Annotation[, Score]]. NO HEADER."
        ),
    )
    ap.add_argument(
        "--mutations_bottom", "-m2",
        default=None,
        help=(
            "Optional bottom cohort mutation file. Same format as --mutations. "
            "If provided, this cohort is plotted below the protein line."
        ),
    )
    ap.add_argument(
        "--mutations-name",
        default="Cohort 1",
        help="Label for the primary/top mutation group. Default: 'Cohort 1'.",
    )
    ap.add_argument(
        "--mutations-bottom-name",
        default="Cohort 2",
        help="Label for the bottom mutation group (if provided). Default: 'Cohort 2'.",
    )
    ap.add_argument(
        "-a", "--architecture",
        default=None,
        help=(
            "Optional protein architecture file. Tab-delimited, 3 columns "
            "with header: architecture_name, start_site, end_site. "
            "If omitted, domains are not drawn (and faceting is disabled)."
        ),
    )
    ap.add_argument(
        "-p", "--posttranslational",
        default=None,
        help=(
            "Optional post-translational modification file. "
            "Tab-delimited, one column 'site' with header. "
            "If omitted, PTM sites are not drawn."
        ),
    )
    ap.add_argument(
        "-l", "--length",
        type=float,
        default=100,
        help="Protein length (REQUIRED).",
    )
    ap.add_argument(
        "-n", "--name",
        default="Test",
        help="Name of your query/study. Default: 'Test'.",
    )
    ap.add_argument(
        "-t", "--ticksize",
        type=float,
        default=10,
        help="Size of ticks on x-axis. This is dynamic with protein size but can be set by the user. Default: 10.",
    )
    ap.add_argument(
        "-s", "--showlabels",
        choices=["yes", "no"],
        default="no",
        help="Option to show mutation labels (yes/no). Default: no.",
    )
    ap.add_argument(
        "-z", "--zoom",
        choices=["yes", "no"],
        default="no",
        help="Option to zoom in somewhere in the protein (yes/no). Default: no.",
    )
    ap.add_argument(
        "-b", "--zoomstart",
        type=float,
        default=1,
        help="Starting AA position for zoom. Used if --zoom yes. Default: 1.",
    )
    ap.add_argument(
        "-c", "--zoomend",
        type=float,
        default=10,
        help="Ending AA position for zoom. Used if --zoom yes. Default: 10.",
    )
    ap.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Output filename. If it has .pdf/.png/.svg, that determines format. "
            "Otherwise the extension is added based on --format / default."
        ),
    )
    ap.add_argument(
        "--format",
        choices=["pdf", "png", "svg"],
        default=None,
        help="Output format: pdf, png, or svg. Default: infer from --output or 'pdf'.",
    )
    ap.add_argument(
        "--hide-architecture",
        action="store_true",
        help="Hide protein architecture domains.",
    )
    ap.add_argument(
        "--hide-ptms",
        action="store_true",
        help="Hide post-translational modification sites.",
    )
    ap.add_argument(
        "--facet-domains",
        action="store_true",
        help="Facet by domain region from the architecture file (one panel per domain).",
    )
    ap.add_argument(
        "--color-by",
        choices=["auto", "annotation", "score", "cohort"],
        default="auto",
        help=(
            "How to color mutation points: "
            "'auto' (default: use Annotation if present, else by cohort), "
            "'annotation', 'score' (requires Score column), or 'cohort'."
        ),
    )
    ap.add_argument(
        "--include-annotations",
        nargs="+",
        default=None,
        help="Only plot mutations whose Annotation is in this list (e.g. damaging LoF missense).",
    )
    ap.add_argument(
        "--min-score",
        type=float,
        default=None,
        help="Only plot mutations with Score >= this value.",
    )
    ap.add_argument(
        "--theme",
        choices=["light", "dark"],
        default="light",
        help='Plot theme: "light" (default) or "dark".',
    )
    ap.add_argument(
        "--palette",
        choices=["default", "colorblind"],
        default="default",
        help='Color palette: "default" or "colorblind".',
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for raster outputs (PNG). Default: 300.",
    )
    ap.add_argument(
        "--jitter",
        choices=["auto", "off"],
        default="auto",
        help="Vertical jitter for nearby mutations: 'auto' (default) or 'off'. This is helpful when mutations are very close to each other.",
    )
    ap.add_argument(
        "--jitter-window",
        type=float,
        default=5.0,
        help="Window (AA) within which mutations are considered 'nearby' for jitter (default: 5).",
    )
    ap.add_argument(
        "--jitter-amplitude",
        type=float,
        default=0.005,
        help="Maximum vertical jitter offset (default: 0.005).",
    )
    ap.add_argument(
        "--grid",
        action="store_true",
        help="Show light vertical grid lines on the x-axis.",
    )
    ap.add_argument(
        "--point-size",
        type=float,
        default=30.0,
        help="Base point size for mutation markers (default: 30).",
    )
    ap.add_argument(
        "--title",
        default=None,
        help="Override main plot title. If not set, a gene/protein-based title is used by default.",
    )
    ap.add_argument(
        "--annotation-colors-out",
        default=None,
        help="Optional TSV to write Annotation -> color mapping when annotation-based colors are used.",
    )
    ap.add_argument(
        "--score-colors-out",
        default=None,
        help="Optional TSV to write Score colormap bins (score_min/max/color) when color mode is 'score'.",
    )
    ap.add_argument(
        "--architecture-labels",
        choices=["yes", "no"],
        default="no",
        help="Show text labels for domains on the main plot (default: no).",
    )
    ap.add_argument(
        "--version",
        action="version",
        version="plot_protein 4.0.0",
    )

    args = ap.parse_args(argv)

    if args.length <= 0:
        ap.error("Protein length must be positive.")

    return args

# --- theme / palette configuration ---
def configure_style(args):
    global COLOR_PROTEIN_LINE, COLOR_TOP_COHORT, COLOR_BOTTOM_COHORT
    global COLOR_PTM, COLOR_DOMAIN, COLOR_TEXT, COLOR_AXIS_SPINE, COLOR_POINT_EDGE

    if args.theme == "dark":
        dark_style = {
            "figure.facecolor": "#111111",
            "axes.facecolor": "#111111",
            "axes.edgecolor": "#666666",
            "xtick.color": "#dddddd",
            "ytick.color": "#dddddd",
            "axes.labelcolor": "#ffffff",
            "axes.titlecolor": "#ffffff",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "Helvetica", "sans-serif"],
        }
        plt.rcParams.update(dark_style)
        COLOR_TEXT = "#f0f0f0"
        COLOR_PROTEIN_LINE = "#e0e0e0"
        COLOR_AXIS_SPINE = "#888888"
        COLOR_POINT_EDGE = "#ffffff"
    else:
        plt.rcParams.update(MODERN_STYLE)
        COLOR_TEXT = "#000000"
        COLOR_PROTEIN_LINE = "#4a4a4a"
        COLOR_AXIS_SPINE = "#000000"
        COLOR_POINT_EDGE = "#d0d0d0"

    if args.palette == "colorblind":
        COLOR_TOP_COHORT = "#0072B2"
        COLOR_BOTTOM_COHORT = "#D55E00"
        COLOR_PTM = "#CC79A7"
        COLOR_DOMAIN = "#009E73"
    else:
        COLOR_TOP_COHORT = "#1f77b4"
        COLOR_BOTTOM_COHORT = "#d62728"
        COLOR_PTM = "#E69F00"
        COLOR_DOMAIN = "#009E73"
        if args.theme == "light":
            COLOR_BOTTOM_COHORT = "#b54b3a"

# --- IO helpers ---
def read_mutation_file(path, label_desc):
    """
    Read a mutation file with 5, 6, or 7 columns, no header:
      5 cols: ProteinId, GeneName, Position, RefAA, AltAA
      6 cols: ProteinId, GeneName, Position, RefAA, AltAA, Annotation
      7 cols: ProteinId, GeneName, Position, RefAA, AltAA, Annotation, Score
    """
    cols_5 = ["ProteinId", "GeneName", "Position", "RefAA", "AltAA"]
    cols_6 = ["ProteinId", "GeneName", "Position", "RefAA", "AltAA", "Annotation"]
    cols_7 = ["ProteinId", "GeneName", "Position", "RefAA", "AltAA", "Annotation", "Score"]

    if path is None:
        return pd.DataFrame(columns=cols_7)

    path = Path(path)
    if not path.is_file():
        sys.exit(f"{label_desc}: file not found: {path}")

    try:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            engine="python",
            header=None,
            dtype=str,
        )
    except Exception as e:
        sys.exit(f"Error reading {label_desc}: {e}")

    if df.shape[1] == 5:
        df.columns = cols_5
        df["Annotation"] = pd.NA
        df["Score"] = pd.NA
    elif df.shape[1] == 6:
        df.columns = cols_6
        df["Score"] = pd.NA
    elif df.shape[1] == 7:
        df.columns = cols_7
    else:
        sys.exit(
            f"{label_desc} must have exactly 5, 6, or 7 columns:\n"
            " 5: ProteinId, GeneName, Position, RefAA, AltAA\n"
            " 6: ProteinId, GeneName, Position, RefAA, AltAA, Annotation\n"
            " 7: ProteinId, GeneName, Position, RefAA, AltAA, Annotation, Score"
        )

    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    if "Score" in df.columns:
        df["Score"] = pd.to_numeric(df["Score"], errors="coerce")

    return df


def read_optional_architecture(path):
    cols = ["architecture_name", "start_site", "end_site"]
    if path is None:
        return pd.DataFrame(columns=cols)

    p = Path(path)
    if not p.is_file():
        sys.exit(f"Architecture file not found: {p}")

    try:
        pa = pd.read_csv(
            p,
            sep="\t",
            header=0,
            dtype=str,
        )
    except Exception as e:
        sys.exit(f"Error reading protein architecture file: {e}")

    if not all(col in pa.columns for col in cols):
        sys.exit(
            "Architecture file must contain header columns: "
            + ", ".join(cols)
        )

    pa["start_site"] = pd.to_numeric(pa["start_site"], errors="coerce")
    pa["end_site"] = pd.to_numeric(pa["end_site"], errors="coerce")

    return pa


def read_optional_ptm(path):
    cols = ["site"]
    if path is None:
        return pd.DataFrame(columns=cols)

    p = Path(path)
    if not p.is_file():
        sys.exit(f"Post-translational modification file not found: {p}")

    try:
        pt = pd.read_csv(
            p,
            sep="\t",
            header=0,
            dtype=str,
        )
    except Exception as e:
        sys.exit(f"Error reading post-translational modification file: {e}")

    if "site" not in pt.columns:
        sys.exit("Post-translational modification file must have a 'site' column.")

    pt["site"] = pd.to_numeric(pt["site"], errors="coerce")

    return pt


def read_inputs(args):
    mut_top = read_mutation_file(args.mutations, "Top/primary cohort mutation file (--mutations)")
    mut_bottom = read_mutation_file(args.mutations_bottom, "Bottom cohort mutation file (--mutations_bottom)")
    pa = read_optional_architecture(args.architecture)
    pt = read_optional_ptm(args.posttranslational)
    return mut_top, mut_bottom, pa, pt


# --- concat helper to avoid pandas FutureWarning ---
def concat_mutation_dfs(mut_top, mut_bottom):
    """
    Safely concatenate mutation dataframes, skipping empty ones,
    to avoid pandas FutureWarning about empty/all-NA entries.
    """
    frames = []
    for df in (mut_top, mut_bottom):
        if df is not None and not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# --- filtering helpers ---
def filter_mutations(df, args):
    """
    Filter mutations by Annotation and Score, if requested.
    After filtering, reset the index so .loc[0] is always valid if not empty.
    """
    if df is None or df.empty:
        return df

    out = df.copy()

    if args.include_annotations is not None and "Annotation" in out.columns:
        allowed = {str(a) for a in args.include_annotations}
        out = out[out["Annotation"].astype(str).isin(allowed)]

    if args.min_score is not None and "Score" in out.columns:
        out = out[pd.to_numeric(out["Score"], errors="coerce") >= args.min_score]

    # Ensure contiguous 0..N-1 indices so code using .loc[0] works
    return out.reset_index(drop=True)


def write_annotation_color_map(path, color_map):
    """
    Write Annotation -> color mapping to a TSV file.

    If color_map is empty, we still write a header-only file.
    """
    if not path:
        return
    try:
        with open(path, "w") as fh:
            fh.write("Annotation\tColor\n")
            for ann, col in sorted(color_map.items()):
                fh.write(f"{ann}\t{col}\n")
    except Exception as e:
        sys.stderr.write(f"Warning: could not write annotation color map to {path}: {e}\n")


def write_score_colormap(path, score_norm, score_cmap, n_bins=20):
    """
    Write score colormap bins to TSV.

    Columns:
      bin  score_min  score_max  score_mid  color_hex
    """
    if not path or score_norm is None or score_cmap is None:
        return

    vmin = score_norm.vmin
    vmax = score_norm.vmax
    if vmin is None or vmax is None:
        return

    try:
        with open(path, "w") as fh:
            fh.write("bin\tscore_min\tscore_max\tscore_mid\tcolor_hex\n")

            if vmax == vmin:
                s = float(vmin)
                rgba = score_cmap(score_norm(s))
                hex_color = mcolors.to_hex(rgba)
                fh.write(f"0\t{s}\t{s}\t{s}\t{hex_color}\n")
                return

            edges = np.linspace(vmin, vmax, n_bins + 1)
            for i in range(n_bins):
                s_min = float(edges[i])
                s_max = float(edges[i + 1])
                s_mid = 0.5 * (s_min + s_max)
                rgba = score_cmap(score_norm(s_mid))
                hex_color = mcolors.to_hex(rgba)
                fh.write(f"{i}\t{s_min}\t{s_max}\t{s_mid}\t{hex_color}\n")
    except Exception as e:
        sys.stderr.write(f"Warning: could not write score colormap to {path}: {e}\n")


# --- range, out-of-range warnings ---
def compute_ranges(args, mut_top, mut_bottom):
    protein_length = float(args.length)

    if args.zoom == "yes":
        z_start = float(args.zoomstart)
        z_end = float(args.zoomend)
        if z_start > z_end:
            z_start, z_end = z_end, z_start
        z_start = max(1.0, z_start)
        z_end = min(protein_length, z_end)

        span = max(z_end - z_start, 1.0)
        pad_left = max(2.0, span * 0.05)
        pad_right = max(1.0, span * 0.02)
        xlim_region = (
            max(0.0, z_start - pad_left),
            min(protein_length, z_end + pad_right),
        )
    else:
        pad_left = max(5.0, protein_length * 0.05)
        pad_right = max(2.0, protein_length * 0.02)
        xlim_region = (
            max(0.0, 1.0 - pad_left),
            protein_length + pad_right,
        )

    combined = concat_mutation_dfs(mut_top, mut_bottom)
    combined = combined[combined["Position"].notna()] if not combined.empty else combined
    if not combined.empty:
        out_of_range = combined[
            (combined["Position"] < 1) | (combined["Position"] > protein_length)
        ]
        if not out_of_range.empty:
            rows = (out_of_range.index + 1).tolist()
            sys.stderr.write(
                "Warning: some mutations fall outside the protein length "
                f"[1, {int(protein_length)}]: combined rows {rows}\n"
            )

    return protein_length, xlim_region


# --- output path / format ---
def get_output_path_and_format(args, mut_top):
    fmt = args.format

    if args.output:
        out_path = Path(args.output)
        ext = out_path.suffix.lower()
        if ext in [".pdf", ".png", ".svg"]:
            fmt = fmt or ext.lstrip(".")
        else:
            fmt = fmt or "pdf"
            out_path = out_path.with_suffix(f".{fmt}")
    else:
        if mut_top.empty:
            gene_name = "protein"
        else:
            gene_name = str(mut_top.loc[0, "GeneName"])
            if gene_name in [None, "", "nan"]:
                gene_name = "protein"
        fmt = fmt or "pdf"
        out_path = Path(f"{gene_name}_protein_plot.{fmt}")

    if fmt not in ["pdf", "png", "svg"]:
        fmt = "pdf"

    return out_path, fmt


# --- recurrence, annotation colors, score normalizer ---
def split_by_recurrence(df):
    if df.empty or "Position" not in df.columns:
        return df, df

    pos = pd.to_numeric(df["Position"], errors="coerce")
    counts = pos.value_counts(dropna=False)
    rec_mask = pos.map(lambda v: counts.get(v, 0) > 1)

    nonrec = df[~rec_mask]
    rec = df[rec_mask]
    return nonrec, rec


def compute_cluster_jitter(df, base_y, window, amplitude):
    """
    Compute vertical jitter for mutations in a cohort:

      - Positions that are isolated stay at base_y.
      - Positions that have neighbors within 'window' amino acids
        are treated as a cluster and get small vertical offsets.
      - All mutations at the exact same position share the same y.
    """
    if df is None or df.empty or "Position" not in df.columns:
        return {idx: base_y for idx in df.index} if df is not None else {}

    y_map = {idx: base_y for idx in df.index}

    pos_series = pd.to_numeric(df["Position"], errors="coerce")
    valid = pos_series.notna()
    if valid.sum() == 0:
        return y_map

    grouped = {}
    for idx, pos in pos_series[valid].items():
        p = float(pos)
        grouped.setdefault(p, []).append(idx)

    positions = sorted(grouped.keys())
    if not positions:
        return y_map

    clusters = []
    current_cluster = [positions[0]]
    prev_pos = positions[0]

    for p in positions[1:]:
        if abs(p - prev_pos) <= window:
            current_cluster.append(p)
        else:
            clusters.append(current_cluster)
            current_cluster = [p]
        prev_pos = p
    clusters.append(current_cluster)

    for cluster in clusters:
        if len(cluster) == 1:
            continue

        k = len(cluster)
        offsets = np.linspace(-amplitude, amplitude, k)
        for p, off in zip(cluster, offsets):
            for idx in grouped[p]:
                y_map[idx] = base_y + off

    return y_map


def build_annotation_color_map(mut_top, mut_bottom):
    combined = concat_mutation_dfs(mut_top, mut_bottom)
    if combined.empty or "Annotation" not in combined.columns:
        return {}

    anns = combined["Annotation"].dropna().unique()
    anns = [a for a in anns if str(a).strip() != ""]
    if not anns:
        return {}

    palette = plt.rcParams.get("axes.prop_cycle", None)
    if palette is not None:
        colors = palette.by_key()["color"]
    else:
        colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

    color_map = {ann: colors[i % len(colors)] for i, ann in enumerate(anns)}
    return color_map


def build_score_normalizer(mut_top, mut_bottom):
    combined = concat_mutation_dfs(mut_top, mut_bottom)
    if combined.empty or "Score" not in combined.columns:
        return None, None

    vals = pd.to_numeric(combined["Score"], errors="coerce").dropna()
    if vals.empty:
        return None, None

    norm = mcolors.Normalize(vmin=float(vals.min()), vmax=float(vals.max()))
    cmap = cm.get_cmap("viridis")
    return norm, cmap


def get_point_colors(df, default_color, annotation_color_map,
                     score_norm, score_cmap, color_mode):
    if df.empty:
        return []

    colors = []
    for _, row in df.iterrows():
        if color_mode == "score":
            if score_norm is None or score_cmap is None:
                colors.append(default_color)
                continue
            s = row.get("Score", None)
            if pd.isna(s):
                colors.append(default_color)
            else:
                try:
                    val = float(s)
                    colors.append(score_cmap(score_norm(val)))
                except Exception:
                    colors.append(default_color)

        elif color_mode in ("annotation", "auto"):
            ann = row.get("Annotation", pd.NA)
            if pd.isna(ann) or str(ann).strip() == "" or not annotation_color_map:
                colors.append(default_color)
            else:
                colors.append(annotation_color_map.get(ann, default_color))

        elif color_mode == "cohort":
            colors.append(default_color)
        else:
            colors.append(default_color)

    return colors


# --- core drawing for one axis (no left panel) ---
def draw_main_axis(
    ax,
    args,
    mut_top,
    mut_bottom,
    pa,
    pt,
    protein_length,
    xlim_region,
    annotation_color_map,
    score_norm,
    score_cmap,
    color_mode,
    show_legend=True,
    show_domains=True,
    show_domain_labels=True,
    title=None,
    legend_outside=False,
):
    """
    Draw the main protein plot on a single axis.
    The first cohort (mut_top) is always plotted above the protein line.
    """
    protein_y = -2.0
    top_y = protein_y + 0.15
    bottom_y = protein_y - 0.15 if not mut_bottom.empty else None

    start_x = max(1, int(xlim_region[0]))
    end_x = min(int(protein_length), int(xlim_region[1]))
    if end_x < start_x:
        start_x, end_x = 1, int(protein_length)
    x_vals = list(range(start_x, end_x + 1))
    y_vals = [protein_y] * len(x_vals)
    ax.plot(
        x_vals,
        y_vals,
        linewidth=5,
        color=COLOR_PROTEIN_LINE,
        solid_capstyle="round",
    )

    # Title logic: priority = explicit title (facets) > --title > auto
    if title is not None:
        ax.set_title(title, fontsize=12, color=COLOR_TEXT)
    elif args.title:
        ax.set_title(args.title, fontsize=15, color=COLOR_TEXT)
    else:
        if not mut_top.empty:
            protein_id = str(mut_top.loc[0, "ProteinId"])
            gene_name = str(mut_top.loc[0, "GeneName"])
        else:
            protein_id = "NA"
            gene_name = "NA"
        ax.set_title(
            f"Amino Acid Changes in {gene_name} ({protein_id})",
            fontsize=15,
            color=COLOR_TEXT,
        )

    ax.set_xlabel("Amino Acid Position", fontsize=12, color=COLOR_TEXT)
    ax.set_ylabel("")

    if bottom_y is not None:
        y_min = bottom_y - 0.05
    else:
        y_min = protein_y - 0.25
    y_max = protein_y + 0.35
    ax.set_ylim(y_min, y_max)

    ax.set_xlim(*xlim_region)
    ax.set_yticks([])

    span = xlim_region[1] - xlim_region[0]
    max_labels = 25
    base_step = max(int(args.ticksize), 1)
    min_step_needed = span / max_labels if max_labels > 0 else base_step
    factor = max(1, math.ceil(min_step_needed / base_step))
    step = int(base_step * factor)

    start_tick = int(step * math.floor(xlim_region[0] / step))
    all_ticks = list(range(start_tick, int(protein_length) + 1, step))
    ticks = [x for x in all_ticks if xlim_region[0] <= x <= xlim_region[1]]
    ax.set_xticks(ticks)

    label_size = 9 if span > 800 else 10
    ax.tick_params(axis="x", labelrotation=90, labelsize=label_size)

    if args.grid:
        ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.3)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(COLOR_AXIS_SPINE)

    top_nonrec, top_rec = split_by_recurrence(mut_top)
    bottom_nonrec, bottom_rec = split_by_recurrence(mut_bottom)

    jw = max(args.jitter_window, 0.0)
    ja = max(args.jitter_amplitude, 0.0)

    if args.jitter == "off":
        y_top_map = {idx: top_y for idx in mut_top.index}
        y_bottom_map = {idx: bottom_y for idx in mut_bottom.index} if bottom_y is not None else {}
    else:
        y_top_map = compute_cluster_jitter(mut_top, base_y=top_y, window=jw, amplitude=ja)
        y_bottom_map = (
            compute_cluster_jitter(mut_bottom, base_y=bottom_y, window=jw, amplitude=ja)
            if bottom_y is not None else {}
        )

    top_nonrec_colors = get_point_colors(
        top_nonrec,
        default_color=COLOR_TOP_COHORT,
        annotation_color_map=annotation_color_map,
        score_norm=score_norm,
        score_cmap=score_cmap,
        color_mode=color_mode,
    )
    top_rec_colors = get_point_colors(
        top_rec,
        default_color=COLOR_TOP_COHORT,
        annotation_color_map=annotation_color_map,
        score_norm=score_norm,
        score_cmap=score_cmap,
        color_mode=color_mode,
    )
    bottom_nonrec_colors = get_point_colors(
        bottom_nonrec,
        default_color=COLOR_BOTTOM_COHORT,
        annotation_color_map=annotation_color_map,
        score_norm=score_norm,
        score_cmap=score_cmap,
        color_mode=color_mode,
    )
    bottom_rec_colors = get_point_colors(
        bottom_rec,
        default_color=COLOR_BOTTOM_COHORT,
        annotation_color_map=annotation_color_map,
        score_norm=score_norm,
        score_cmap=score_cmap,
        color_mode=color_mode,
    )

    base_point_size = max(args.point_size, 1.0)
    POINT_SIZE = base_point_size
    POINT_SIZE_REC = base_point_size * 1.3

    if not top_nonrec.empty:
        y_top_nonrec = [y_top_map.get(idx, top_y) for idx in top_nonrec.index]
        ax.scatter(
            top_nonrec["Position"],
            y_top_nonrec,
            s=POINT_SIZE,
            c=top_nonrec_colors,
            marker="o",
            edgecolor=COLOR_POINT_EDGE,
            linewidth=0.6,
            zorder=3,
        )
    if not top_rec.empty:
        y_top_rec = [y_top_map.get(idx, top_y) for idx in top_rec.index]
        ax.scatter(
            top_rec["Position"],
            y_top_rec,
            s=POINT_SIZE_REC,
            c=top_rec_colors,
            marker="D",
            edgecolor=COLOR_POINT_EDGE,
            linewidth=0.6,
            zorder=4,
        )

    if bottom_y is not None:
        if not bottom_nonrec.empty:
            y_bottom_nonrec = [y_bottom_map.get(idx, bottom_y) for idx in bottom_nonrec.index]
            ax.scatter(
                bottom_nonrec["Position"],
                y_bottom_nonrec,
                s=POINT_SIZE,
                c=bottom_nonrec_colors,
                marker="o",
                edgecolor=COLOR_POINT_EDGE,
                linewidth=0.6,
                zorder=3,
            )
        if not bottom_rec.empty:
            y_bottom_rec = [y_bottom_map.get(idx, bottom_y) for idx in bottom_rec.index]
            ax.scatter(
                bottom_rec["Position"],
                y_bottom_rec,
                s=POINT_SIZE_REC,
                c=bottom_rec_colors,
                marker="D",
                edgecolor=COLOR_POINT_EDGE,
                linewidth=0.6,
                zorder=4,
            )

    if args.showlabels == "yes":
        for _, row in mut_top.iterrows():
            pos = row["Position"]
            if pd.isna(pos):
                continue
            label = f"{row['RefAA']}{int(pos)}{row['AltAA']}"
            ax.text(
                pos,
                top_y + 0.01,
                label,
                color=COLOR_TEXT,
                fontsize=8,
                rotation=90,
                ha="center",
                va="bottom",
            )
        if bottom_y is not None:
            for _, row in mut_bottom.iterrows():
                pos = row["Position"]
                if pd.isna(pos):
                    continue
                label = f"{row['RefAA']}{int(pos)}{row['AltAA']}"
                ax.text(
                    pos,
                    bottom_y + 0.04,
                    label,
                    color=COLOR_TEXT,
                    fontsize=8,
                    rotation=90,
                    ha="center",
                    va="top",
                )

    if (not args.hide_ptms) and (pt is not None) and (not pt.empty):
        for _, row in pt.iterrows():
            site = row["site"]
            if pd.isna(site):
                continue
            ax.plot(
                [site, site],
                [protein_y, protein_y + 0.06],
                color=COLOR_PROTEIN_LINE,
                linewidth=1.4,
                zorder=2,
            )
            ax.scatter(
                site,
                protein_y + 0.06,
                s=30,
                color=COLOR_PTM,
                edgecolor=COLOR_TEXT,
                linewidth=0.6,
                alpha=0.9,
                zorder=3,
            )

    if show_domains and (not args.hide_architecture) and (pa is not None) and (not pa.empty):
        # Bars
        for _, row in pa.iterrows():
            start = row["start_site"]
            end = row["end_site"]
            if pd.isna(start) or pd.isna(end):
                continue
            width = end - start

            rect = Rectangle(
                (start, protein_y - 0.03),
                width,
                0.06,
                facecolor=COLOR_DOMAIN,
                edgecolor=COLOR_TEXT,
                linewidth=0.6,
                alpha=0.9,
                zorder=1,
            )
            ax.add_patch(rect)

        # Labels: only if both show_domain_labels AND CLI flag says yes
        if show_domain_labels and getattr(args, "architecture_labels", "no") == "yes":
            for _, row in pa.iterrows():
                start = row["start_site"]
                end = row["end_site"]
                if pd.isna(start) or pd.isna(end):
                    continue
                mid = (start + end) / 2.0
                ax.text(
                    mid,
                    protein_y + 0.05,
                    str(row["architecture_name"]),
                    fontsize=6,
                    color=COLOR_TEXT,
                    ha="center",
                    va="bottom",
                    rotation=45,
                    rotation_mode="anchor",
                    clip_on=False,
                )

    ymin, ymax = ax.get_ylim()

    def y_to_axfrac(y):
        return (y - ymin) / (ymax - ymin) if ymax != ymin else 0.5

    if not mut_top.empty:
        gene_name = str(mut_top.loc[0, "GeneName"])
        protein_id = str(mut_top.loc[0, "ProteinId"])
        protein_label = gene_name if gene_name not in ["", "nan", "None"] else protein_id
    else:
        protein_label = "Protein"

    x_label_axes = -0.03

    ax.text(
        x_label_axes,
        y_to_axfrac(top_y),
        args.mutations_name,
        color=COLOR_TOP_COHORT,
        fontsize=12,
        ha="right",
        va="center",
        transform=ax.transAxes,
        clip_on=False,
    )

    ax.text(
        x_label_axes,
        y_to_axfrac(protein_y),
        protein_label,
        color=COLOR_TEXT,
        fontsize=12,
        fontweight="bold",
        ha="right",
        va="center",
        transform=ax.transAxes,
        clip_on=False,
    )

    if bottom_y is not None:
        ax.text(
            x_label_axes,
            y_to_axfrac(bottom_y),
            args.mutations_bottom_name,
            color=COLOR_BOTTOM_COHORT,
            fontsize=12,
            ha="right",
            va="center",
            transform=ax.transAxes,
        )

    if show_legend:
        handles = []
        labels = []

        if color_mode in ("annotation", "auto") and annotation_color_map:
            for ann, col in annotation_color_map.items():
                handles.append(
                    Line2D(
                        [0], [0],
                        marker="o",
                        color="w",
                        markerfacecolor=col,
                        markersize=7,
                        linestyle="None",
                    )
                )
                labels.append(str(ann))

        if (not args.hide_architecture) and (pa is not None) and (not pa.empty):
            handles.append(Rectangle((0, 0), 1, 1, color=COLOR_DOMAIN))
            labels.append("Protein Domain")
        if (not args.hide_ptms) and (pt is not None) and (not pt.empty):
            handles.append(
                Line2D(
                    [0], [0],
                    marker="o",
                    color="w",
                    markerfacecolor=COLOR_PTM,
                    markersize=7,
                    linestyle="None",
                )
            )
            labels.append("Post-Translational Modification")

        handles.append(
            Line2D(
                [0], [0],
                marker="o",
                color=COLOR_TEXT,
                markerfacecolor=COLOR_TEXT,
                markersize=6,
                linestyle="None",
            )
        )
        labels.append("Single occurrence")

        handles.append(
            Line2D(
                [0], [0],
                marker="D",
                color=COLOR_TEXT,
                markerfacecolor=COLOR_TEXT,
                markersize=6,
                linestyle="None",
            )
        )
        labels.append("Recurrent (≥2)")

        if handles:
            legend_kwargs = dict(
                loc="upper right",
                frameon=True,
                framealpha=1.0,
                facecolor="white" if args.theme == "light" else "#222222",
                edgecolor="white" if args.theme == "light" else "#444444",
                fontsize=9,
            )
            if legend_outside:
                legend_kwargs["loc"] = "upper left"
                legend_kwargs["bbox_to_anchor"] = (1.02, 1.0)
                legend_kwargs["borderaxespad"] = 0.0

            leg = ax.legend(handles, labels, **legend_kwargs)
            for txt in leg.get_texts():
                txt.set_color(COLOR_TEXT)


# --- full-protein plot (single panel) ---
def plot_protein_full(
    args,
    mut_top,
    mut_bottom,
    pa,
    pt,
    protein_length,
    xlim_region,
    output_path,
    fmt,
    annotation_color_map,
    score_norm,
    score_cmap,
    color_mode,
):
    # Use constrained_layout to handle spacing (no tight_layout)
    fig, ax_main = plt.subplots(1, 1, figsize=(10, 7.5), constrained_layout=True)

    draw_main_axis(
        ax=ax_main,
        args=args,
        mut_top=mut_top,
        mut_bottom=mut_bottom,
        pa=pa,
        pt=pt,
        protein_length=protein_length,
        xlim_region=xlim_region,
        annotation_color_map=annotation_color_map,
        score_norm=score_norm,
        score_cmap=score_cmap,
        color_mode=color_mode,
        show_legend=True,
        show_domains=True,
        show_domain_labels=True,  # labels further gated by --architecture-labels
        title=None,
        legend_outside=False,
    )

    # Normal-sized colorbar for full plot (if score mode)
    if color_mode == "score" and score_norm is not None and score_cmap is not None:
        sm = cm.ScalarMappable(norm=score_norm, cmap=score_cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_main, label="Score")
        cbar.ax.yaxis.label.set_color(COLOR_TEXT)
        for t in cbar.ax.get_yticklabels():
            t.set_color(COLOR_TEXT)

    save_kwargs = {}
    if fmt == "png":
        save_kwargs["dpi"] = args.dpi
    fig.savefig(output_path, format=fmt, **save_kwargs)
    plt.close(fig)


# --- faceted-by-domain plot ---
def plot_protein_faceted(
    args,
    mut_top,
    mut_bottom,
    pa,
    pt,
    protein_length,
    output_path,
    fmt,
    annotation_color_map,
    score_norm,
    score_cmap,
    color_mode,
):
    n_panels = len(pa)
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=(10, max(2.5 * n_panels, 5)),
        sharey=True,
        constrained_layout=True,  # let matplotlib handle layout
    )

    if n_panels == 1:
        axes = [axes]

    if not mut_top.empty:
        protein_id = str(mut_top.loc[0, "ProteinId"])
        gene_name = str(mut_top.loc[0, "GeneName"])
    else:
        protein_id = "NA"
        gene_name = "NA"
    fig.suptitle(
        f"{args.name} – Faceted by domain ({gene_name} {protein_id})",
        fontsize=14,
        color=COLOR_TEXT,
    )

    for i, (_, row) in enumerate(pa.iterrows()):
        ax = axes[i]
        start = row["start_site"]
        end = row["end_site"]
        if pd.isna(start) or pd.isna(end):
            continue

        span = max(end - start, 1.0)
        pad_left = max(8.0, span * 0.15)
        pad_right = max(4.0, span * 0.08)
        xlim_region = (
            max(0.0, start - pad_left),
            min(protein_length, end + pad_right),
        )

        pa_sub = pa.iloc[[i]]
        title = f"{row['architecture_name']} ({int(start)}–{int(end)})"
        draw_main_axis(
            ax=ax,
            args=args,
            mut_top=mut_top,
            mut_bottom=mut_bottom,
            pa=pa_sub,
            pt=pt,
            protein_length=protein_length,
            xlim_region=xlim_region,
            annotation_color_map=annotation_color_map,
            score_norm=score_norm,
            score_cmap=score_cmap,
            color_mode=color_mode,
            show_legend=(i == n_panels - 1),
            show_domains=True,
            show_domain_labels=False,  # no inline labels in facets
            title=title,
            legend_outside=True,
        )

    # SUPER TINY colorbar ONLY for the faceted plot (score mode)
    if color_mode == "score" and score_norm is not None and score_cmap is not None:
        sm = cm.ScalarMappable(norm=score_norm, cmap=score_cmap)
        sm.set_array([])

        cbar = fig.colorbar(
            sm,
            ax=axes,
            orientation="horizontal",  # bottom, across panels
            label="Score",
            fraction=0.02,   # thickness of the bar (very thin)
            pad=0.06,        # distance from bottom subplot
            shrink=0.4,      # make it shorter than total width
            aspect=40,
        )
        cbar.ax.xaxis.label.set_color(COLOR_TEXT)
        cbar.ax.xaxis.label.set_size(8)
        cbar.ax.tick_params(labelsize=6)
        for t in cbar.ax.get_xticklabels():
            t.set_color(COLOR_TEXT)

    save_kwargs = {}
    if fmt == "png":
        save_kwargs["dpi"] = args.dpi
    fig.savefig(output_path, format=fmt, **save_kwargs)
    plt.close(fig)


# --- main ---
def main(argv=None):
    args = parse_args(argv)
    configure_style(args)
    mut_top, mut_bottom, pa, pt = read_inputs(args)

    # Apply filters
    mut_top = filter_mutations(mut_top, args)
    mut_bottom = filter_mutations(mut_bottom, args)

    if mut_top.empty and mut_bottom.empty:
        sys.stderr.write(
            "Note: no mutations to plot (after filtering); plotting protein architecture/PTMs only.\n"
        )

    protein_length, xlim_region = compute_ranges(args, mut_top, mut_bottom)
    output_path, fmt = get_output_path_and_format(args, mut_top)

    if pa is not None and not pa.empty:
        bad_pa = pa[
            (pa["start_site"] < 1) |
            (pa["end_site"] > protein_length)
        ]
        if not bad_pa.empty:
            sys.stderr.write(
                "Warning: some architecture regions fall outside [1, "
                f"{int(protein_length)}]. They will still be plotted but may look odd.\n"
            )

    if pt is not None and not pt.empty:
        bad_pt = pt[
            (pt["site"] < 1) |
            (pt["site"] > protein_length)
        ]
        if not bad_pt.empty:
            sys.stderr.write(
                "Warning: some PTM sites fall outside [1, "
                f"{int(protein_length)}].\n"
            )

    annotation_color_map = build_annotation_color_map(mut_top, mut_bottom)

    score_norm = score_cmap = None
    color_mode = args.color_by

    if color_mode == "score":
        score_norm, score_cmap = build_score_normalizer(mut_top, mut_bottom)
        if score_norm is None or score_cmap is None:
            sys.stderr.write(
                "Warning: --color-by score requested, but no valid numeric Score values found. "
                "Falling back to cohort-based colors.\n"
            )
            color_mode = "cohort"

    if color_mode == "annotation" and not annotation_color_map:
        sys.stderr.write(
            "Warning: --color-by annotation requested, but no Annotation values found. "
            "Falling back to cohort-based colors.\n"
        )
        color_mode = "cohort"

    if color_mode == "auto":
        if annotation_color_map:
            color_mode = "annotation"
        else:
            color_mode = "cohort"

    # Optionally write annotation color map (only in annotation mode)
    if args.annotation_colors_out:
        if color_mode == "annotation":
            if annotation_color_map:
                write_annotation_color_map(args.annotation_colors_out, annotation_color_map)
            else:
                sys.stderr.write(
                    "Note: --annotation-colors-out was given, but no Annotation values were found; "
                    "writing a header-only annotation color file.\n"
                )
                write_annotation_color_map(args.annotation_colors_out, {})
        else:
            sys.stderr.write(
                f"Note: --annotation-colors-out is only used when color mode is 'annotation' "
                f"(or --color-by auto resolves to annotation). Current mode is '{color_mode}', "
                "so no annotation color file will be written.\n"
            )

    # Optionally write score colormap bins (only in score mode)
    if args.score_colors_out:
        if color_mode == "score" and score_norm is not None and score_cmap is not None:
            write_score_colormap(args.score_colors_out, score_norm, score_cmap)
        else:
            sys.stderr.write(
                "Note: --score-colors-out is only used when color mode is 'score' with valid scores; "
                f"current mode is '{color_mode}', so no score colormap file will be written.\n"
            )

    if args.facet_domains:
        if pa is None or pa.empty:
            sys.exit(
                "Error: --facet-domains was requested, but no architecture data is available.\n"
                "Provide an architecture file with columns: architecture_name, start_site, end_site."
            )

        plot_protein_faceted(
            args,
            mut_top,
            mut_bottom,
            pa,
            pt,
            protein_length,
            output_path,
            fmt,
            annotation_color_map,
            score_norm,
            score_cmap,
            color_mode,
        )
    else:
        plot_protein_full(
            args,
            mut_top,
            mut_bottom,
            pa,
            pt,
            protein_length,
            xlim_region,
            output_path,
            fmt,
            annotation_color_map,
            score_norm,
            score_cmap,
            color_mode,
        )

    print(
        f"Protein plot written to: {output_path} "
        f"(format={fmt}, theme={args.theme}, palette={args.palette}, dpi={args.dpi})"
    )


if __name__ == "__main__":
    main()

