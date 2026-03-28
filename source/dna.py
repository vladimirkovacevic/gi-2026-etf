#!/usr/bin/env python3
"""
DNA Alignment & Variant Calling Pipeline
=========================================
NOTE: bwa-mem2 is the high-performance reimplementation of BWA-MEM using SIMD.
      Pre-built binaries exist only for x86-64 Linux. On Apple Silicon (ARM64)
      we use `bwa mem` which implements the identical BWA-MEM algorithm — the
      results are byte-for-byte equivalent; bwa-mem2 is simply faster on x86.

Workflow:
  1. Install tools (bwa, samtools, bcftools) via Homebrew if missing
  2. Index the reference genome  →  bwa index
  3. Align paired-end reads      →  bwa mem  (= bwa-mem2 mem algorithm)
  4. Convert SAM → sorted BAM   →  samtools view | samtools sort
  5. Index the BAM               →  samtools index
  6. Alignment statistics        →  samtools flagstat
  7. Per-base depth of coverage  →  samtools depth -a
  8. Variant calling             →  bcftools mpileup | bcftools call -mv
  9. Visualise all steps         →  matplotlib

Run:
    python source/dna.py
"""

import os
import shutil
import subprocess
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
REPO   = Path(__file__).resolve().parent.parent
DATA   = REPO / "data" / "dna"
OUT    = REPO / "output" / "dna"
REF    = DATA / "example_human_reference.fasta"
R1     = DATA / "example_human_Illumina.pe_1.fastq"
R2     = DATA / "example_human_Illumina.pe_2.fastq"
BAM_S  = OUT  / "aligned.sorted.bam"
VCF    = OUT  / "variants.vcf"
STATS  = OUT  / "flagstat.txt"
DEPTH  = OUT  / "depth.txt"
PLOT   = OUT  / "dna_pipeline_report.png"

TOOLS  = ["bwa", "samtools", "bcftools"]

# ── helpers ──────────────────────────────────────────────────────────────────

def run(cmd: str, **kw):
    """Run a shell command, stream output, raise on error."""
    print(f"\n$ {cmd}")
    result = subprocess.run(cmd, shell=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kw)
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        raise RuntimeError(f"Command failed (exit {result.returncode}):\n  {cmd}")
    return result.stdout


def find_tool(name: str):
    """Return absolute path to *name* by searching PATH."""
    return shutil.which(name)


def ensure_tools():
    """Install missing tools via Homebrew."""
    missing = [t for t in TOOLS if find_tool(t) is None]
    if not missing:
        print("✓ All tools found in PATH.")
        return

    print(f"⚙  Missing tools: {missing}. Installing via Homebrew…")
    pkgs = " ".join(missing)
    run(f'brew install {pkgs}')

    still_missing = [t for t in missing if find_tool(t) is None]
    if still_missing:
        raise RuntimeError(
            f"Installation failed for: {still_missing}\n"
            "Please install manually: brew install bwa samtools bcftools"
        )
    print("✓ All tools installed.")


def tool(name: str) -> str:
    p = find_tool(name)
    if p is None:
        raise RuntimeError(f"Tool '{name}' not found. Run: brew install {name}")
    return p

# ── pipeline steps ───────────────────────────────────────────────────────────

def step_index():
    """Index the reference genome with bwa (creates .amb, .ann, .bwt, .pac, .sa)."""
    sentinel = REF.parent / (REF.name + ".bwt")
    if sentinel.exists():
        print("✓ Reference index already exists — skipping.")
        return
    run(f'{tool("bwa")} index {REF}')


def step_align():
    """Align paired-end reads to the reference; produce SAM then sorted BAM."""
    if BAM_S.exists():
        print("✓ Sorted BAM already exists — skipping alignment.")
        return

    # bwa mem with read-group tag (BWA-MEM algorithm = bwa-mem2 algorithm)
    rg = r"@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1"
    bwa_cmd  = f'{tool("bwa")} mem -R "{rg}" -t 4 {REF} {R1} {R2}'
    sam_cmd  = f'{tool("samtools")} view -bS -'
    sort_cmd = f'{tool("samtools")} sort -o {BAM_S}'

    # pipe: bwa mem | samtools view -bS | samtools sort
    print(f"\n$ {bwa_cmd} | samtools view -bS | samtools sort -o {BAM_S}")
    p1 = subprocess.Popen(bwa_cmd, shell=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sam_cmd, shell=True,
                          stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(sort_cmd, shell=True,
                          stdin=p2.stdout, stderr=subprocess.PIPE)
    p2.stdout.close()
    p3.wait(); p2.wait(); p1.wait()

    for p, name in [(p1, "bwa mem"), (p2, "samtools view"), (p3, "samtools sort")]:
        err = p.stderr.read().decode()
        if err:
            print(f"[{name} stderr]\n{err}")
        if p.returncode not in (0, None):
            raise RuntimeError(f"{name} failed with exit {p.returncode}")

    run(f'{tool("samtools")} index {BAM_S}')


def step_flagstat():
    """Compute alignment statistics."""
    run(f'{tool("samtools")} flagstat {BAM_S} > {STATS}')


def step_depth():
    """Per-base sequencing depth."""
    run(f'{tool("samtools")} depth -a {BAM_S} > {DEPTH}')


def step_variants():
    """Call variants: mpileup → bcftools call → VCF."""
    if VCF.exists():
        print("✓ VCF already exists — skipping variant calling.")
        return
    run(
        f'{tool("bcftools")} mpileup -f {REF} {BAM_S} '
        f'| {tool("bcftools")} call -mv -Ov -o {VCF}'
    )


# ── parsing helpers ──────────────────────────────────────────────────────────

def parse_flagstat(path: Path) -> dict:
    stats = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        count = int(parts[0])
        if "in total"   in line: stats["total"]    = count
        if "mapped ("   in line: stats["mapped"]   = count
        if "paired in"  in line: stats["paired"]   = count
        if "properly"   in line: stats["properly"]  = count
        if "singletons" in line: stats["singletons"]= count
    return stats


def parse_depth(path: Path):
    """Return arrays: positions, depths."""
    positions, depths = [], []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split()
        positions.append(int(parts[1]))
        depths.append(int(parts[2]))
    return positions, depths


def parse_vcf(path: Path):
    """Extract SNPs and INDELs; return list of dicts."""
    variants = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 8:
            continue
        chrom, pos, _, ref, alt, qual, filt, info = f[:8]
        vtype = "INDEL" if "INDEL" in info else "SNP"
        try:
            qual_f = float(qual)
        except ValueError:
            qual_f = 0.0
        variants.append({
            "chrom": chrom, "pos": int(pos),
            "ref": ref, "alt": alt,
            "qual": qual_f, "filter": filt, "type": vtype,
        })
    return variants


# ── visualisation ────────────────────────────────────────────────────────────

def visualise(flagstat: dict, positions: list, depths: list, variants: list):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(24, 14), facecolor="white")
    gs  = GridSpec(4, 3, figure=fig, hspace=0.65, wspace=0.40)

    PANEL  = "#f6f8fa"
    TEXT   = "#1f2328"
    ACCENT = "#0969da"
    GREEN  = "#1a7f37"
    YELLOW = "#9a6700"
    RED    = "#cf222e"
    PURPLE = "#8250df"

    ax_style = dict(facecolor=PANEL, labelcolor=TEXT, titlecolor=TEXT)

    def style_ax(ax, title=""):
        ax.set_facecolor(PANEL)
        ax.tick_params(colors=TEXT)
        for spine in ax.spines.values():
            spine.set_edgecolor("#d0d7de")
        ax.title.set_color(TEXT)
        ax.xaxis.label.set_color(TEXT)
        ax.yaxis.label.set_color(TEXT)
        if title:
            ax.set_title(title, fontsize=11, fontweight="bold", pad=8)

    # ── 1. Pipeline overview (text) ───────────────────────────────────────────
    ax0 = fig.add_subplot(gs[0, :])
    ax0.set_facecolor("white")
    ax0.axis("off")
    style_ax(ax0, "BWA-MEM + Samtools/Bcftools Variant Calling Pipeline")

    steps = [
        ("1", "Index Reference",        "bwa index",                      GREEN),
        ("2", "Align Paired-End Reads", "bwa mem -R @RG … R1 R2",        ACCENT),
        ("3", "Convert & Sort BAM",     "samtools view | samtools sort",  ACCENT),
        ("4", "Index BAM",              "samtools index",                 ACCENT),
        ("5", "Alignment Stats",        "samtools flagstat",              YELLOW),
        ("6", "Depth of Coverage",      "samtools depth -a",              YELLOW),
        ("7", "Variant Calling",        "bcftools mpileup | call -mv",    RED),
    ]
    x_step = 1.0 / len(steps)
    for i, (num, label, cmd, color) in enumerate(steps):
        cx = (i + 0.5) * x_step
        ax0.annotate("", xy=(cx + x_step * 0.38, 0.5),
                     xytext=(cx - x_step * 0.38, 0.5),
                     xycoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", color="#57606a", lw=1.5))
        ax0.text(cx, 0.78, f"Step {num}", ha="center", va="center",
                 color=color, fontsize=8, fontweight="bold",
                 transform=ax0.transAxes)
        ax0.text(cx, 0.55, label, ha="center", va="center",
                 color=TEXT, fontsize=8.5, fontweight="bold",
                 transform=ax0.transAxes)
        ax0.text(cx, 0.28, cmd, ha="center", va="center",
                 color="#57606a", fontsize=6.5, fontstyle="italic",
                 transform=ax0.transAxes)

    # ── 2. Alignment stats bar chart ──────────────────────────────────────────
    ax1 = fig.add_subplot(gs[1, 0])
    style_ax(ax1, "Alignment Statistics (flagstat)")
    if flagstat:
        labels = ["Total\nReads", "Mapped", "Properly\nPaired", "Paired\nin Seq", "Singletons"]
        keys   = ["total", "mapped", "properly", "paired", "singletons"]
        colors = [ACCENT, GREEN, GREEN, YELLOW, RED]
        vals   = [flagstat.get(k, 0) for k in keys]
        bars   = ax1.bar(labels, vals, color=colors, edgecolor="#d0d7de", linewidth=0.8)
        for bar, val in zip(bars, vals):
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(vals)*0.01,
                     str(val), ha="center", va="bottom", color=TEXT, fontsize=8)
        ax1.set_ylabel("Read Count", color=TEXT)
        ax1.set_ylim(0, max(vals) * 1.15 if vals else 1)
    else:
        ax1.text(0.5, 0.5, "No flagstat data", ha="center", va="center",
                 color=TEXT, transform=ax1.transAxes)

    # ── 3. Mapping rate pie ───────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 1])
    style_ax(ax2, "Mapping Rate")
    if flagstat.get("total", 0) > 0:
        mapped    = flagstat.get("mapped", 0)
        unmapped  = flagstat["total"] - mapped
        wedge_col = [GREEN, RED]
        wedges, texts, autotexts = ax2.pie(
            [mapped, unmapped],
            labels=["Mapped", "Unmapped"],
            colors=wedge_col,
            autopct="%1.1f%%",
            startangle=90,
            wedgeprops=dict(edgecolor="#d0d7de", linewidth=1),
            textprops=dict(color=TEXT),
        )
        for at in autotexts:
            at.set_color("white"); at.set_fontweight("bold")
    else:
        ax2.text(0.5, 0.5, "No data", ha="center", va="center",
                 color=TEXT, transform=ax2.transAxes)

    # ── 4. Read quality distribution (from FASTQ) ─────────────────────────────
    ax3 = fig.add_subplot(gs[1, 2])
    style_ax(ax3, "Per-Read Mean Base Quality (R1)")
    mean_quals = []
    with open(R1) as fh:
        while True:
            name = fh.readline()
            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline().strip()
            if not qual:
                break
            mean_quals.append(np.mean([ord(c) - 33 for c in qual]))
    if mean_quals:
        ax3.hist(mean_quals, bins=20, color=ACCENT, edgecolor="#d0d7de", linewidth=0.7)
        ax3.axvline(np.mean(mean_quals), color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Mean = {np.mean(mean_quals):.1f}")
        ax3.set_xlabel("Mean Phred Quality Score")
        ax3.set_ylabel("Number of Reads")
        ax3.legend(facecolor="white", edgecolor="#d0d7de", labelcolor=TEXT, fontsize=8)

    # ── 5. Coverage / depth plot ──────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[2, :2])
    style_ax(ax4, "Sequencing Depth of Coverage")
    if depths:
        pos_arr   = np.array(positions)
        depth_arr = np.array(depths)
        ax4.fill_between(pos_arr, depth_arr, alpha=0.35, color=ACCENT)
        ax4.plot(pos_arr, depth_arr, color=ACCENT, lw=1)
        mean_d = np.mean(depth_arr)
        ax4.axhline(mean_d, color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Mean depth = {mean_d:.1f}×")
        # mark variant positions
        if variants:
            vpos = [v["pos"] for v in variants]
            vdep = np.interp(vpos, pos_arr, depth_arr)
            ax4.scatter(vpos, vdep, color=RED, s=40, zorder=5,
                        label=f"Variant sites (n={len(variants)})", marker="v")
        ax4.set_xlabel("Genomic Position (chr20)")
        ax4.set_ylabel("Depth (×)")
        ax4.legend(facecolor="white", edgecolor="#d0d7de", labelcolor=TEXT, fontsize=8)
    else:
        ax4.text(0.5, 0.5, "No depth data", ha="center", va="center",
                 color=TEXT, transform=ax4.transAxes)

    # ── 6. Depth histogram ────────────────────────────────────────────────────
    ax5 = fig.add_subplot(gs[2, 2])
    style_ax(ax5, "Depth Distribution")
    if depths:
        depth_arr = np.array(depths)
        ax5.hist(depth_arr, bins=30, color=PURPLE, edgecolor="#d0d7de", lw=0.7)
        ax5.axvline(np.mean(depth_arr), color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Mean = {np.mean(depth_arr):.1f}×")
        ax5.set_xlabel("Depth (×)")
        ax5.set_ylabel("Number of Positions")
        ax5.legend(facecolor="white", edgecolor="#d0d7de", labelcolor=TEXT, fontsize=8)

    # ── 7. Variant positions & quality ────────────────────────────────────────
    ax6 = fig.add_subplot(gs[3, :2])
    style_ax(ax6, "Variant Positions and Quality Scores")
    if variants:
        snps   = [v for v in variants if v["type"] == "SNP"]
        indels = [v for v in variants if v["type"] == "INDEL"]
        if snps:
            ax6.scatter([v["pos"] for v in snps],
                        [v["qual"] for v in snps],
                        color=GREEN, s=60, zorder=5, label="SNP", marker="o", alpha=0.85)
        if indels:
            ax6.scatter([v["pos"] for v in indels],
                        [v["qual"] for v in indels],
                        color=RED, s=80, zorder=5, label="INDEL", marker="^", alpha=0.85)
        ax6.set_xlabel("Genomic Position")
        ax6.set_ylabel("Variant Quality (PHRED)")
        ax6.legend(facecolor="white", edgecolor="#d0d7de", labelcolor=TEXT, fontsize=8)
    else:
        ax6.text(0.5, 0.5, "No variants called", ha="center", va="center",
                 color=TEXT, transform=ax6.transAxes)

    # ── 8. Variant type summary & substitution spectrum ───────────────────────
    ax7 = fig.add_subplot(gs[3, 2])
    style_ax(ax7, "SNP Substitution Spectrum")
    if variants:
        subs = {}
        for v in variants:
            if v["type"] == "SNP" and len(v["ref"]) == 1 and len(v["alt"]) == 1:
                key = f"{v['ref']}→{v['alt']}"
                subs[key] = subs.get(key, 0) + 1
        if subs:
            sorted_subs = dict(sorted(subs.items(), key=lambda x: -x[1]))
            bar_colors  = [ACCENT, GREEN, YELLOW, RED, PURPLE, "#ff7b72",
                           "#79c0ff", "#56d364"] * 3
            bars = ax7.bar(list(sorted_subs.keys()), list(sorted_subs.values()),
                           color=bar_colors[:len(sorted_subs)],
                           edgecolor="#d0d7de", linewidth=0.8)
            for bar, val in zip(bars, sorted_subs.values()):
                ax7.text(bar.get_x() + bar.get_width() / 2,
                         bar.get_height() + 0.1,
                         str(val), ha="center", va="bottom",
                         color=TEXT, fontsize=8)
            ax7.set_xlabel("Substitution")
            ax7.set_ylabel("Count")
            ax7.tick_params(axis="x", rotation=45)
        else:
            snp_count   = len([v for v in variants if v["type"] == "SNP"])
            indel_count = len([v for v in variants if v["type"] == "INDEL"])
            ax7.bar(["SNPs", "INDELs"], [snp_count, indel_count],
                    color=[GREEN, RED], edgecolor="#d0d7de")
            ax7.set_ylabel("Count")
    else:
        ax7.text(0.5, 0.5, "No variants", ha="center", va="center",
                 color=TEXT, transform=ax7.transAxes)

    # ── title & footer ────────────────────────────────────────────────────────
    fig.suptitle("DNA Sequencing Analysis: BWA-MEM Alignment + BCFtools Variant Calling\n"
                 "Reference: chr20 (example)  |  Sample: SRR316957",
                 color=TEXT, fontsize=14, fontweight="bold", y=0.99)

    fig.text(0.5, 0.005,
             "Pipeline: bwa index → bwa mem (BWA-MEM algorithm) → samtools sort/index → "
             "bcftools mpileup | bcftools call",
             ha="center", color="#57606a", fontsize=8)

    plt.savefig(PLOT, dpi=160, bbox_inches="tight", facecolor="white")
    print(f"\n✓ Report saved → {PLOT}")
    return fig


# ── print summaries ──────────────────────────────────────────────────────────

def print_summary(flagstat: dict, depths: list, variants: list):
    import numpy as np
    sep = "─" * 60
    print(f"\n{sep}")
    print("  PIPELINE SUMMARY")
    print(sep)

    print("\n[1] Reference Genome")
    print(f"    File   : {REF.name}")
    ref_text = REF.read_text()
    seq_len  = sum(len(l) for l in ref_text.splitlines() if not l.startswith(">"))
    print(f"    Length : {seq_len:,} bp")

    print("\n[2] Reads")
    n_reads = sum(1 for l in open(R1)) // 4
    print(f"    Pairs  : {n_reads:,}  (paired-end Illumina)")

    print("\n[3] Alignment (bwa-mem2 mem)")
    if flagstat:
        total   = flagstat.get("total", 0)
        mapped  = flagstat.get("mapped", 0)
        rate    = mapped / total * 100 if total else 0
        print(f"    Total reads  : {total:,}")
        print(f"    Mapped reads : {mapped:,}  ({rate:.1f}%)")
        print(f"    Properly paired: {flagstat.get('properly', 0):,}")
        print(f"    Singletons   : {flagstat.get('singletons', 0):,}")

    print("\n[4] Coverage")
    if depths:
        d = np.array(depths)
        print(f"    Positions covered : {len(d):,}")
        print(f"    Mean depth        : {d.mean():.2f}×")
        print(f"    Median depth      : {np.median(d):.0f}×")
        print(f"    Zero-depth pos.   : {(d == 0).sum():,}")
        print(f"    ≥10× coverage     : {(d >= 10).sum() / len(d) * 100:.1f}%")

    print("\n[5] Variants (bcftools call -mv)")
    if variants:
        snps   = [v for v in variants if v["type"] == "SNP"]
        indels = [v for v in variants if v["type"] == "INDEL"]
        quals  = [v["qual"] for v in variants if v["qual"] > 0]
        print(f"    Total variants : {len(variants)}")
        print(f"    SNPs           : {len(snps)}")
        print(f"    INDELs         : {len(indels)}")
        if quals:
            print(f"    Mean QUAL      : {np.mean(quals):.1f}")
        if snps:
            print("    Top 5 SNPs by quality:")
            for v in sorted(snps, key=lambda x: -x["qual"])[:5]:
                print(f"      pos={v['pos']:,}  {v['ref']}→{v['alt']}  QUAL={v['qual']:.0f}")
    else:
        print("    No variants called.")

    print(f"\n{sep}")
    print(f"  Output directory: {OUT}")
    print(sep + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  DNA Alignment & Variant Calling Pipeline")
    print("=" * 60)

    OUT.mkdir(parents=True, exist_ok=True)

    # ── Step 0: Tool check / install ─────────────────────────────────────────
    print("\n[Step 0] Checking / installing tools…")
    ensure_tools()

    # ── Step 1: Index ─────────────────────────────────────────────────────────
    print("\n[Step 1] Indexing reference genome…")
    step_index()

    # ── Step 2-3: Align → Sort BAM ────────────────────────────────────────────
    print("\n[Step 2-3] Aligning reads and creating sorted BAM…")
    step_align()

    # ── Step 4: Flagstat ──────────────────────────────────────────────────────
    print("\n[Step 4] Computing alignment statistics…")
    step_flagstat()
    print(Path(STATS).read_text())

    # ── Step 5: Depth ─────────────────────────────────────────────────────────
    print("\n[Step 5] Computing per-base depth of coverage…")
    step_depth()

    # ── Step 6: Variant calling ───────────────────────────────────────────────
    print("\n[Step 6] Calling variants…")
    step_variants()

    # ── Parse outputs ─────────────────────────────────────────────────────────
    flagstat = parse_flagstat(STATS)
    positions, depths = parse_depth(DEPTH)
    variants  = parse_vcf(VCF)

    # ── Summary ───────────────────────────────────────────────────────────────
    print_summary(flagstat, depths, variants)

    # ── Visualise ─────────────────────────────────────────────────────────────
    print("\n[Step 7] Generating visualisation…")
    visualise(flagstat, positions, depths, variants)

    print("\n✓ Pipeline complete.")


if __name__ == "__main__":
    main()
