#!/usr/bin/env python3
"""
BWA-MEM2 Alignment Pipeline & Visualisation
============================================
NOTE: bwa-mem2 pre-built binaries target x86-64 Linux.  On Apple Silicon
      (ARM64 / macOS) this script automatically falls back to `bwa mem`,
      which implements the identical BWA-MEM algorithm — results are
      algorithmically equivalent; bwa-mem2 is simply faster on x86.

Workflow:
  1. Check / install tools  (bwa-mem2 or bwa, samtools)
  2. Index reference genome  →  bwa-mem2 index  (or bwa index)
  3. Align paired-end reads  →  bwa-mem2 mem     (or bwa mem)
  4. SAM → sorted BAM        →  samtools view | samtools sort
  5. Index BAM               →  samtools index
  6. Alignment statistics    →  samtools flagstat
  7. Per-base depth          →  samtools depth -a
  8. Per-position pileup     →  samtools mpileup  (base composition)
  9. Visualise               →  matplotlib

Run:
    python source/aln.py
"""

import os
import re
import shutil
import subprocess
import platform
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
REPO   = Path(__file__).resolve().parent.parent
DATA   = REPO / "data" / "dna"
OUT    = REPO / "output" / "aln"
REF    = DATA / "example_human_reference.fasta"
R1     = DATA / "example_human_Illumina.pe_1.fastq"
R2     = DATA / "example_human_Illumina.pe_2.fastq"
BAM_S  = OUT  / "aligned.sorted.bam"
STATS  = OUT  / "flagstat.txt"
DEPTH  = OUT  / "depth.txt"
PILE   = OUT  / "pileup.txt"
PLOT   = OUT  / "aln_report.png"

# ── helpers ───────────────────────────────────────────────────────────────────

def run(cmd: str, **kw) -> str:
    """Run a shell command, stream output, raise on error."""
    print(f"\n$ {cmd}")
    result = subprocess.run(
        cmd, shell=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kw
    )
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        raise RuntimeError(f"Command failed (exit {result.returncode}):\n  {cmd}")
    return result.stdout


def tool(name: str) -> str:
    p = shutil.which(name)
    if p is None:
        raise RuntimeError(f"Tool '{name}' not found.  Run: brew install {name}")
    return p


def pick_aligner() -> tuple[str, str]:
    """
    Return (aligner_bin, aligner_label).
    Prefer bwa-mem2; fall back to bwa on platforms without bwa-mem2 binaries.
    """
    bwa2 = shutil.which("bwa-mem2")
    if bwa2:
        return bwa2, "bwa-mem2"

    machine = platform.machine().lower()
    if machine in ("arm64", "aarch64"):
        print(
            "ℹ  bwa-mem2 not found.  On Apple Silicon (ARM64) no pre-built binary\n"
            "   exists; falling back to `bwa mem` (identical BWA-MEM algorithm)."
        )
    else:
        print("ℹ  bwa-mem2 not found in PATH; falling back to `bwa mem`.")

    bwa = shutil.which("bwa")
    if bwa is None:
        raise RuntimeError(
            "Neither bwa-mem2 nor bwa found.  Install one:\n"
            "  conda install -c bioconda bwa-mem2   # Linux x86-64\n"
            "  brew install bwa                     # macOS"
        )
    return bwa, "bwa (BWA-MEM algorithm)"


def ensure_samtools():
    if shutil.which("samtools") is None:
        print("⚙  samtools not found — installing via Homebrew…")
        run("brew install samtools")
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools install failed.  Run: brew install samtools")
    print("✓ samtools found.")


# ── pipeline steps ────────────────────────────────────────────────────────────

def step_index(aligner: str, label: str):
    """Index the reference genome."""
    # bwa-mem2 creates .0123 / .amb / .ann / .pac / .bwt.2bit.64 etc.
    # bwa creates .amb / .ann / .bwt / .pac / .sa
    sentinel_bwa2 = REF.parent / (REF.name + ".0123")
    sentinel_bwa  = REF.parent / (REF.name + ".bwt")

    if sentinel_bwa2.exists() or sentinel_bwa.exists():
        print(f"✓ Reference index already exists — skipping ({label} index).")
        return

    run(f'"{aligner}" index {REF}')


def step_align(aligner: str, label: str):
    """Align paired-end reads; produce sorted, indexed BAM."""
    if BAM_S.exists():
        print("✓ Sorted BAM already exists — skipping alignment.")
        return

    rg = r"@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1"
    aln_cmd  = f'"{aligner}" mem -R "{rg}" -t 4 {REF} {R1} {R2}'
    view_cmd = f'{tool("samtools")} view -bS -'
    sort_cmd = f'{tool("samtools")} sort -o {BAM_S}'

    print(f"\n$ {aln_cmd} | samtools view -bS | samtools sort -o {BAM_S}")
    p1 = subprocess.Popen(aln_cmd,  shell=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(view_cmd, shell=True,
                          stdin=p1.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(sort_cmd, shell=True,
                          stdin=p2.stdout, stderr=subprocess.PIPE)
    p2.stdout.close()
    p3.wait(); p2.wait(); p1.wait()

    for p, name in [(p1, label), (p2, "samtools view"), (p3, "samtools sort")]:
        err = p.stderr.read().decode()
        if err:
            print(f"[{name} stderr]\n{err}")
        if p.returncode not in (0, None):
            raise RuntimeError(f"{name} failed (exit {p.returncode})")

    run(f'{tool("samtools")} index {BAM_S}')


def step_flagstat():
    run(f'{tool("samtools")} flagstat {BAM_S} > {STATS}')


def step_depth():
    run(f'{tool("samtools")} depth -a {BAM_S} > {DEPTH}')


def step_pileup():
    """Per-position base composition (A/T/G/C counts via mpileup)."""
    run(
        f'{tool("samtools")} mpileup --no-BAQ -d 10000 -f {REF} {BAM_S} '
        f'> {PILE}'
    )


# ── parse samtools view for per-read metrics ──────────────────────────────────

def parse_read_metrics(max_reads: int = 200_000) -> dict:
    """
    Stream samtools view output and collect:
      - insert_sizes  : list[int]   (|TLEN| for properly-paired reads)
      - mapqs         : list[int]   (MAPQ for all primary alignments)
      - soft_clips    : list[int]   (total soft-clip bp per read)
      - read_lengths  : list[int]
      - mismatches    : list[int]   (NM tag value per read)
    """
    cmd = f'{tool("samtools")} view -F 0x900 {BAM_S}'   # primary alignments only
    proc = subprocess.Popen(cmd, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)

    insert_sizes, mapqs, soft_clips, read_lengths, mismatches = [], [], [], [], []
    count = 0

    for line in proc.stdout:
        if count >= max_reads:
            break
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        flag  = int(fields[1])
        mapq  = int(fields[4])
        cigar = fields[5]
        tlen  = int(fields[8])
        seq   = fields[9]

        mapqs.append(mapq)
        read_lengths.append(len(seq))

        # soft clips from CIGAR
        sc = sum(int(n) for n, op in re.findall(r'(\d+)([A-Z])', cigar) if op == 'S')
        soft_clips.append(sc)

        # insert size: properly-paired reads (flag 0x2), avoid double-counting
        if (flag & 0x2) and not (flag & 0x10) and tlen > 0:
            insert_sizes.append(tlen)

        # NM tag (edit distance)
        nm_match = re.search(r'\tNM:i:(\d+)', line)
        if nm_match:
            mismatches.append(int(nm_match.group(1)))

        count += 1

    proc.stdout.close()
    proc.wait()

    return {
        "insert_sizes": insert_sizes,
        "mapqs":         mapqs,
        "soft_clips":    soft_clips,
        "read_lengths":  read_lengths,
        "mismatches":    mismatches,
    }


# ── parse helpers ─────────────────────────────────────────────────────────────

def parse_flagstat(path: Path) -> dict:
    stats = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        count = int(parts[0])
        if "in total"   in line: stats["total"]     = count
        if "mapped ("   in line: stats["mapped"]    = count
        if "paired in"  in line: stats["paired"]    = count
        if "properly"   in line: stats["properly"]  = count
        if "singletons" in line: stats["singletons"]= count
        if "secondary"  in line: stats["secondary"] = count
        if "supplementary" in line: stats["supplementary"] = count
    return stats


def parse_depth(path: Path):
    positions, depths = [], []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split()
        positions.append(int(parts[1]))
        depths.append(int(parts[2]))
    return positions, depths


def parse_pileup(path: Path):
    """
    Return per-position base composition from mpileup output.
    Returns: positions, ref_bases, a_counts, t_counts, g_counts, c_counts
    """
    positions, ref_bases = [], []
    a_cnt, t_cnt, g_cnt, c_cnt = [], [], [], []

    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        pos        = int(parts[1])
        ref_base   = parts[2].upper()
        bases_raw  = parts[4].upper()

        # strip read-start/end markers and indel annotations
        bases_clean = re.sub(r'\^.', '', bases_raw)   # ^X (read start + mapq)
        bases_clean = re.sub(r'\$',  '', bases_clean) # $ (read end)
        bases_clean = re.sub(r'[+-]\d+[ACGTNacgtn]+', '', bases_clean)  # indels

        a = bases_clean.count('A')
        t = bases_clean.count('T')
        g = bases_clean.count('G')
        c = bases_clean.count('C')

        positions.append(pos)
        ref_bases.append(ref_base)
        a_cnt.append(a); t_cnt.append(t)
        g_cnt.append(g); c_cnt.append(c)

    return positions, ref_bases, a_cnt, t_cnt, g_cnt, c_cnt


# ── visualisation ─────────────────────────────────────────────────────────────

def visualise(label: str, flagstat: dict, positions: list, depths: list,
              pile_data: tuple, metrics: dict):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    pile_pos, pile_ref, a_cnt, t_cnt, g_cnt, c_cnt = pile_data

    PANEL  = "#f6f8fa"
    TEXT   = "#1f2328"
    ACCENT = "#0969da"
    GREEN  = "#1a7f37"
    YELLOW = "#9a6700"
    RED    = "#cf222e"
    PURPLE = "#8250df"
    GREY   = "#57606a"

    def style(ax, title=""):
        ax.set_facecolor(PANEL)
        ax.tick_params(colors=TEXT)
        for sp in ax.spines.values():
            sp.set_edgecolor("#d0d7de")
        ax.title.set_color(TEXT)
        ax.xaxis.label.set_color(TEXT)
        ax.yaxis.label.set_color(TEXT)
        if title:
            ax.set_title(title, fontsize=10, fontweight="bold", pad=7)

    fig = plt.figure(figsize=(26, 18), facecolor="white")
    gs  = GridSpec(4, 3, figure=fig, hspace=0.60, wspace=0.38)

    # ── 0. Pipeline overview ──────────────────────────────────────────────────
    ax0 = fig.add_subplot(gs[0, :])
    ax0.set_facecolor("white")
    ax0.axis("off")
    style(ax0, f"BWA-MEM2 Alignment Pipeline  ·  aligner used: {label}")

    steps = [
        ("1", "Index\nReference",     f"{label} index",         GREEN),
        ("2", "Align\nPaired-End",    f"{label} mem -R @RG",   ACCENT),
        ("3", "SAM → sorted\nBAM",    "samtools view | sort",  ACCENT),
        ("4", "Index BAM",            "samtools index",         ACCENT),
        ("5", "Alignment\nStats",     "samtools flagstat",      YELLOW),
        ("6", "Depth of\nCoverage",   "samtools depth -a",      YELLOW),
        ("7", "Base\nComposition",    "samtools mpileup",        RED),
    ]
    xs = 1.0 / len(steps)
    for i, (num, lbl, cmd, col) in enumerate(steps):
        cx = (i + 0.5) * xs
        if i < len(steps) - 1:
            ax0.annotate("", xy=(cx + xs * 0.40, 0.50),
                         xytext=(cx + xs * 0.10, 0.50),
                         xycoords="axes fraction",
                         arrowprops=dict(arrowstyle="->", color=GREY, lw=1.4))
        ax0.text(cx, 0.88, f"Step {num}", ha="center", va="center",
                 color=col, fontsize=7.5, fontweight="bold",
                 transform=ax0.transAxes)
        ax0.text(cx, 0.60, lbl, ha="center", va="center",
                 color=TEXT, fontsize=8.5, fontweight="bold",
                 transform=ax0.transAxes, linespacing=1.3)
        ax0.text(cx, 0.22, cmd, ha="center", va="center",
                 color=GREY, fontsize=6.5, fontstyle="italic",
                 transform=ax0.transAxes)

    # ── 1. Alignment stats bar ────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[1, 0])
    style(ax1, "Alignment Statistics (flagstat)")
    if flagstat:
        lbls   = ["Total", "Mapped", "Properly\nPaired", "Singletons",
                  "Secondary", "Supplementary"]
        keys   = ["total", "mapped", "properly", "singletons",
                  "secondary", "supplementary"]
        colors = [ACCENT, GREEN, GREEN, YELLOW, GREY, PURPLE]
        vals   = [flagstat.get(k, 0) for k in keys]
        bars   = ax1.bar(lbls, vals, color=colors,
                         edgecolor="#d0d7de", linewidth=0.8)
        for bar, val in zip(bars, vals):
            ax1.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + max(vals) * 0.01,
                     f"{val:,}", ha="center", va="bottom",
                     color=TEXT, fontsize=7.5)
        ax1.set_ylabel("Read Count")
        ax1.set_ylim(0, max(vals) * 1.18 if vals else 1)
        ax1.tick_params(axis="x", labelsize=8)

    # ── 2. Mapping rate pie ───────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 1])
    style(ax2, "Mapping Rate")
    total  = flagstat.get("total", 0)
    mapped = flagstat.get("mapped", 0)
    if total > 0:
        wedges, texts, autotexts = ax2.pie(
            [mapped, total - mapped],
            labels=["Mapped", "Unmapped"],
            colors=[GREEN, RED],
            autopct="%1.1f%%",
            startangle=90,
            wedgeprops=dict(edgecolor="#d0d7de", linewidth=1),
            textprops=dict(color=TEXT),
        )
        for at in autotexts:
            at.set_color("white"); at.set_fontweight("bold")

    # ── 3. Mapping quality distribution ──────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 2])
    style(ax3, "Mapping Quality Distribution")
    mqs = metrics.get("mapqs", [])
    if mqs:
        mqs_arr = [q for q in mqs if q > 0]   # exclude unmapped (MAPQ=0)
        import numpy as np
        ax3.hist(mqs_arr, bins=range(0, 62, 2), color=ACCENT,
                 edgecolor="#d0d7de", linewidth=0.6)
        ax3.axvline(np.mean(mqs_arr), color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Mean MQ = {np.mean(mqs_arr):.1f}")
        q20 = sum(1 for q in mqs_arr if q >= 20) / len(mqs_arr) * 100
        ax3.text(0.97, 0.95, f"MQ≥20: {q20:.1f}%",
                 ha="right", va="top", transform=ax3.transAxes,
                 color=GREEN, fontsize=8.5, fontweight="bold")
        ax3.set_xlabel("MAPQ Score")
        ax3.set_ylabel("Number of Reads")
        ax3.legend(facecolor="white", edgecolor="#d0d7de",
                   labelcolor=TEXT, fontsize=8)

    # ── 4. Depth of coverage ──────────────────────────────────────────────────
    import numpy as np
    ax4 = fig.add_subplot(gs[2, :2])
    style(ax4, "Sequencing Depth of Coverage Along Reference")
    if depths:
        pos_a   = np.array(positions)
        dep_a   = np.array(depths)
        ax4.fill_between(pos_a, dep_a, alpha=0.30, color=ACCENT)
        ax4.plot(pos_a, dep_a, color=ACCENT, lw=0.8)
        mean_d = dep_a.mean()
        ax4.axhline(mean_d, color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Mean depth = {mean_d:.1f}×")
        ax4.set_xlabel("Genomic Position")
        ax4.set_ylabel("Depth (×)")
        ax4.legend(facecolor="white", edgecolor="#d0d7de",
                   labelcolor=TEXT, fontsize=8)

    # ── 5. Depth histogram ────────────────────────────────────────────────────
    ax5 = fig.add_subplot(gs[2, 2])
    style(ax5, "Depth Distribution")
    if depths:
        dep_a = np.array(depths)
        ax5.hist(dep_a, bins=40, color=PURPLE,
                 edgecolor="#d0d7de", linewidth=0.6)
        ax5.axvline(dep_a.mean(), color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Mean = {dep_a.mean():.1f}×")
        ax5.set_xlabel("Depth (×)")
        ax5.set_ylabel("Number of Positions")
        ax5.legend(facecolor="white", edgecolor="#d0d7de",
                   labelcolor=TEXT, fontsize=8)

    # ── 6. Stacked base composition (pileup) ──────────────────────────────────
    ax6 = fig.add_subplot(gs[3, :2])
    style(ax6, "Per-Position Base Composition (mpileup)")
    if pile_pos:
        pp   = np.array(pile_pos,  dtype=float)
        a_a  = np.array(a_cnt,     dtype=float)
        t_a  = np.array(t_cnt,     dtype=float)
        g_a  = np.array(g_cnt,     dtype=float)
        c_a  = np.array(c_cnt,     dtype=float)
        tot  = a_a + t_a + g_a + c_a
        tot  = np.where(tot == 0, 1, tot)       # avoid /0
        a_f, t_f = a_a / tot, t_a / tot
        g_f, c_f = g_a / tot, c_a / tot

        # smooth with a running window to keep plot readable
        win = max(1, len(pp) // 500)
        def smooth(x):
            return np.convolve(x, np.ones(win) / win, mode='same')

        ax6.stackplot(pp, smooth(a_f), smooth(t_f),
                          smooth(g_f), smooth(c_f),
                      labels=["A", "T", "G", "C"],
                      colors=["#56d364", "#ff7b72", "#ffa657", "#79c0ff"],
                      alpha=0.85)
        ax6.set_xlabel("Genomic Position")
        ax6.set_ylabel("Fraction of Bases")
        ax6.set_ylim(0, 1)
        ax6.legend(loc="upper right", facecolor="white",
                   edgecolor="#d0d7de", labelcolor=TEXT, fontsize=8,
                   ncol=4)

    # ── 7. Insert size & soft-clip distributions ──────────────────────────────
    ax7 = fig.add_subplot(gs[3, 2])
    style(ax7, "Insert Size Distribution")
    isizes = metrics.get("insert_sizes", [])
    if isizes:
        is_arr = np.array(isizes)
        p99 = np.percentile(is_arr, 99)
        ax7.hist(is_arr[is_arr <= p99], bins=60,
                 color=GREEN, edgecolor="#d0d7de", linewidth=0.6)
        ax7.axvline(np.median(is_arr), color=YELLOW, linestyle="--", lw=1.5,
                    label=f"Median = {int(np.median(is_arr))} bp")
        ax7.set_xlabel("Insert Size (bp)")
        ax7.set_ylabel("Number of Read Pairs")
        ax7.legend(facecolor="white", edgecolor="#d0d7de",
                   labelcolor=TEXT, fontsize=8)
    else:
        ax7.text(0.5, 0.5, "No properly-paired\nreads found",
                 ha="center", va="center", color=TEXT,
                 transform=ax7.transAxes)

    # ── title & footer ────────────────────────────────────────────────────────
    fig.suptitle(
        f"BWA-MEM2 Alignment Report  ·  aligner: {label}\n"
        "Reference: example_human_reference.fasta  ·  Sample: example_human_Illumina",
        color=TEXT, fontsize=13, fontweight="bold", y=0.995,
    )
    fig.text(
        0.5, 0.003,
        f"Pipeline: {label} index → {label} mem → samtools sort/index "
        "→ flagstat / depth -a / mpileup",
        ha="center", color=GREY, fontsize=8,
    )

    plt.savefig(PLOT, dpi=160, bbox_inches="tight", facecolor="white")
    print(f"\n✓ Report saved → {PLOT}")
    return fig


# ── text summary ──────────────────────────────────────────────────────────────

def print_summary(label: str, flagstat: dict, depths: list, metrics: dict):
    import numpy as np
    sep = "─" * 60
    print(f"\n{sep}")
    print("  ALIGNMENT SUMMARY")
    print(sep)

    print(f"\n  Aligner : {label}")

    print("\n[1] Input")
    ref_bp = sum(len(l) for l in REF.read_text().splitlines()
                 if not l.startswith(">"))
    n_reads = sum(1 for _ in open(R1)) // 4
    print(f"    Reference : {REF.name}  ({ref_bp:,} bp)")
    print(f"    Reads (R1): {n_reads:,} pairs  (paired-end Illumina)")

    print("\n[2] Alignment (flagstat)")
    total  = flagstat.get("total",  0)
    mapped = flagstat.get("mapped", 0)
    prop   = flagstat.get("properly", 0)
    sing   = flagstat.get("singletons", 0)
    rate   = mapped / total * 100 if total else 0
    print(f"    Total reads      : {total:,}")
    print(f"    Mapped reads     : {mapped:,}  ({rate:.2f}%)")
    print(f"    Properly paired  : {prop:,}  ({prop/total*100:.2f}%)" if total else "")
    print(f"    Singletons       : {sing:,}")

    print("\n[3] Coverage (samtools depth -a)")
    if depths:
        d = np.array(depths)
        print(f"    Positions covered  : {len(d):,}")
        print(f"    Mean depth         : {d.mean():.2f}×")
        print(f"    Median depth       : {np.median(d):.0f}×")
        print(f"    Zero-depth pos.    : {(d == 0).sum():,}")
        print(f"    ≥10× positions     : {(d >= 10).sum() / len(d) * 100:.1f}%")
        print(f"    ≥20× positions     : {(d >= 20).sum() / len(d) * 100:.1f}%")

    print("\n[4] Read metrics (sample)")
    mqs = metrics.get("mapqs", [])
    iss = metrics.get("insert_sizes", [])
    scs = metrics.get("soft_clips", [])
    nms = metrics.get("mismatches", [])
    if mqs:
        mqs_a = np.array(mqs)
        print(f"    Mean MAPQ          : {mqs_a.mean():.1f}")
        print(f"    Reads with MQ≥20   : {(mqs_a >= 20).sum() / len(mqs_a) * 100:.1f}%")
    if iss:
        iss_a = np.array(iss)
        print(f"    Median insert size : {np.median(iss_a):.0f} bp")
        print(f"    Mean insert size   : {iss_a.mean():.0f} bp")
    if scs and metrics.get("read_lengths"):
        rl_a  = np.array(metrics["read_lengths"])
        sc_a  = np.array(scs)
        sc_pct = sc_a.sum() / rl_a.sum() * 100
        print(f"    Soft-clip fraction : {sc_pct:.2f}% of all bases")
    if nms:
        nm_a = np.array(nms)
        print(f"    Mean mismatches/read (NM): {nm_a.mean():.2f}")

    print(f"\n{sep}")
    print(f"  Output : {OUT}")
    print(f"  Plot   : {PLOT}")
    print(sep + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  BWA-MEM2 Alignment Pipeline")
    print("=" * 60)

    OUT.mkdir(parents=True, exist_ok=True)

    print("\n[Step 0] Checking tools…")
    aligner, label = pick_aligner()
    ensure_samtools()
    print(f"✓ Aligner : {label}  ({aligner})")

    print(f"\n[Step 1] Indexing reference with {label}…")
    step_index(aligner, label)

    print(f"\n[Step 2-3] Aligning reads with {label} mem → sorted BAM…")
    step_align(aligner, label)

    print("\n[Step 4] Alignment statistics (flagstat)…")
    step_flagstat()
    print(STATS.read_text())

    print("\n[Step 5] Per-base depth of coverage…")
    step_depth()

    print("\n[Step 6] Per-position base composition (mpileup)…")
    step_pileup()

    print("\n[Step 7] Collecting per-read metrics from BAM…")
    metrics = parse_read_metrics()

    # ── parse outputs ─────────────────────────────────────────────────────────
    flagstat          = parse_flagstat(STATS)
    positions, depths = parse_depth(DEPTH)
    pile_data         = parse_pileup(PILE)

    print_summary(label, flagstat, depths, metrics)

    print("\n[Step 8] Generating visualisation…")
    visualise(label, flagstat, positions, depths, pile_data, metrics)

    print("\n✓ Pipeline complete.")


if __name__ == "__main__":
    main()
