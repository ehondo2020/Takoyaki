import json
import sys
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd

GARD_JSON_PATH = "*.json" # Output_file by Hyphy GARD
FASTA_FILE_PATH = "*.fa"
PBR_TXT_PATH = "*.txt" # Output_file by PBR_nonPBR_generator_F.py
OUTPUT_IMAGE_PATH = "*.png"

TARGET_SPECIES_KEYWORD = "Homo_sapiens" # Specify species name on the fasta file headers

DOMAIN_CONFIGS = [
    {"name": r"$\alpha1$ Domain", "start": 1, "end": 270, "color": "#4C72B0"},
    {"name": r"$\alpha2$ Domain", "start": 271, "end": 546, "color": "#55A868"},
    {"name": r"$\alpha3$ Domain", "start": 547, "end": 822, "color": "#C44E52"},
    {"name": "TM / Cytoplasmic", "start": 823, "end": None, "color": "#8172B0"},
]

FIG_SIZE = (32, 14.5)
DPI = 600
TITLE_TEXT = "Human Mature MHC Class I with GARD Breakpoints & Bjorkman PBR Sites"
TITLE_FONT_SIZE = 34
X_LABEL_TEXT = "MHC Class I Mature Sequence Position (DNA bp)"
X_LABEL_FONT_SIZE = 24
X_TICK_FONT_SIZE = 18

TRACK_Y = 0.55
TRACK_HEIGHT = 0.15
TRACK_BORDER_COLOR = "black"
TRACK_ALPHA = 0.85
DOMAIN_LABEL_FONT_SIZE = 26

BP_LINE_COLOR = "#D62728"
BP_LINE_STYLE = "--"
BP_LINE_WIDTH = 1.2
BP_BOX_BG_COLOR = "#FFE6E6"
BP_LABEL_FONT_SIZE = 34
BP_TEXT_BASE_Y = TRACK_Y + 0.32
BP_TEXT_Y_OFFSET = 0.12

PBR_TRACK_Y = TRACK_Y - 0.12
PBR_TRACK_HEIGHT = 0.06
PBR_COLOR = "#FF7F0E"
PBR_ALPHA = 0.9
PBR_BORDER_COLOR = "#D95F02"

SHOW_PBR_LABELS = True
PBR_LABEL_FONT_SIZE = 16
PBR_LABEL_BASE_Y = PBR_TRACK_Y - 0.06
PBR_LABEL_Y_OFFSET = 0.2 

def load_gard_breakpoints(json_path):
    with open(json_path, "r", encoding="utf-8") as f:
        gard_data = json.load(f)
    improvements = gard_data.get("improvements", {})
    latest_step = max(improvements.keys(), key=int)
    return [bp[0] for bp in improvements[latest_step]["breakpoints"]]


def load_target_sequence(fasta_path, species_keyword):
    aligned_seq = ""
    with open(fasta_path, "r", encoding="utf-8") as f:
        header = ""
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if species_keyword in header:
                    aligned_seq = "".join(seq_lines)
                    break
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if species_keyword in header and not aligned_seq:
            aligned_seq = "".join(seq_lines)

    if not aligned_seq:
        raise ValueError(
            f"Failed to find sequences matching '{species_keyword}' in the FASTA file."
        )

    ungapped_seq = aligned_seq.replace("-", "")
    return aligned_seq, ungapped_seq


def convert_aligned_to_ungapped_pos(aligned_seq, aligned_pos_1based):
    sub_seq = aligned_seq[:aligned_pos_1based]
    return len(sub_seq.replace("-", ""))


def load_and_verify_pbr_table(pbr_txt_path, ungapped_dna):
    try:
        pbr_df = pd.read_csv(
            pbr_txt_path, sep=r"\s+|\t", engine="python", header=None
        )
    except Exception as e:
        raise FileNotFoundError(
            f"Failed to load the PBR table file '{pbr_txt_path}'"
        )

    pbr_dna_positions = []
    for idx, row in pbr_df.iterrows():
        try:
            aa_pos = int(row.iloc[0])
        except ValueError:
            continue

        expected_codon = (
            str(row.iloc[3]).strip().upper() if len(row) >= 4 else None
        )

        dna_start = (aa_pos - 1) * 3 + 1
        dna_end = aa_pos * 3

        if dna_end > len(ungapped_dna):
            continue

        actual_codon = ungapped_dna[dna_start - 1 : dna_end].upper()

        if expected_codon and expected_codon not in ["-", "NAN"]:
            if actual_codon != expected_codon:
                raise ValueError(
                    f"PBR mapping mismatch error: Codon mismatch at AA position {aa_pos}."
                )

        pbr_dna_positions.append(
            {
                "aa_pos": aa_pos,
                "dna_start": dna_start,
                "dna_end": dna_end,
                "codon": actual_codon,
            }
        )

    return pbr_dna_positions


def main():
    raw_breakpoints = load_gard_breakpoints(GARD_JSON_PATH)
    aligned_seq, ungapped_seq = load_target_sequence(
        FASTA_FILE_PATH, TARGET_SPECIES_KEYWORD
    )

    corrected_breakpoints = [
        convert_aligned_to_ungapped_pos(aligned_seq, bp)
        for bp in raw_breakpoints
    ]

    total_dna_len = len(ungapped_seq)
    pbr_list = load_and_verify_pbr_table(PBR_TXT_PATH, ungapped_seq)

    print("=== [Validation Log] Coordinate transformation results ===")
    print(f"Full-length mature DNA: {total_dna_len} bp")
    for orig, corr in zip(raw_breakpoints, corrected_breakpoints):
        print(
            f"Breakpoint -> Alignment pos: {orig} Col | Rendered mature DNA pos: {corr} bp"
        )
    print("===============================")

    fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)

    for dom in DOMAIN_CONFIGS:
        start = dom["start"]
        end = dom["end"] if dom["end"] is not None else total_dna_len
        box_start_x = start - 0.5
        width = (end - start) + 1.0

        rect = patches.Rectangle(
            (box_start_x, TRACK_Y),
            width,
            TRACK_HEIGHT,
            linewidth=1,
            edgecolor=TRACK_BORDER_COLOR,
            facecolor=dom["color"],
            alpha=TRACK_ALPHA,
        )
        ax.add_patch(rect)
        ax.text(
            box_start_x + width / 2,
            TRACK_Y + TRACK_HEIGHT / 2,
            dom["name"],
            color="white",
            weight="bold",
            fontsize=DOMAIN_LABEL_FONT_SIZE,
            ha="center",
            va="center",
        )

    ax.add_patch(
        patches.Rectangle(
            (0.5, PBR_TRACK_Y),
            total_dna_len,
            PBR_TRACK_HEIGHT,
            linewidth=0.5,
            edgecolor="#CCCCCC",
            facecolor="#F0F0F0",
            zorder=1,
        )
    )

    for idx, item in enumerate(pbr_list):
        dna_start = item["dna_start"]
        dna_end = item["dna_end"]

        box_start_x = dna_start - 0.5
        width = 3.0

        pbr_rect = patches.Rectangle(
            (box_start_x, PBR_TRACK_Y),
            width,
            PBR_TRACK_HEIGHT,
            linewidth=0.5,
            edgecolor=PBR_BORDER_COLOR,
            facecolor=PBR_COLOR,
            alpha=PBR_ALPHA,
            zorder=2,
        )
        ax.add_patch(pbr_rect)

        if SHOW_PBR_LABELS:
            label_y = PBR_LABEL_BASE_Y - (idx % 2) * PBR_LABEL_Y_OFFSET
            label_text = f"{item['aa_pos']} ({item['dna_start']}-{item['dna_end']})"

            center_x = float(dna_start + 1)

            ax.text(
                center_x,
                label_y,
                label_text,
                fontsize=PBR_LABEL_FONT_SIZE,
                ha="center",
                va="top",
                color="#222222",
                rotation=90,
            )
            ax.plot(
                [center_x, center_x],
                [PBR_TRACK_Y, label_y + 0.01],
                color="#AAAAAA",
                lw=0.4,
                zorder=1,
            )

    ax.text(
        -10,
        PBR_TRACK_Y + PBR_TRACK_HEIGHT / 2,
        "Bjorkman PBR Sites",
        fontsize=24,
        weight="bold",
        ha="right",
        va="center",
        color="#333333",
    )

    for i, bp_ungapped in enumerate(corrected_breakpoints):
        orig_bp = raw_breakpoints[i]
        text_y = BP_TEXT_BASE_Y + (i % 2) * BP_TEXT_Y_OFFSET

        bp_x = float(bp_ungapped)

        ax.annotate(
            f"Breakpoint {i+1}\nPositon: {bp_ungapped} bp", # f"Breakpoint {i+1}\nUngapped: {bp_ungapped} bp\n(Aligned: {orig_bp})",
            xy=(bp_x, TRACK_Y + TRACK_HEIGHT),
            xytext=(bp_x, text_y),
            arrowprops=dict(
                facecolor=BP_LINE_COLOR,
                edgecolor=BP_LINE_COLOR,
                shrink=0.05,
                width=1.8,
                headwidth=7,
            ),
            fontsize=BP_LABEL_FONT_SIZE,
            weight="bold",
            color="#8C0000",
            ha="center",
            va="bottom",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor=BP_BOX_BG_COLOR,
                edgecolor=BP_LINE_COLOR,
                lw=1,
            ),
            zorder=4,
        )

        ax.plot(
            [bp_x, bp_x],
            [PBR_TRACK_Y, text_y],
            color=BP_LINE_COLOR,
            linestyle=BP_LINE_STYLE,
            linewidth=BP_LINE_WIDTH,
            zorder=3,
        )

    ax.set_xlim(-15, total_dna_len + 30)

    tick_step = 100
    tick_values = [1] + list(range(100, total_dna_len + 1, tick_step))
    ax.set_xticks(tick_values)
    ax.set_xticklabels([str(v) for v in tick_values])

    ax.tick_params(axis="x", labelsize=X_TICK_FONT_SIZE)

    ax.set_ylim(0, 1.25)
    ax.set_xlabel(X_LABEL_TEXT, fontsize=X_LABEL_FONT_SIZE, labelpad=10)
    ax.set_yticks([])
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    plt.title(TITLE_TEXT, fontsize=TITLE_FONT_SIZE, weight="bold", pad=20)
    plt.tight_layout()

    plt.savefig(OUTPUT_IMAGE_PATH, dpi=DPI)
    print(
        f"[SUCCESS] Rendering complete: Saved to '{OUTPUT_IMAGE_PATH}' using the geometric centroid axis model."
    )
    plt.show()


if __name__ == "__main__":
    main()