import os
import sys

AA_INPUT_FILE = "*.fa" #Peptide file by Macse alignment software
DNA_INPUT_FILE = "*.fa" #DNA file by Macse alignment software

AA_PBR_OUT = "*.fasta"
AA_NON_PBR_OUT = "*.fasta"
DNA_PBR_OUT = "*.fasta"
DNA_NON_PBR_OUT = "*.fasta"

# PBR coordinate mapping table
PBR_MAPPING_OUT = "*.txt"

# Header ID of reference sequences
REF_SEQ_ID = "Homo_sapiens_GCF_000001405"

# Mature protein offset (e.g., set to 24 if the reference sequence contains a 24 aa signal peptide)
MATURE_PROTEIN_OFFSET = 0

# ==============================================================================
# 2. PBR (Peptide-Binding Residues) definition setting
# ==============================================================================
# PBR definition set: "BJORKMAN" or "HUGHES_NEI"
PBR_MODE = "BJORKMAN"

PBR_DEFINITIONS = {
    # Bjorkman et al. (1987)
    "BJORKMAN": {
        # alpha1 domain (18 AA residues)
        5,
        7,
        9,
        22,
        24,
        45,
        59,
        62,
        63,
        66,
        67,
        70,
        73,
        74,
        77,
        80,
        81,
        84,
        # alpha2 domain (16 AA residues)
        95,
        97,
        99,
        114,
        116,
        123,
        143,
        146,
        147,
        150,
        152,
        156,
        159,
        163,
        167,
        171,
    },
    # Hughes & Nei (1988) 
    "HUGHES_NEI": {
        # alpha1 domain (27 AA residues)
        5,
        7,
        9,
        22,
        24,
        26,
        30,
        43,
        45,
        50,
        52,
        53,
        55,
        58,
        59,
        62,
        63,
        66,
        67,
        70,
        73,
        74,
        77,
        80,
        81,
        84,
        88,
        # alpha2 domain (30 AA residues)
        95,
        97,
        99,
        114,
        116,
        123,
        124,
        133,
        138,
        143,
        146,
        147,
        150,
        151,
        152,
        156,
        158,
        159,
        162,
        163,
        166,
        167,
        170,
        171,
        174,
        175,
        177,
        178,
        180,
        184,
    },
}

if PBR_MODE in PBR_DEFINITIONS:
    PBR_POSITIONS = PBR_DEFINITIONS[PBR_MODE]
else:
    raise ValueError(
        f"Invalid PBR_MODE: '{PBR_MODE}'. Must be 'BJORKMAN' or 'HUGHES_NEI'."
    )


def read_fasta(filepath):
    sequences = {}
    current_id = None
    current_seq = []

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(filepath, sequences):
    with open(filepath, "w", encoding="utf-8") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def main():
    print(f"PBR definition mode: {PBR_MODE} (Designated residues: {len(PBR_POSITIONS)})")
    print("Reading fasta files...")
    aa_dict = read_fasta(AA_INPUT_FILE)
    dna_dict = read_fasta(DNA_INPUT_FILE)

    first_aa_id = next(iter(aa_dict))
    alignment_aa_len = len(aa_dict[first_aa_id])

    pbr_col_indices = set()
    mapping_rows = []

    ref_key = None
    if REF_SEQ_ID:
        for seq_id in aa_dict:
            if REF_SEQ_ID in seq_id:
                ref_key = seq_id
                break
    
    for seq_id, aa_seq in aa_dict.items():
        if seq_id in dna_dict:
            dna_seq = dna_dict[seq_id]
            if len(dna_seq) != len(aa_seq) * 3:
                raise ValueError(
                    f"Fatal error: DNA length ({len(dna_seq)} bp) for {seq_id} "
                    f"is not exactly 3 times the AA length ({len(aa_seq)} aa). "
                    f"The sequence may not be codon-aligned."
                )

    if ref_key:
        ref_aa_clean = aa_dict[ref_key].replace("-", "").replace(".", "")
        print(f"[Check] N-terminal 5 residues of mature reference sequence ({ref_key}): {ref_aa_clean[:5]}")
        if not ref_aa_clean.startswith("G"):
            print(f"[WARNING] First residue of reference sequence is not 'G' (Glycine). "
                  f"Signal peptide may still be present, or the N-terminus might be truncated. "
                  f"Please check the value of MATURE_PROTEIN_OFFSET.")

    if ref_key:
        ref_seq = aa_dict[ref_key]
        ref_dna = dna_dict[ref_key]

        res_count = 0

        for col_idx, aa_char in enumerate(ref_seq):
            if aa_char.isalpha():
                res_count += 1
                hla_residue_num = res_count - MATURE_PROTEIN_OFFSET

                if hla_residue_num in PBR_POSITIONS:
                    pbr_col_indices.add(col_idx)
                    codon = ref_dna[col_idx * 3 : col_idx * 3 + 3]
                    mapping_rows.append((hla_residue_num, col_idx, aa_char, codon))

        print(
            f"Identified {len(pbr_col_indices)} PBR columns based on the reference sequence ({ref_key})."
        )
    else:
        if REF_SEQ_ID:
            raise KeyError(
                f"Error: Reference sequence ID '{REF_SEQ_ID}' was not found in the input FASTA file. "
                f"Please check the ID formatting in your FASTA headers."
            )
        pbr_col_indices = {pos - 1 for pos in PBR_POSITIONS if pos <= alignment_aa_len}
        print(
            f"Specified {len(pbr_col_indices)} columns as PBR using alignment column indices directly."
        )

    all_col_indices = list(range(alignment_aa_len))
    sorted_pbr_cols = sorted(list(pbr_col_indices))
    sorted_non_pbr_cols = [c for c in all_col_indices if c not in pbr_col_indices]

    print("--------------------------------------------------")
    print(f"Selected PBR Definition : {PBR_MODE}")
    print(f"Alignment AA Length      : {alignment_aa_len} aa")
    print(f"PBR columns              : {len(sorted_pbr_cols)} codons")
    print(f"nonPBR columns           : {len(sorted_non_pbr_cols)} codons")
    print("--------------------------------------------------")

    aa_pbr_dict = {}
    aa_non_pbr_dict = {}
    dna_pbr_dict = {}
    dna_non_pbr_dict = {}

    for seq_id, aa_seq in aa_dict.items():
        if seq_id not in dna_dict:
            print(
                f"Warning: {seq_id} was not found in the DNA file and will be skipped."
            )
            continue

        dna_seq = dna_dict[seq_id]

        if len(dna_seq) != len(aa_seq) * 3:
            raise ValueError(
                f"Fatal error: DNA length ({len(dna_seq)} bp) and AA length ({len(aa_seq)} aa) for {seq_id} do not have a 1:3 ratio."
            )

        aa_pbr = "".join([aa_seq[i] for i in sorted_pbr_cols])
        aa_non_pbr = "".join([aa_seq[i] for i in sorted_non_pbr_cols])

        dna_pbr = "".join([dna_seq[3 * i : 3 * i + 3] for i in sorted_pbr_cols])
        dna_non_pbr = "".join([dna_seq[3 * i : 3 * i + 3] for i in sorted_non_pbr_cols])

        dna_pbr = dna_pbr.replace("!", "-")
        dna_non_pbr = dna_non_pbr.replace("!", "-")

        aa_pbr = aa_pbr.replace("!", "-").replace("*", "X")
        aa_non_pbr = aa_non_pbr.replace("!", "-").replace("*", "X")

        aa_pbr_dict[seq_id] = aa_pbr
        aa_non_pbr_dict[seq_id] = aa_non_pbr
        dna_pbr_dict[seq_id] = dna_pbr
        dna_non_pbr_dict[seq_id] = dna_non_pbr

    write_fasta(AA_PBR_OUT, aa_pbr_dict)
    write_fasta(AA_NON_PBR_OUT, aa_non_pbr_dict)
    write_fasta(DNA_PBR_OUT, dna_pbr_dict)
    write_fasta(DNA_NON_PBR_OUT, dna_non_pbr_dict)

    with open(PBR_MAPPING_OUT, "w", encoding="utf-8") as f:

        f.write(
            "Mature_Position\tAlignment_Column(0-based)\tAA\tCodon\n"
        )

        for pos, col, aa, codon in sorted(mapping_rows):

            f.write(
                f"{pos}\t{col}\t{aa}\t{codon}\n"
            )

    print(f"Exported PBR coordinate mapping table: {PBR_MAPPING_OUT}")

    print("Analysis successfully completed.")

if __name__ == "__main__":
    main()