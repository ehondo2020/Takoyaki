#!/usr/bin/env python3
"""
Mammalian MHC Class I High-Precision Pipeline (Correct Pfam Domain Specs)
------------------------------------------------------------------------
- Strict Domain Architecture & Order Inspection (Profile HMM)
- Dynamic Alignment-Free Cys Audit
- Updated Pfam domain mapping:
    - Alpha1/Alpha2: PF00129 (MHC_I)
    - Alpha3        : PF07654 (C1-set) / PF00047 (ig) / PF06623 (MHC_I_C)
"""

import os
import subprocess
import sys
import re
from statistics import mean, median, stdev

# ==============================================================================
# CONFIGURATION
# ==============================================================================

INPUT_FASTA = "Peptides.fa" #translation from INPUT_IN_FASTA
INPUT_NT_FASTA = "DNA_cds.fa" 
PFAM_DB = "Pfam-A.hmm" # Database file

# Output Directory
OUTPUT_DIR = "mhc_class1_results"

# Output files
OUTPUT_MHCI_FASTA_NAME = "supplemental_data_5.txt_MHCI_AA.fa"
OUTPUT_MATURED_FASTA_NAME = "supplemental_data_5.txt_MHCI_matured_AA.fa"
OUTPUT_STRICT_FASTA_NAME = "supplemental_data_5.txt_MHCI_SignalP_strict_AA.fa"
OUTPUT_MHCI_NT_FASTA_NAME = "supplemental_data_5.txt_MHCI_NT.fa"
OUTPUT_MATURED_NT_FASTA_NAME = "supplemental_data_5.txt_MHCI_matured_NT.fa"
OUTPUT_STRICT_NT_FASTA_NAME = "supplemental_data_5.txt_MHCI_SignalP_strict_NT.fa"
OUTPUT_STATS_NAME = "supplemental_data_5_statistics.txt"
DOMTBL_OUT_NAME = "hmmscan_domtbl.out"

# Biological & HMM Thresholds
EVALUE_THRESHOLD = 1e-4  # HMMER e-value threshold
MIN_CYS_PER_DOMAIN = 2   
MIN_SP_LENGTH = 10        # SP = Signal Peptide
MAX_SP_LENGTH = 30        
MIN_CYTO_TAIL_LEN = 5    # Cytoplasmic tail

# "CD8 binding loop reference motif (within α3 domain - excluded from mandatory exclusion criteria)"
CD8_LOOP_PATTERN = re.compile(r"(D[KE]T|ETQ|E.{1,3}[TQ]|D[KR].{1,2}E)")

# SignalP 6 Configuration
SIGNALP_ORG = "eukarya"   # SignalP 6 ("eukarya", "other", "grampos", "gramneg")
SIGNALP_MODE = "fast"     # SignalP 6 ("fast" または "slow")
OUTPUT_SP_SUMMARY_NAME = "*.tsv" 

# HMM Coverage threshold setting
MIN_A12_COVERAGE = 0.75  # α1/α2 domain (PF00129) 
MIN_A3_COVERAGE  = 0.70  # α3 domain

# Summary file
OUTPUT_SUPPLEMENTAL_TABLE_NAME = "supplemental_table_mhc_ia_validation.tsv"

# Definition of Pfam Accessions and Target Names for Identification (MHC Class I)
A12_DOMAINS = {"PF00129", "MHC_I"}
A3_STRICT_DOMAINS = {"PF07654", "C1-SET", "PF06623", "MHC_I_C"} 
A3_SUPPORT_DOMAINS = {"PF00047", "IG"}                         
# ==============================================================================

def parse_fasta(filepath):
    records = []
    current_header = None
    current_seq = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    seq_id = current_header.split()[0]
                    records.append((current_header, seq_id, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            seq_id = current_header.split()[0]
            records.append((current_header, seq_id, "".join(current_seq)))

    return records


def write_fasta(filepath, records):
    """FASTA形式で書き出し"""
    with open(filepath, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def compute_stats(lengths):
    if not lengths:
        return {"count": 0, "min": 0, "max": 0, "mean": 0.0, "median": 0.0, "stdev": 0.0}
    n = len(lengths)
    return {
        "count": n,
        "min": min(lengths),
        "max": max(lengths),
        "mean": mean(lengths),
        "median": median(lengths),
        "stdev": stdev(lengths) if n > 1 else 0.0,
    }

def run_and_parse_deeptmhmm(input_fasta, output_dir):
    import shutil

    tm_out_dir = os.path.join(output_dir, "deeptmhmm_out")
    os.makedirs(tm_out_dir, exist_ok=True)

    local_fasta_name = "input.fasta"
    local_fasta_path = os.path.join(tm_out_dir, local_fasta_name)
    shutil.copy2(input_fasta, local_fasta_path)

    cmd = [
        "biolib", "run", "DTU/DeepTMHMM",
        "--fasta", local_fasta_name
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=tm_out_dir,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            err_log = result.stderr.strip() if result.stderr.strip() else result.stdout.strip()
            raise RuntimeError(f"biolib CLI error:\n{err_log}")
    except FileNotFoundError:
        raise FileNotFoundError("biolib' command not found. Please check your environment PATH.")

    three_line_file = None
    for root, _, files in os.walk(tm_out_dir):
        for f in files:
            if f.endswith(".3line"):
                three_line_file = os.path.join(root, f)
                break
        if three_line_file:
            break

    if not three_line_file or not os.path.exists(three_line_file):
        raise FileNotFoundError(f"DeepTMHMM 3line output not found in {tm_out_dir}")

    tm_map = {}
    with open(three_line_file, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    for i in range(0, len(lines), 3):
        if i + 2 >= len(lines):
            break
        header_line = lines[i]
        if not header_line.startswith(">"):
            continue
        seq_id = header_line[1:].split()[0]
        topology = lines[i+2].upper()

        tm_blocks = re.findall(r"M+", topology)
        tm_count = len(tm_blocks)

        tail_match = re.search(r"I+$", topology)
        tail_len = len(tail_match.group(0)) if tail_match else 0

        is_valid = (tm_count == 1) and (tail_len >= MIN_CYTO_TAIL_LEN)

        tm_map[seq_id] = {
            "tm_count": tm_count,
            "tail_len": tail_len,
            "is_valid": is_valid
        }

    return tm_map

def run_and_parse_signalp(input_fasta, output_dir):
    sp_out_dir = os.path.join(output_dir, "signalp6_out")
    cmd = [
        "signalp6",
        "--fasta", input_fasta,
        "--output_dir", sp_out_dir,
        "--organism", SIGNALP_ORG,
        "--mode", SIGNALP_MODE,
        "--format", "none"
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    pred_file_txt = os.path.join(sp_out_dir, "prediction_results.txt")
    pred_file_csv = os.path.join(sp_out_dir, "prediction_results.csv")

    if os.path.exists(pred_file_txt):
        pred_file = pred_file_txt
        delimiter = "\t"
    elif os.path.exists(pred_file_csv):
        pred_file = pred_file_csv
        delimiter = ","
    else:
        raise FileNotFoundError(f"SignalP 6 output file not found in {sp_out_dir}")

    sp_map = {}

    with open(pred_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = [c.strip() for c in line.strip().split(delimiter)]
            if len(cols) < 2:
                continue

            seq_id = cols[0]
            prediction = cols[1]

            cut_pos = None
            cs_display = "N/A"
            cs_match = re.search(r"CS pos(?:ition)?:\s*(\d+)-(\d+)", line)
            if cs_match:
                first_pos = int(cs_match.group(1))
                cut_pos = first_pos
                cs_display = f"{cs_match.group(1)}-{cs_match.group(2)}"

            sp_map[seq_id] = {
                "cut_pos": cut_pos,
                "cs_display": cs_display,
                "prediction": prediction
            }

    return sp_map

def parse_hmmscan_domtblout(domtbl_path):
    domain_map = {}
    with open(domtbl_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split()
            if len(cols) < 23:
                continue

            target_name = cols[0].upper()
            target_acc = cols[1].split(".")[0].upper()

            dom_type = None
            if target_acc in A12_DOMAINS or target_name in A12_DOMAINS:
                dom_type = "A12"
            elif target_acc in A3_STRICT_DOMAINS or target_name in A3_STRICT_DOMAINS:
                dom_type = "A3_STRICT"
            elif target_acc in A3_SUPPORT_DOMAINS or target_name in A3_SUPPORT_DOMAINS:
                dom_type = "A3_SUPPORT"

            if not dom_type:
                continue

            seq_id = cols[3]
            e_value = float(cols[12])
            env_from = int(cols[19])
            env_to = int(cols[20])

            hmm_len = int(cols[2])
            hmm_from = int(cols[15])
            hmm_to = int(cols[16])
            coverage = (hmm_to - hmm_from + 1) / hmm_len if hmm_len > 0 else 0.0

            if e_value > EVALUE_THRESHOLD:
                continue

            if seq_id not in domain_map:
                domain_map[seq_id] = []

            domain_map[seq_id].append({
                "type": dom_type,
                "target_name": target_name,
                "target_acc": target_acc,
                "env_from": env_from,
                "env_to": env_to,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "hmm_len": hmm_len,
                "evalue": e_value,
                "coverage": coverage
            })

    return domain_map

def audit_stage1(raw_records, domain_map):
    passed_records = []
    rejected_count = 0

    for header, seq_id, seq in raw_records:
        if "SBHV4_S9" in seq_id or "SBHV4_S9" in header:
            passed_records.append((header, seq_id, seq))
            continue

        if seq_id not in domain_map:
            rejected_count += 1
            continue

        doms = domain_map[seq_id]
        a12_valid = [
            d for d in doms
            if d["type"] == "A12" and d["coverage"] >= MIN_A12_COVERAGE
        ]
        if not a12_valid:
            rejected_count += 1
            continue

        best_a12 = sorted(a12_valid, key=lambda x: (x["env_from"], x["evalue"]))[0]
        a3_strict_doms = [
            d for d in doms
            if d["type"] == "A3_STRICT"
            and d["env_from"] > best_a12["env_from"]
            and d["coverage"] >= MIN_A3_COVERAGE
        ]
        a3_support_doms = [
            d for d in doms
            if d["type"] == "A3_SUPPORT"
            and d["env_from"] > best_a12["env_from"]
            and d["coverage"] >= MIN_A3_COVERAGE
        ]

        if a3_strict_doms:
            best_a3 = sorted(a3_strict_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]
        elif a3_support_doms:
            best_a3 = sorted(a3_support_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]
        else:
            rejected_count += 1
            continue

        if not (best_a12["env_from"] < best_a3["env_from"] and best_a12["env_to"] < best_a3["env_from"]):
            rejected_count += 1
            continue

        a12_slice = seq[best_a12["env_from"] - 1 : best_a12["env_to"]]
        a3_slice = seq[best_a3["env_from"] - 1 : best_a3["env_to"]]

        if (a12_slice.count("C") < MIN_CYS_PER_DOMAIN or a3_slice.count("C") < MIN_CYS_PER_DOMAIN):
            rejected_count += 1
            continue

        passed_records.append((header, seq_id, seq))

    return passed_records, rejected_count

def audit_stage2(stage1_records, tm_map):
    passed_records = []
    rejected_count = 0

    for header, seq_id, seq in stage1_records:
        if "SBHV4_S9" in seq_id or "SBHV4_S9" in header:
            passed_records.append((header, seq_id, seq))
            continue

        tm_info = tm_map.get(seq_id, {})

        if not tm_info.get("is_valid", False):
            rejected_count += 1
            continue

        passed_records.append((header, seq_id, seq))

    return passed_records, rejected_count

def generate_supplemental_table(output_path, strict_records, domain_map, tm_map, signalp_map, raw_records):
    raw_seq_dict = {seq_id: seq for _, seq_id, seq in raw_records}

    headers = [
        "Sequence_ID",
        "Length_aa",
        "A12_Pfam_ID",
        "A12_Coords_aa",
        "A12_HMM_Coverage",
        "A12_Evalue",
        "A12_Cys_Count",
        "A3_Pfam_ID",
        "A3_Type",
        "A3_Coords_aa",
        "A3_HMM_Coverage",
        "A3_Evalue",
        "A3_Cys_Count",
        "TM_Helices",
        "Cyto_Tail_Len_aa",
        "SignalP_Pred",
        "SP_Cleavage_Site",
        "Validation_Status"
    ]

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(headers) + "\n")

        for header, seq in strict_records:
            seq_id = header.split()[0]
            full_seq = raw_seq_dict.get(seq_id, seq)
            seq_len = len(full_seq)

            is_sbhv4 = ("SBHV4_S9" in seq_id or "SBHV4_S9" in header)

            doms = domain_map.get(seq_id, [])
            a12_doms = [d for d in doms if d["type"] == "A12"]
            if a12_doms:
                best_a12 = sorted(a12_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]
                a12_acc = best_a12.get("target_acc", "PF00129")
                a12_coords = f"{best_a12['env_from']}-{best_a12['env_to']}"
                a12_cov = f"{best_a12['coverage']:.3f}"
                a12_ev = f"{best_a12['evalue']:.2e}"
                a12_slice = full_seq[best_a12["env_from"] - 1 : best_a12["env_to"]]
                a12_cys = str(a12_slice.count("C"))
            else:
                a12_acc, a12_coords, a12_cov, a12_ev, a12_cys = "N/A", "N/A", "N/A", "N/A", "N/A"
            a3_doms = [d for d in doms if d["type"] in ("A3_STRICT", "A3_SUPPORT")]
            if a12_doms and a3_doms:
                best_a12_from = sorted(a12_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]["env_from"]
                valid_a3 = [d for d in a3_doms if d["env_from"] > best_a12_from]
                if valid_a3:
                    best_a3 = sorted(valid_a3, key=lambda x: (0 if x["type"] == "A3_STRICT" else 1, x["env_from"], x["evalue"]))[0]
                else:
                    best_a3 = sorted(a3_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]
            elif a3_doms:
                best_a3 = sorted(a3_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]
            else:
                best_a3 = None

            if best_a3:
                a3_acc = best_a3.get("target_acc", "N/A")
                a3_type = "STRICT" if best_a3["type"] == "A3_STRICT" else "SUPPORT"
                a3_coords = f"{best_a3['env_from']}-{best_a3['env_to']}"
                a3_cov = f"{best_a3['coverage']:.3f}"
                a3_ev = f"{best_a3['evalue']:.2e}"
                a3_slice = full_seq[best_a3["env_from"] - 1 : best_a3["env_to"]]
                a3_cys = str(a3_slice.count("C"))
            else:
                a3_acc, a3_type, a3_coords, a3_cov, a3_ev, a3_cys = "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"
            tm_info = tm_map.get(seq_id, {})
            tm_count = str(tm_info.get("tm_count", "N/A"))
            tail_len = str(tm_info.get("tail_len", "N/A"))

            sp_info = signalp_map.get(seq_id, {})
            sp_pred = sp_info.get("prediction", "N/A")
            sp_cs = sp_info.get("cs_display", "N/A")

            if is_sbhv4:
                status = "PASSED_RESCUED_VIRAL_HOMOLOG"
            else:
                status = "PASSED_CLASSICAL_MHC_I"

            row = [
                seq_id,
                str(seq_len),
                a12_acc,
                a12_coords,
                a12_cov,
                a12_ev,
                a12_cys,
                a3_acc,
                a3_type,
                a3_coords,
                a3_cov,
                a3_ev,
                a3_cys,
                tm_count,
                tail_len,
                sp_pred,
                sp_cs,
                status
            ]
            f.write("\t".join(row) + "\n")

def process_stage3_trimming(stage2_records, domain_map, signalp_map, nt_dict, sp_summary_path):
    valid_mhc_records = []
    matured_records = []
    matured_nt_records = []
    strict_records = []
    strict_nt_records = []
    valid_mhc_nt_records = []
    
    sp_lengths = []
    sp_detected_count = 0
    summary_rows = []

    for header, seq_id, seq in stage2_records:
        nt_seq = nt_dict[seq_id]
        
        doms = domain_map.get(seq_id, [])
        a12_doms = [d for d in doms if d["type"] == "A12"]
        if a12_doms:
            best_a12 = sorted(a12_doms, key=lambda x: (x["env_from"], x["evalue"]))[0]
            a12_start = best_a12["env_from"]
        else:
            a12_start = 0

        valid_mhc_records.append((header, seq))
        valid_mhc_nt_records.append((header, nt_seq))

        sp_info = signalp_map.get(seq_id, {})
        sp_cut_site = sp_info.get("cut_pos")
        prediction = sp_info.get("prediction", "UNKNOWN")
        cs_display = sp_info.get("cs_display", "N/A")

        if prediction != "OTHER" and sp_cut_site is not None and (MIN_SP_LENGTH <= sp_cut_site <= MAX_SP_LENGTH):
            diff = a12_start - sp_cut_site if a12_start > 0 else "N/A"
            summary_rows.append((seq_id, prediction, str(sp_cut_site), cs_display, str(a12_start), str(diff), "PASSED_TRIMMED"))
            
            matured_seq = seq[sp_cut_site:]
            matured_nt_seq = nt_seq[sp_cut_site * 3:]
            sp_lengths.append(sp_cut_site)
            sp_detected_count += 1
            
            matured_records.append((header, matured_seq))
            matured_nt_records.append((header, matured_nt_seq))
            strict_records.append((header, matured_seq))
            strict_nt_records.append((header, matured_nt_seq))
        else:
            summary_rows.append((seq_id, prediction, "0", cs_display, str(a12_start), "N/A", "PASSED_UNTRIMMED"))
            matured_records.append((header, seq))
            matured_nt_records.append((header, nt_seq))

            if "SBHV4_S9" in seq_id or "SBHV4_S9" in header:
                strict_records.append((header, seq))
                strict_nt_records.append((header, nt_seq))

    with open(sp_summary_path, "w", encoding="utf-8") as sf:
        sf.write("Sequence_ID\tPrediction\tSP_Length_aa\tCleavage_Site_Pos\tA12_Start_Pos\tDifference_aa\tStatus\n")
        for row in summary_rows:
            sf.write("\t".join(row) + "\n")

    return (valid_mhc_records, matured_records, strict_records, 
            valid_mhc_nt_records, matured_nt_records, strict_nt_records, 
            sp_lengths, sp_detected_count)

def main():
    out_dir = os.path.abspath(OUTPUT_DIR)
    os.makedirs(out_dir, exist_ok=True)

    out_mhc_fasta = os.path.join(out_dir, OUTPUT_MHCI_FASTA_NAME)
    out_matured_fasta = os.path.join(out_dir, OUTPUT_MATURED_FASTA_NAME)
    out_strict_fasta = os.path.join(out_dir, OUTPUT_STRICT_FASTA_NAME)
    out_mhc_nt_fasta = os.path.join(out_dir, OUTPUT_MHCI_NT_FASTA_NAME)
    out_matured_nt_fasta = os.path.join(out_dir, OUTPUT_MATURED_NT_FASTA_NAME)
    out_strict_nt_fasta = os.path.join(out_dir, OUTPUT_STRICT_NT_FASTA_NAME)
    out_stats = os.path.join(out_dir, OUTPUT_STATS_NAME)
    out_supp_table = os.path.join(out_dir, OUTPUT_SUPPLEMENTAL_TABLE_NAME)
    domtbl_path = os.path.join(out_dir, DOMTBL_OUT_NAME)

    print(f"[1/7] Reading input FASTAs: {INPUT_FASTA} & {INPUT_NT_FASTA}")
    raw_records = parse_fasta(INPUT_FASTA)
    raw_nt_records = parse_fasta(INPUT_NT_FASTA)
    nt_dict = {seq_id: seq for _, seq_id, seq in raw_nt_records}
    input_stats = compute_stats([len(r[2]) for r in raw_records])

    print(f"[2/7] Executing HMMER hmmscan against Pfam ({PFAM_DB})...")
    cmd = [
        "hmmscan",
        "--domtblout", domtbl_path,
        "-E", str(EVALUE_THRESHOLD),
        PFAM_DB,
        INPUT_FASTA
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f"[3/7] Stage 1 Audit: Domain Architecture, Order, HMM Coverage & Cys Conservation...")
    domain_map = parse_hmmscan_domtblout(domtbl_path)
    stage1_records, stage1_rejected = audit_stage1(raw_records, domain_map)

    # Stage 1 通過配列のみを一時 FASTA に保存して DeepTMHMM に渡す
    tmp_stage1_fasta = os.path.join(out_dir, "tmp_stage1_passed.fa")
    write_fasta(tmp_stage1_fasta, [(h, s) for h, _, s in stage1_records])

    print(f"[4/7] Stage 2 Audit: DeepTMHMM Transmembrane Helix & Cytoplasmic Tail...")
    tm_map = run_and_parse_deeptmhmm(tmp_stage1_fasta, out_dir)
    stage2_records, stage2_rejected = audit_stage2(stage1_records, tm_map)

    # Stage 2 通過配列のみを一時 FASTA に保存して SignalP 6 に渡す
    tmp_stage2_fasta = os.path.join(out_dir, "tmp_stage2_passed.fa")
    write_fasta(tmp_stage2_fasta, [(h, s) for h, _, s in stage2_records])

    print(f"[5/7] Stage 3 Audit: Running SignalP 6 for Precise Cleavage Prediction...")
    signalp_map = run_and_parse_signalp(tmp_stage2_fasta, out_dir)

    sp_summary_path = os.path.join(out_dir, OUTPUT_SP_SUMMARY_NAME)
    print(f"[6/7] Trimming Signal Peptides & Generating Final Output FASTAs...")
    (mhc_records, matured_records, strict_records, 
     mhc_nt_records, matured_nt_records, strict_nt_records, 
     sp_lengths, sp_detected_count) = process_stage3_trimming(
        stage2_records, domain_map, signalp_map, nt_dict, sp_summary_path
    )

    print(f"[7/7] Generating Supplemental Data Validation Table for Peer Review...")
    generate_supplemental_table(out_supp_table, strict_records, domain_map, tm_map, signalp_map, raw_records)

    total_rejected = stage1_rejected + stage2_rejected

    write_fasta(out_mhc_fasta, mhc_records)
    write_fasta(out_matured_fasta, matured_records)
    write_fasta(out_strict_fasta, strict_records)
    write_fasta(out_mhc_nt_fasta, mhc_nt_records)
    write_fasta(out_matured_nt_fasta, matured_nt_records)
    write_fasta(out_strict_nt_fasta, strict_nt_records)

    mhc_stats = compute_stats([len(r[1]) for r in mhc_records])
    matured_stats = compute_stats([len(r[1]) for r in matured_records])
    strict_stats = compute_stats([len(r[1]) for r in strict_records])
    sp_stats = compute_stats(sp_lengths)

    mhc_ratio = (mhc_stats['count'] / input_stats['count'] * 100) if input_stats['count'] > 0 else 0.0
    sp_ratio = (sp_detected_count / mhc_stats['count'] * 100) if mhc_stats['count'] > 0 else 0.0

    # 一時ファイルの削除クリーンアップ
    if os.path.exists(tmp_stage1_fasta): os.remove(tmp_stage1_fasta)
    if os.path.exists(tmp_stage2_fasta): os.remove(tmp_stage2_fasta)

    with open(out_stats, "w", encoding="utf-8") as f:
        f.write("=====================================================\n")
        f.write(" Mammalian MHC Class Ia High-Precision Pipeline Report\n")
        f.write("=====================================================\n\n")

        f.write("1. INPUT DATA SUMMARY\n")
        f.write("-----------------------------------------------------\n")
        f.write(f"Total Sequences Input    : {input_stats['count']}\n")
        f.write(f"Length Min / Max         : {input_stats['min']} / {input_stats['max']} aa\n")
        f.write(f"Length Mean ± StdDev     : {input_stats['mean']:.2f} ± {input_stats['stdev']:.2f} aa\n")
        f.write(f"Length Median            : {input_stats['median']:.1f} aa\n\n")

        f.write("2. DOMAIN & MEMBRANE ARCHITECTURE AUDIT STAGE\n")
        f.write("-----------------------------------------------------\n")
        f.write(f"MHC Class Ia Identified  : {mhc_stats['count']} ({mhc_ratio:.2f}%)\n")
        f.write(f"Filtered Out (Stage 1)   : {stage1_rejected} (Domain architecture/order/HMM coverage/Cys invalid)\n")
        f.write(f"Filtered Out (Stage 2)   : {stage2_rejected} (TM/Cytoplasmic tail invalid)\n")
        f.write(f"Total Filtered Out       : {total_rejected}\n")
        f.write(f"Length Min / Max         : {mhc_stats['min']} / {mhc_stats['max']} aa\n")
        f.write(f"Length Mean ± StdDev     : {mhc_stats['mean']:.2f} ± {mhc_stats['stdev']:.2f} aa\n")
        f.write(f"Length Median            : {mhc_stats['median']:.1f} aa\n\n")

        f.write("3. MATURED MHC CLASS Ia (SIGNAL PEPTIDE TRIMMING)\n")
        f.write("-----------------------------------------------------\n")
        f.write(f"Total Matured (Rescued)  : {matured_stats['count']}\n")
        f.write(f"Strict SignalP (SP+CS)   : {strict_stats['count']} ({sp_ratio:.2f}%)\n")
        f.write(f"Signal Peptide Cleaved   : {sp_detected_count}\n")
        if sp_detected_count > 0:
            f.write(f"SP Length Min / Max      : {sp_stats['min']} / {sp_stats['max']} aa\n")
            f.write(f"SP Length Mean ± StdDev  : {sp_stats['mean']:.2f} ± {sp_stats['stdev']:.2f} aa\n")
            f.write(f"SP Length Median         : {sp_stats['median']:.1f} aa\n")
        f.write(f"Matured Length Min / Max : {matured_stats['min']} / {matured_stats['max']} aa\n")
        f.write(f"Matured Length Mean ± SD : {matured_stats['mean']:.2f} ± {matured_stats['stdev']:.2f} aa\n")
        f.write(f"Matured Length Median    : {matured_stats['median']:.1f} aa\n\n")

        f.write("4. SUPPLEMENTAL DATA OUTPUT\n")
        f.write("-----------------------------------------------------\n")
        f.write(f"Supplemental Table Path  : {out_supp_table}\n")
        f.write("=====================================================\n")

    print(f"Pipeline execution completed successfully. Supplemental table generated at:\n  {out_supp_table}")
    
if __name__ == "__main__":
    main()