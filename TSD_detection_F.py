import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def load_fasta(file_path):
    print(f"Loading FASTA file: {file_path}")
    with open(file_path, "r") as file:
        return str(next(SeqIO.parse(file, "fasta")).seq)

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

REF_FASTA_PATH = "SBHV4_genome.fa" #viral genome sequence
QUERY_FASTA_PATH = "SBHV4_S9_full_length.fa" # viral target sequence
LOG_FILE_PATH = "TSD_analysis_execution.log" # output logfile

FLANKING_LENGTH = 100 

TSD_MIN_LEN = 8              # TSD minimum length (bp)
TSD_MAX_LEN = 20             # TSD longest length (bp)
TSD_STEP = 1                 # Subsequence step size (bp)
MAX_BOUNDARY_OFFSET = 5      # Maximum allowed distance from junction to TSD (bp)
MAX_MISMATCH_ALLOWED = 1     # 0 or 1

ALIGN_MATCH_SCORE = 1.0      # match score
ALIGN_MISMATCH_SCORE = -1.0  # mismatch penalty
ALIGN_OPEN_GAP_SCORE = -2.0  # gap open penalty
ALIGN_EXTEND_GAP_SCORE = -0.5 # gap extension penalty

NUM_PERMUTATIONS = 1000      # permutation iterations (Search-level test)
FDR_THRESHOLD = 0.05         # FDR threashold
STRICT_P_VAL_THRESHOLD = 0.01 # Secondary strict p-value threshold for filtering strict poly(A) output files

PAS_MOTIFS_PLUS = [
    "AATAAA", "ATTAAA", 
    "TATAAA", "AGTAAA", "AATATA", "CATAAA", "GATAAA", "AATAGA", "ACTAAA"
]                            # PAS motif
POLYA_MIN_LEN = 15           # Minimum length of poly(A) tract (bp) (stringent: 10 -> 15 bp)
POLYA_MAX_LEN = 30           # Maximum length of poly(A) tract (bp)
POLYA_MIN_A_RATIO = 0.90     # Minimum A-content ratio within window (stringent: 0.70 -> 0.90)
POLYA_MIN_CONSECUTIVE = 8    # Minimum required consecutive A bases (stringent: 4 -> 8 bp)

PAS_REQUIRE_PAS = True       # Require PAS to demonstrate authentic poly(A) (True)
PAS_SEARCH_OFFSET_MIN = 10   # Minimum distance from poly(A) to PAS (bp)
PAS_SEARCH_OFFSET_MAX = 40   # Maximum distance from poly(A) to PAS (bp)

def sliding_window_match(ref, query):
    print("Performing sliding window matching...")
    max_match, best_pos, best_variant = 0, (0, 0), None
    query_variants = [("original", query), ("reverse", query[::-1]), ("reverse_complement", reverse_complement(query))]

    for name, q in query_variants:
        for i in range(len(ref) - len(q) + 1):
            if i % 10000 == 0:
                print(f"Processing reference position: {i}")
            match = sum(1 for a, b in zip(ref[i:i+len(q)], q) if a == b)
            match_percentage = match / len(q)
            if match_percentage > max_match:
                max_match, best_pos, best_variant = match_percentage, (i, i + len(q)), name

    print(f"Best match found at: {best_pos} with match percentage: {max_match:.2f} ({best_variant})")
    return best_pos[0], best_pos[1], best_variant, max_match

def extract_flanking_sequences(ref, start, end, length=FLANKING_LENGTH):
    print(f"Extracting flanking sequences (length={length} bp)...")
    seq_f_start = max(0, start - length)
    seq_f_end = start
    seq_r_start = end
    seq_r_end = min(len(ref), end + length)
    seq_f = ref[seq_f_start:seq_f_end]
    seq_r = ref[seq_r_start:seq_r_end]
    return seq_f, seq_r, (seq_f_start, seq_f_end), (seq_r_start, seq_r_end)

def generate_subsequences(seq, start_pos, min_len=TSD_MIN_LEN, max_len=TSD_MAX_LEN, step=TSD_STEP, max_offset=MAX_BOUNDARY_OFFSET, side='end'):
    print(f"Generating boundary-constrained subsequences (min_len={min_len}, max_len={max_len}, max_offset={max_offset}bp, side={side})...")
    subsequences = []
    seq_len = len(seq)
    for length in range(min_len, max_len + 1):
        for i in range(0, seq_len - length + 1, step):
            offset_from_boundary = (seq_len - (i + length)) if side == 'end' else i
            if offset_from_boundary <= max_offset:
                abs_start = start_pos + i
                abs_end = abs_start + length
                subsequences.append((f"ID_{i}_{length}", seq[i:i+length], abs_start, abs_end))
    return subsequences

def calculate_hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Sequences must be of the same length for Hamming distance calculation.")
    return sum(1 for char1, char2 in zip(s1, s2) if char1 != char2)

def detect_strict_polya(ref_seq, start, end, strand='+'):
    pas_minus = [reverse_complement(m) for m in PAS_MOTIFS_PLUS]
    
    if strand == '+':
        win_start = max(0, end - 60)
        win_end = min(len(ref_seq), end + 40)
        region_seq = ref_seq[win_start:win_end].upper()
        
        for l in range(POLYA_MAX_LEN, POLYA_MIN_LEN - 1, -1):
            for i in range(len(region_seq) - l + 1):
                sub = region_seq[i:i+l]
                a_ratio = sub.count('A') / float(l)
                consecutive_a = 'A' * POLYA_MIN_CONSECUTIVE
                
                if a_ratio >= POLYA_MIN_A_RATIO and consecutive_a in sub:
                    pas_search_start = max(0, i - PAS_SEARCH_OFFSET_MAX)
                    pas_search_end = max(0, i - PAS_SEARCH_OFFSET_MIN)
                    upstream_region = region_seq[pas_search_start:pas_search_end]
                    
                    has_pas = any(pas in upstream_region for pas in PAS_MOTIFS_PLUS)
                    
                    if PAS_REQUIRE_PAS:
                        if has_pas:
                            return sub, True
                    else:
                        return sub, True
                        
        return "None", False

    elif strand == '-':
        win_start = max(0, start - 40)
        win_end = min(len(ref_seq), start + 60)
        region_seq = ref_seq[win_start:win_end].upper()
        
        for l in range(POLYA_MAX_LEN, POLYA_MIN_LEN - 1, -1):
            for i in range(len(region_seq) - l + 1):
                sub = region_seq[i:i+l]
                t_ratio = sub.count('T') / float(l)
                consecutive_t = 'T' * POLYA_MIN_CONSECUTIVE
                
                if t_ratio >= POLYA_MIN_A_RATIO and consecutive_t in sub:
                    pas_search_start = min(len(region_seq), i + l + PAS_SEARCH_OFFSET_MIN)
                    pas_search_end = min(len(region_seq), i + l + PAS_SEARCH_OFFSET_MAX)
                    downstream_region = region_seq[pas_search_start:pas_search_end]
                    
                    has_pas = any(pas in downstream_region for pas in pas_minus)
                    
                    if PAS_REQUIRE_PAS:
                        if has_pas:
                            polya_on_minus = reverse_complement(sub)
                            return polya_on_minus, True
                    else:
                        polya_on_minus = reverse_complement(sub)
                        return polya_on_minus, True
                        
        return "None", False
        
    return "None", False

def calculate_empirical_p_value(seq_f, seq_r, search_max_offset, observed_score, aligner, num_permutations=NUM_PERMUTATIONS, target_side='start'):
    count_equal_or_better = 0
    seq_r_list = list(seq_r)
    target_len = len(seq_r)
    query_len = len(seq_f)
    
    if target_side == 'start':
        min_pos = 0
        max_pos = min(search_max_offset + 1, target_len - query_len + 1)
    else:
        max_pos = target_len - query_len + 1
        min_pos = max(0, target_len - query_len - search_max_offset)
        
    for _ in range(num_permutations):
        random.shuffle(seq_r_list)
        shuffled_seq_r = "".join(seq_r_list)
        
        max_perm_score = float('-inf')
        for i in range(min_pos, max_pos):
            alignments = aligner.align(seq_f, shuffled_seq_r[i:i+query_len])
            if alignments:
                score = max(aln.score for aln in alignments)
                if score > max_perm_score:
                    max_perm_score = score
        
        if max_perm_score >= observed_score:
            count_equal_or_better += 1
            
    p_val = (count_equal_or_better + 1.0) / (num_permutations + 1.0)
    return p_val

def calculate_fdr(p_values):
    n = len(p_values)
    if n == 0:
        return []
    sorted_indices = sorted(range(n), key=lambda i: p_values[i])
    fdr = [1.0] * n
    cum_min = 1.0
    
    for rank_idx in range(n - 1, -1, -1):
        orig_idx = sorted_indices[rank_idx]
        rank = rank_idx + 1
        adj_p = (p_values[orig_idx] * n) / float(rank)
        cum_min = min(cum_min, adj_p)
        fdr[orig_idx] = min(1.0, cum_min)
        
    return fdr

def align_sequences(seq_f_parts, seq_r, seq_f_range, seq_r_range, num_permutations=NUM_PERMUTATIONS, target_side='start', log_func=print):
    log_func(f"Performing boundary-constrained sequence alignment (target_side={target_side}) with Search-Level Permutation Testing...")
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = ALIGN_MATCH_SCORE
    aligner.mismatch_score = ALIGN_MISMATCH_SCORE
    aligner.open_gap_score = ALIGN_OPEN_GAP_SCORE
    aligner.extend_gap_score = ALIGN_EXTEND_GAP_SCORE
    
    raw_results = []
    target_len = len(seq_r)
    
    count_generated = len(seq_f_parts)
    count_aligned = 0
    count_passed_nogap = 0
    count_passed_mismatch = 0
    min_hamming_dist = float('inf')
    
    for idx, (id_f, seq_f, seq_f_start, seq_f_end) in enumerate(seq_f_parts):
        query_len = len(seq_f)
        best_score, best_alignment, best_range = float('-inf'), None, (0, 0)
        
        if target_side == 'start':
            min_pos = 0
            max_pos = min(MAX_BOUNDARY_OFFSET + 1, target_len - query_len + 1)
        else:
            max_pos = target_len - query_len + 1
            min_pos = max(0, target_len - query_len - MAX_BOUNDARY_OFFSET)
            
        for i in range(min_pos, max_pos):
            alignments = aligner.align(seq_f, seq_r[i:i+query_len])
            if alignments:
                score = max(aln.score for aln in alignments)
                if score > best_score:
                    best_score = score
                    best_alignment = alignments[0]
                    best_range = (i, i+query_len)
        
        if best_alignment:
            count_aligned += 1
            aligned_seq_f, aligned_seq_r = best_alignment.sequences
            current_tsd_f = seq_f
            current_tsd_r = seq_r[best_range[0]:best_range[1]]
            
            hamming_dist = calculate_hamming_distance(current_tsd_f, current_tsd_r)
            gaps = aligned_seq_f.count('-') + aligned_seq_r.count('-')
            
            if gaps == 0:
                count_passed_nogap += 1
                if hamming_dist < min_hamming_dist:
                    min_hamming_dist = hamming_dist
                    
                if hamming_dist <= MAX_MISMATCH_ALLOWED:
                    count_passed_mismatch += 1
                    identity_percentage = (len(current_tsd_f) - hamming_dist) / float(len(current_tsd_f)) if len(current_tsd_f) > 0 else 0.0
                    
                    p_val = calculate_empirical_p_value(
                        seq_f=current_tsd_f,
                        seq_r=seq_r,
                        search_max_offset=MAX_BOUNDARY_OFFSET,
                        observed_score=best_score,
                        aligner=aligner,
                        num_permutations=num_permutations,
                        target_side=target_side
                    )
                    
                    raw_results.append({
                        'id_f': id_f,
                        'score': best_score,
                        'tsd_f': current_tsd_f,
                        'tsd_r': current_tsd_r,
                        'f_range': (seq_f_start, seq_f_end),
                        'r_range': (seq_r_range[0] + best_range[0], seq_r_range[0] + best_range[1]),
                        'gaps': gaps,
                        'identity': identity_percentage,
                        'p_value': p_val
                    })
            
        if idx % 50 == 0:
            log_func(f"Aligned and statistically evaluated {idx}/{len(seq_f_parts)} subsequences...")
            
    p_vals = [r['p_value'] for r in raw_results]
    fdr_vals = calculate_fdr(p_vals)
    
    final_results = []
    count_significant = 0
    for r, fdr in zip(raw_results, fdr_vals):
        is_sig = "YES" if fdr < FDR_THRESHOLD else "NO"
        if is_sig == "YES":
            count_significant += 1
        final_results.append((
            r['id_f'], r['score'], r['tsd_f'], r['tsd_r'],
            r['f_range'], r['r_range'], r['gaps'], r['identity'],
            r['p_value'], fdr, is_sig
        ))
        
    min_dist_str = str(min_hamming_dist) if min_hamming_dist != float('inf') else "N/A"
    
    log_func("--- Filtering Statistics Breakdown ---")
    log_func(f"  Candidates Generated               : {count_generated}")
    log_func(f"  Successfully Aligned               : {count_aligned}")
    log_func(f"  Passed No-Gap Filter               : {count_passed_nogap}")
    log_func(f"  Minimum Hamming Distance Observed  : {min_dist_str}")
    log_func(f"  Passed Mismatch Filter (<= {MAX_MISMATCH_ALLOWED} bp)     : {count_passed_mismatch}")
    log_func(f"  Candidates Reaching Permutation Test: {count_passed_mismatch}")
    log_func(f"  Passed Permutation Test            : {count_significant}")
    log_func(f"  Final Reported                     : {len(final_results)}")
    log_func("--------------------------------------")
        
    final_results.sort(key=lambda x: (0 if x[10] == "YES" else 1, x[9], -x[7], -x[1]))
    return final_results

def write_results(file_name, results, ref_range, query_variant, match_percentage, polya_seq, strict_p_val=None):
    print(f"Writing results to {file_name}...")
    
    filtered_results = results
    if strict_p_val is not None:
        if polya_seq == "None":
            filtered_results = []
        else:
            filtered_results = [res for res in results if res[8] < strict_p_val]
            
    with open(file_name, "w") as file:
        file.write(f"Reference Range: {ref_range[0]}-{ref_range[1]}, Query Variant Used: {query_variant}, Max Match Percentage: {match_percentage:.2f}, Detected PolyA: {polya_seq}\n")
        file.write("ID\tScore\tTSD_F\tTSD_R\tTSD_F_Range\tTSD_R_Range\tGaps\tIdentity\tPolyA\tP_value\tFDR\tSignificance\n")
        for res in filtered_results:
            file.write(f"{res[0]}\t{res[1]:.2f}\t{res[2]}\t{res[3]}\t"
                       f"{res[4][0]}-{res[4][1]}\t{res[5][0]}-{res[5][1]}\t{res[6]}\t{res[7]:.2f}\t"
                       f"{polya_seq}\t{res[8]:.4f}\t{res[9]:.4f}\t{res[10]}\n")

def main():
    log_file = open(LOG_FILE_PATH, "w")
    def log_print(msg):
        print(msg)
        log_file.write(msg + "\n")
        log_file.flush()

    log_print("Starting rigorous statistical TSD analysis with Search-Level Permutation & Logging...")
    ref_seq = load_fasta(REF_FASTA_PATH)
    query_seq = load_fasta(QUERY_FASTA_PATH)
    
    start, end, query_variant, match_percentage = sliding_window_match(ref_seq, query_seq)
    ref_range = (start, end)
    seq_f, seq_r, seq_f_range, seq_r_range = extract_flanking_sequences(ref_seq, start, end, length=FLANKING_LENGTH)
    
    log_print(f"Insertion position detected: {start}-{end}, Variant: {query_variant}, Match: {match_percentage:.2f}")
    
    polya_plus, has_polya_plus = detect_strict_polya(ref_seq, start, end, strand='+')
    polya_minus, has_polya_minus = detect_strict_polya(ref_seq, start, end, strand='-')
    
    log_print(f"Poly(A) Plus: {polya_plus} (Found: {has_polya_plus})")
    log_print(f"Poly(A) Minus: {polya_minus} (Found: {has_polya_minus})")
    
    log_print("Processing TSD_plus...")
    seq_f_parts = generate_subsequences(seq_f, seq_f_range[0], min_len=TSD_MIN_LEN, max_len=TSD_MAX_LEN, step=TSD_STEP, side='end')
    results = align_sequences(seq_f_parts, seq_r, seq_f_range, seq_r_range, num_permutations=NUM_PERMUTATIONS, target_side='start', log_func=log_print)
    
    write_results("TSD_plus_full_length.txt", results, ref_range, query_variant, match_percentage, polya_plus)
    write_results("TSD_plus_strict_polyA.txt", results, ref_range, query_variant, match_percentage, polya_plus, strict_p_val=STRICT_P_VAL_THRESHOLD)
    
    log_print("Processing TSD_minus...")
    seq_r_parts = generate_subsequences(seq_r, seq_r_range[0], min_len=TSD_MIN_LEN, max_len=TSD_MAX_LEN, step=TSD_STEP, side='start')
    results_reverse = align_sequences(seq_r_parts, seq_f, seq_r_range, seq_f_range, num_permutations=NUM_PERMUTATIONS, target_side='end', log_func=log_print)
    
    write_results("TSD_minus_full_length.txt", results_reverse, ref_range, query_variant, match_percentage, polya_minus)
    write_results("TSD_minus_strict_polyA.txt", results_reverse, ref_range, query_variant, match_percentage, polya_minus, strict_p_val=STRICT_P_VAL_THRESHOLD)
    
    log_print(f"Statistical TSD Analysis complete! Log saved to {LOG_FILE_PATH}")
    log_file.close()

if __name__ == "__main__":
    main()