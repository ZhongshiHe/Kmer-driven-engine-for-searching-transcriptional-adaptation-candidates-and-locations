import os
import csv
import math
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def get_base_name(seq_id):
    base = seq_id.replace("TR_", "").replace("gDNA_", "")
    if base.startswith("PAP"): return base.replace("PAP", "DPP")
    return base

def main():
    parser = argparse.ArgumentParser(description="Map k-mers using strict Hamming distance.")
    parser.add_argument('-k', '--kmer', type=int, required=True, help="K-mer size (e.g., 20)")
    parser.add_argument('-i', '--identity', type=str, required=True, help="Identity threshold (e.g., 80 or 80%)")
    args = parser.parse_args()

    k = args.kmer
    id_str = args.identity.replace('%', '')
    identity = float(id_str) / 100.0

    min_matches = math.ceil(k * identity)
    max_mismatches = k - min_matches
    if max_mismatches < 0: max_mismatches = 0

    print(f"\n=== Hamming Distance MAPPER ===")
    print(f"K-mer size         : {k}")
    print(f"Target Identity    : {identity * 100}%")
    print(f"Allowed Mismatches : {max_mismatches}")
    print(f"=================================\n")

    cds_recs = list(SeqIO.parse("transcripts_renamed.fa", "fasta"))
    gdna_recs = list(SeqIO.parse("gDNA_withUTR.fasta", "fasta"))

    gdna_recs.sort(key=lambda x: get_base_name(x.id))
    cds_recs.sort(key=lambda x: get_base_name(x.id))
    gdna_names = [get_base_name(rec.id) for rec in gdna_recs]

    print("Pre-computing gDNA k-mer coordinate pools (Forward & Reverse)...")
    gdna_pools = {}
    for gdna in gdna_recs:
        base_name = get_base_name(gdna.id)
        fwd_seq = str(gdna.seq).upper().replace("U", "T")
        rev_seq = str(gdna.seq.reverse_complement()).upper().replace("U", "T")
        
        pool = []
        L = len(fwd_seq)
        
        for i in range(L - k + 1):
            pool.append((fwd_seq[i:i+k], i + 1, '+')) 
            
        for i in range(L - k + 1):
            physical_pos = L - i - k + 1 
            pool.append((rev_seq[i:i+k], physical_pos, '-'))
            
        gdna_pools[base_name] = pool

    raw_matrix, norm_matrix = [], []
    map_records = [] 

    for i, cds in enumerate(cds_recs):
        base_name = get_base_name(cds.id)
        print(f"  -> Mapping {i+1}/{len(cds_recs)}: {base_name}")
        
        raw_row, norm_row = [base_name], [base_name]
        
        t_seq = str(cds.seq).upper().replace("U", "T")
        total_kmers = len(t_seq) - k + 1
        if total_kmers <= 0: continue

        for g_name in gdna_names:
            g_kmers_pool = gdna_pools[g_name]
            hits = 0
            
            for t_idx in range(total_kmers):
                tk = t_seq[t_idx : t_idx+k]
                t_pos = t_idx + 1 
                
                matched = False
                for gk_seq, g_pos, g_strand in g_kmers_pool:
                    mismatches = 0
                    
                    for idx in range(k):
                        if tk[idx] != gk_seq[idx]:
                            mismatches += 1
                            if mismatches > max_mismatches:
                                break 
                                
                    if mismatches <= max_mismatches:
                        matched = True
                        map_records.append([
                            base_name, t_pos, tk, 
                            g_name, g_pos, g_strand, mismatches
                        ])
                
                # FIXED: Properly indented to trigger after evaluating EACH k-mer!
                if matched:
                    hits += 1

            raw_row.append(str(hits))
            norm_row.append(str(round((hits / total_kmers) * 100, 2)))

        raw_matrix.append(raw_row)
        norm_matrix.append(norm_row)

    id_label = int(identity * 100)
    raw_out = f"pure_hamming_matrix_k{k}_i{id_label}_RAW.csv"
    norm_out = f"pure_hamming_matrix_k{k}_i{id_label}_NORM.csv"
    map_out = f"pure_hamming_mapping_k{k}_i{id_label}.csv" 

    with open(raw_out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Gene'] + gdna_names)
        writer.writerows(raw_matrix)

    with open(norm_out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Gene'] + [f"{n}_%" for n in gdna_names])
        writer.writerows(norm_matrix)

    with open(map_out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Transcript_Gene', 'Transcript_Pos_1Based', 'Kmer_Sequence', 'Target_gDNA', 'gDNA_Pos_1Based', 'gDNA_Strand', 'Mismatches'])
        writer.writerows(map_records)

    print(f"\n[SUCCESS] Saved RAW Matrix to: {raw_out}")
    print(f"[SUCCESS] Saved NORM Matrix to: {norm_out}")
    print(f"[SUCCESS] Saved Detailed MAPPING TABLE to: {map_out}")

if __name__ == "__main__":
    main()
