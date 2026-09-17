import os
import re
import sys
import glob
import subprocess
import pandas as pd

print("=======================================================")
print("CRRNA DESIGN PIPELINE (YUQING MAO PROTOCOL)")
print("=======================================================\n")

# Set up target output directory
target_dir = os.path.join("CRISPR", "FLASH_crRNA")
os.makedirs(target_dir, exist_ok=True)

raw_fasta = os.path.join(target_dir, "nucleotide_fasta_protein_homolog_model.fasta")
trimmed_fasta = os.path.join(target_dir, "card_homolog_trimmed_yuqing.fasta")

# ---------------------------------------------------------
# STEP 0: UNIT TEST - VERIFY SPACER EXTRACTION & COORDINATES
# ---------------------------------------------------------
print("🔎 Running Unit Test for spacer extraction and coordinate mapping...")
test_seq = "ATGCGTTTNGATCGATCGATCGATCGATC"
i_test = test_seq.find("TTT")
test_spacer = test_seq[i_test+4:i_test+24]
test_start = i_test + 5
test_end = i_test + 24

assert test_spacer == "GATCGATCGATCGATCGATC", "Unit Test Error: Extracted spacer sequence is incorrect!"
assert len(test_spacer) == 20, "Unit Test Error: Spacer length is not 20 nt!"
assert test_start == 10, "Unit Test Error: Spacer Start coordinate mismatch!"
assert test_end == 29, "Unit Test Error: Spacer End coordinate mismatch!"
print("✅ Unit Test Passed: Biological 1-based coordinate logic verified.\n")

# 1. DOWNLOAD & EXTRACT OFFICIAL CARD DATABASE IF NOT PRESENT
if not os.path.exists(raw_fasta):
    print("📥 Downloading official CARD Protein Homolog Model database...")
    tar_path = os.path.join(target_dir, "card-data.tar.gz")
    
    cmd_download = f"wget --no-check-certificate https://card.mcmaster.ca/latest/data -O {tar_path}"
    subprocess.run(cmd_download, shell=True, check=True)
    
    print("📦 Extracting CARD data archive...")
    cmd_extract = f"tar -xvf {tar_path} -C {target_dir}"
    subprocess.run(cmd_extract, shell=True, check=True)
    
    matches = glob.glob(f"{target_dir}/**/nucleotide_fasta_protein_homolog_model.fasta", recursive=True)
    if matches:
        raw_fasta = matches[0]

input_file = raw_fasta if os.path.exists(raw_fasta) else os.path.join(target_dir, "card_homolog_trimmed.fasta")
print(f"📄 Processing source file: {input_file}")

# 2. HEADER STANDARDIZATION (STEP 5 OF YUQING MAO PROTOCOL)
# Extract ARO accession number (e.g., >3002070) to clean FASTA headers
with open(input_file, "r", encoding="utf-8", errors="ignore") as f:
    raw_content = f.readlines()

output_content = []
for line in raw_content:
    line_str = line.strip()
    if line_str.startswith(">"):
        aro_match = re.findall(r'ARO:[0-9]+', line_str)
        if aro_match:
            output_content.append(f">{aro_match[0].split(':')[1]}")
        else:
            output_content.append(line_str.split()[0])
    else:
        output_content.append(line_str)

with open(trimmed_fasta, "w", encoding="utf-8") as f_out:
    f_out.write("\n".join(output_content))

# 3. DUAL-STRAND PAM SCANNING & COORDINATE MAPPING
def rev_comp(seq):
    """Generates reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

output_csv = os.path.join(target_dir, "crRNA_oligo_pool_full_detailed.csv")
output_fasta = os.path.join(target_dir, "crRNA_oligo_pool_full_detailed.fasta")

with open(trimmed_fasta, "r", encoding="utf-8") as f:
    content = f.read().split(">")

total_genes = len(content) - 1
print(f"📊 Total CARD gene entries detected: {total_genes}")

oligos = []
seen_spacers = set()
gene_count = 0

# Fixed functional components for IVT oligo synthesis
t7_promoter = "TAATACGACTCACTATAG"      # 18 nt T7 promoter
crrna_handle = "GTTTTAGAGCTATGCTGTTTTG" # 22 nt crRNA scaffold handle

for block in content:
    if not block.strip():
        continue
    gene_count += 1
    lines = block.strip().split("\n")
    header = lines[0]
    seq = "".join(lines[1:]).upper()
    seq_len = len(seq)

    # --- Scan Plus Strand (+) ---
    for i in range(seq_len - 23):
        if seq[i:i+3] == "TTT": # Cas12a/Cas9 PAM motif (TTTN / TTTV)
            sp = seq[i+4:i+24]  # Extract 20 nt spacer
            if len(sp) == 20 and "N" not in sp and sp not in seen_spacers:
                seen_spacers.add(sp)
                start_pos = i + 5  # 1-based start coordinate
                end_pos = i + 24    # 1-based end coordinate
                full_oligo = f"{t7_promoter}{sp}{crrna_handle}"
                oligos.append((header, sp, start_pos, end_pos, len(sp), "+", full_oligo))

    # --- Scan Minus Strand (-) ---
    rc_seq = rev_comp(seq)
    for i in range(seq_len - 23):
        if rc_seq[i:i+3] == "TTT":
            sp = rc_seq[i+4:i+24]
            if len(sp) == 20 and "N" not in sp and sp not in seen_spacers:
                seen_spacers.add(sp)
                # Map coordinates back to reference (+) strand
                orig_start = seq_len - (i + 24) + 1
                orig_end = seq_len - (i + 4)
                full_oligo = f"{t7_promoter}{sp}{crrna_handle}"
                oligos.append((header, sp, orig_start, orig_end, len(sp), "-", full_oligo))

    if gene_count % 1000 == 0 or gene_count == total_genes:
        print(f"⏳ Progress: {gene_count}/{total_genes} genes processed | Collected: {len(oligos)} crRNA spacers...")

# 4. EXPORT TABULAR & FASTA RESULTS
cols = ["ARO_ID", "Spacer_20nt", "Spacer_Start", "Spacer_End", "Spacer_Len", "Strand", "Full_oligo_5to3"]
df = pd.DataFrame(oligos, columns=cols)
df.to_csv(output_csv, index=False)

with open(output_fasta, "w", encoding="utf-8") as f_out:
    for idx, row in df.iterrows():
        f_out.write(f">crRNA_{idx+1}_ARO_{row['ARO_ID']}_{row['Spacer_Start']}-{row['Spacer_End']}({row['Strand']})\n{row['Full_oligo_5to3']}\n")

print("\n=======================================================")
print("🎉 CRRNA ENRICHMENT POOL DESIGN COMPLETE!")
print(f"CARD Genes Processed  : {total_genes}")
print(f"Unique Spacers Found  : {len(oligos)}")
print(f"CSV Output File       : {output_csv}")
print(f"FASTA Output File     : {output_fasta}")
print("=======================================================\n")