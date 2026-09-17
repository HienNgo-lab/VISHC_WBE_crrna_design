# Target-Enriched CRISPR crRNA Pool Design Pipeline

This repository contains the Python pipeline for automated extraction, coordinate mapping, and construction of target-enriched crRNA oligo pools targeting Antibiotic Resistance Genes (ARGs) from the Comprehensive Antibiotic Resistance Database (CARD). The experimental methodology and bioinformatic workflow align with the CRISPR-NGS target-enrichment framework developed by Mao et al. (2025) [1] and its step-by-step protocol (Mao & Nguyen, 2024) [2], incorporating the FLASH algorithm for guide RNA selection (Quan et al., 2019) [3], CARD (Alcock et al., 2023) [4], and Cas9-targeted engineering protocols (Liang et al., 2015) [5].

---

## 🛠️ Requirements & Installation

Ensure you have Python 3.8+ and `pandas` installed:

`pip install pandas`

---

## 🚀 Usage

Run the design pipeline directly using Python:

`python3 run_design.py`

The script will automatically download the official CARD protein homolog model database, process header formatting, perform dual-strand PAM scanning, and generate detailed output files[1, 2].

---

## 📊 Outputs

The output files will be saved in `CRISPR/crRNA_FLASH/`:

* **`crRNA_oligo_pool.csv`**: Contains `ARO_ID`, `Spacer_20nt`, `Spacer_Start`, `Spacer_End`, `Spacer_Len`, `Strand`, and `Full_oligo_5to3`[1, 4].
* **`crRNA_oligo_pool.fasta`**: FASTA file formatted for direct commercial High-Throughput Oligo Pool ordering (e.g., Twist Bioscience, Agilent)[1, 4].

---

## 🔬 Technical Explanation & Logic Breakdown

### 1. Oligo Architecture (60 nt Total)
Each generated single-stranded DNA (ssDNA) oligo is designed for downstream *in vitro* T7 transcription (IVT)[1, 4]:

Full Oligo (60 nt) = T7 Promoter (18 nt) + Spacer (20 nt) + crRNA Scaffold (22 nt)

* **18 nt T7 Promoter (`TAATACGACTCACTATAG`)**: Transcription initiation site for T7 RNA Polymerase[1, 4].
* **20 nt Spacer**: Target-specific sequence matching the target ARG[1, 4].
* **22 nt crRNA Handle (`GTTTTAGAGCTATGCTGTTTTG`)**: Conserved scaffold region required for tracrRNA annealing and Cas9 complexation[1, 4].

> **PAM Exclusion Note:** The 3 nt PAM sequence (`TTT`) is used for site localization but is strictly excluded from the synthesized oligo[1, 4]. PAM is recognized directly by the Cas9 protein domain on target DNA, not via RNA-DNA base pairing[1, 2].

### 2. Biological Coordinate Indexing
* **Plus Strand (`+`)**: Coordinates use 1-based indexing[1, 4]. For PAM located at 0-based Python index `i`, `Spacer_Start = i + 5` and `Spacer_End = i + 24`[1, 4].
* **Minus Strand (`-`)**: Targets identified on the reverse complement sequence are mathematically mapped back to the reference `(+)` strand coordinate system for seamless alignment in genome visualization tools (e.g., IGV)[1, 4].

### 3. Global Deduplication
Identical 20 nt spacers across gene variants are tracked using a global hash set (`seen_spacers`) to remove duplicates and streamline synthesis costs[1].

---

## 📚 References

1. **Mao, Y., Shisler, J. L., & Nguyen, T. H. (2025).** Enhanced detection for antibiotic resistance genes in wastewater samples using a CRISPR-enriched metagenomic method. *Water Research*, 274, 123056. https://doi.org/10.1016/j.watres.2024.123056.
2. **Mao, Y., & Nguyen, T. H. (2024).** Multiplexed CRISPR-based target-enriched next-generation sequencing for detecting antibiotic resistance genes in environmental samples (Version 2). *protocols.io*. https://dx.doi.org/10.17504/protocols.io.8epv5xdnjg1b/v2.
3. **Quan, J., Langelier, C., Kuchta, A., Batson, J., Teyssier, N., Lyden, A., ... & Crawford, E. D. (2019).** FLASH: a next-generation CRISPR diagnostic for multiplexed detection of antimicrobial resistance sequences. *Nucleic Acids Research*, 47(14), e83. https://doi.org/10.1093/nar/gkz418.
4. **Alcock, B. P., Huynh, W., Chalil, R., Smith, K. W., Raphenya, A. R., Wlodarski, M. A., ... & McArthur, A. G. (2023).** CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1), D690-D699. https://doi.org/10.1093/nar/gkac920.
5. **Liang, X., Potter, J., Kumar, S., Zou, Y., Quintanilla, R., Sridharan, M., ... & Chesnut, J. D. (2015).** Rapid and highly efficient mammalian cell engineering via Cas9 protein transfection. *Journal of Biotechnology*, 208, 44-53. https://doi.org/10.1016/j.jbiotec.2015.04.024.
