# Bacterial Codon-Usage and Proteome-Wide Analysis Pipeline

This repository contains a Python pipeline for bacterial codon-usage analysis and proteome-wide characterization.

## Repository Contents

- `bacteria.fa` — bacterial coding sequences in FASTA format.
- Python script — performs codon-usage and protein-level analyses.

## Requirements

- Python 3
- Biopython
- NumPy
- pandas
- openpyxl

Install the required packages using:

`pip install biopython numpy pandas openpyxl`

## How to Run

Place `bacteria.fa` in the location expected by the script, open a terminal in the appropriate directory, and run:

`python script_name.py`

Replace `script_name.py` with the actual Python script filename.

## Pipeline Features

The pipeline performs:

- Nucleotide composition and synonymous third-position analysis.
- Codon counting, codon frequencies, and relative synonymous codon usage (RSCU).
- Effective number of codons (Nc), codon adaptation index (CAI), and frequency of optimal codons (Fop).
- Parity analysis and P2.
- RSCU- and relative amino acid usage (RAAU)-based principal component analysis (PCA).
- Protein length, theoretical isoelectric point (pI), molecular weight, GRAVY, and aromaticity.
- Export of per-gene results and summary values to Excel worksheets.

## Publication and Citation

The analytical methods implemented in this pipeline are described in the associated research publication. The final publication details and DOI will be added after publication.

If you use this pipeline in your research, please cite the associated publication.

