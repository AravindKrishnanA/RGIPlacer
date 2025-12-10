# RGIPlacer
RGI-Placer is an analysis framework that automatically identifies novel candidate antimicrobial resistance(AMR) hits from the BLAST-based Resistance Gene Identifier(RGI) tool through the use of phylogenetic placement. 

# Metagenomic AMR Discovery Pipeline

A comprehensive workflow for identifying antimicrobial resistance (AMR) genes in human gut metagenomes using homology-based detection and phylogenetic placement.

## Overview

This pipeline performs the following steps:

1. **Read Download** - Retrieve metagenomic sequencing data
2. **Assembly** - Assemble reads into contigs
3. **ORF Prediction** - Identify open reading frames
4. **CARD/RGI Annotation** - Detect AMR genes using the CARD database
5. **Clustering** - Reduce redundancy in detected sequences
6. **Reference Tree Construction** - Build phylogenetic tree from CARD homologs
7. **Phylogenetic Placement** - Place query sequences onto reference tree
8. **Candidate Identification** - Prioritize novel AMR mechanisms

## Prerequisites

### Required Software

- **SRA Toolkit** - For downloading sequencing data
- **MEGAHIT** - Metagenomic assembler
- **Prodigal** - ORF prediction tool
- **RGI (Resistance Gene Identifier)** - CARD database AMR annotation
- **BLAST+** - Sequence alignment
- **MAFFT** - Multiple sequence alignment
- **FastTree** or **RAxML-NG** - Phylogenetic tree inference
- **CD-HIT** - Sequence clustering
- **pplacer** or **EPA-ng** - Phylogenetic placement

### Environment Setup

```bash
# Create project directory
mkdir -p ~/amr_discovery
cd ~/amr_discovery
```

## Pipeline Steps

### 1. Download Human Gut Metagenome Reads

Download paired-end reads from NCBI SRA. This example uses the first 500k reads as a test subset.

```bash
mkdir -p data/raw_fastq
module load sra-toolkit

# Download first 500k reads (paired-end)
fastq-dump --split-files \
  --maxSpotId 500000 \
  --outdir data/raw_fastq \
  SRR341549
```

**Outputs:**
- `data/raw_fastq/SRR341549_1.fastq`
- `data/raw_fastq/SRR341549_2.fastq`

---

### 2. Assemble Reads with MEGAHIT

#### Build MEGAHIT (one-time setup)

```bash
cd $HOME/data/aiyer51/2180/software
git clone https://github.com/voutcn/megahit.git
cd megahit
git submodule update --init
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
make simple_test

# Add to PATH
export PATH=$HOME/data/aiyer51/2180/software/megahit/build:$PATH
```

#### Run Assembly

```bash
cd /oscar/data/kchellap/aiyer51/2180
mkdir -p data/assembled_out

megahit \
  -1 data/raw_fastq/SRR341549_1.fastq \
  -2 data/raw_fastq/SRR341549_2.fastq \
  -o data/assembled_out
```

**Output:**
- `data/assembled_out/final.contigs.fa` - Assembled contigs

---

### 3. Predict ORFs with Prodigal

Run Prodigal in metagenome mode and clean stop codons (RGI rejects sequences with `*` characters).

```bash
cd /oscar/data/kchellap/aiyer51/2180
module load prodigal

prodigal \
  -i data/assembled_out/final.contigs.fa \
  -a data/assembled_out/final.contigs.proteins.faa \
  -d data/assembled_out/final.contigs.orfs.fna \
  -p meta

# Remove stop codons that cause RGI rejections
sed 's/\*//g' data/assembled_out/final.contigs.proteins.faa \
  > data/assembled_out/final.contigs.proteins.clean.faa
```

**Outputs:**
- `data/assembled_out/final.contigs.proteins.faa` - Predicted proteins (raw)
- `data/assembled_out/final.contigs.proteins.clean.faa` - Cleaned proteins
- `data/assembled_out/final.contigs.orfs.fna` - Nucleotide ORF sequences

---

### 4. Detect AMR Genes with RGI (CARD)

Run RGI on cleaned protein sequences using BLAST mode. Ensure CARD database is loaded and BLAST is available.

```bash
cd /oscar/data/kchellap/aiyer51/2180
source rgienv/bin/activate

# Verify CARD database version
rgi database --version

module load blast
mkdir -p results/gut_microbiome

rgi main \
  -i data/assembled_out/final.contigs.proteins.clean.faa \
  -o results/gut_microbiome/metagenome_rgi_1 \
  -t protein \
  -a BLAST \
  --include_loose --low_quality \
  -n 16 \
  --clean
```

#### Verify RGI Output

```bash
# Count total hits
wc -l results/gut_microbiome/metagenome_rgi_1.txt

# Preview results
head -n 20 results/gut_microbiome/metagenome_rgi_1.txt

# Count strict and loose hits
grep -c "Strict" results/gut_microbiome/metagenome_rgi_1.txt
grep -c "Loose"  results/gut_microbiome/metagenome_rgi_1.txt
```

**Output:**
- `results/gut_microbiome/metagenome_rgi_1.txt` - RGI tabular results

---

### 5. Cluster RGI Hits

Cluster RGI hit protein sequences to reduce redundancy using CD-HIT at 70% and 90% identity thresholds.

```bash
# Run clustering script
python scripts/cluster_rgi_hits.py
```

**Outputs:**
- `results/gut_microbiome/rgi_hits_70.faa` - Representative sequences (70% identity)
- `results/gut_microbiome/rgi_hits_90.faa` - Representative sequences (90% identity)

---

### 6. Build CARD Reference Phylogenetic Tree

#### Download CARD Data

```bash
wget -O card_data.tar.bz2 https://card.mcmaster.ca/latest/data
tar -xjf card_data.tar.bz2
```

#### Align with MAFFT

```bash
mafft --auto --thread 48 \
  protein_fasta_protein_homolog_model.fasta \
  > card_reference_aligned.fasta

# Verify alignment
grep -c "^>" card_reference_aligned.fasta

# Sanitize alignment
python SanitizeMafft.py
```

#### Infer Tree with FastTree

```bash
FastTree -wag -gamma card_reference_aligned.clean.fasta \
  > card_reference_fasttree.nwk
```

**Alternative: RAxML-NG**
```bash
raxml-ng --all --msa card_reference_aligned.clean.fasta \
  --model LG+G8+F --threads 48 --prefix card_reference
```

**Outputs:**
- `card_reference_aligned.clean.fasta` - Sanitized alignment
- `card_reference_fasttree.nwk` - Reference phylogenetic tree

---

### 7. Phylogenetic Placement

Add clustered query sequences to the reference alignment and place them onto the reference tree.

```bash
# Add queries to reference alignment
mafft --add results/gut_microbiome/rgi_hits_70.faa \
  --keeplength \
  --thread 16 \
  card_reference_aligned.clean.fasta \
  > combined_alignment.fasta

# Preprocess and run placement
python query_aligned.py
sbatch epa.sh

# Analyze placement results
python analyze_placements.py
```

**Outputs:**
- `.jplace` or `.epa` files - Placement results
- Analysis tables summarizing placement statistics

---

### 8. Identify Novel Mechanism Candidates

Filter placements and annotation metadata to prioritize novel or suspicious AMR candidates.

#### Filtering Criteria

- Non-flat LWR (likelihood weight ratio) distributions across the tree
- Placements in subtrees lacking matching CARD annotations
- High placement confidence but low sequence identity to annotated homologs
- Consistent evidence across multiple clustering thresholds

```bash
# Example analysis command
python analyze_placements.py \
  --placements results/placements.jplace \
  --rgi results/gut_microbiome/metagenome_rgi_1.txt \
  --out results/novel_candidates.tsv
```

**Output:**
- `results/novel_candidates.tsv` - Prioritized candidate AMR sequences

---

## Directory Structure

```
amr_discovery/
├── data/
│   ├── raw_fastq/              # Downloaded FASTQ files
│   └── assembled_out/          # Assembly and ORF prediction outputs
├── results/
│   └── gut_microbiome/         # RGI results and clustering outputs
├── scripts/                    # Analysis scripts
├── card_reference_aligned.clean.fasta
├── card_reference_fasttree.nwk
└── README.md
```

---

## Notes & Best Practices

- **Path Configuration**: Replace example paths (e.g., `/oscar/data/kchellap/aiyer51/2180`) with your environment-specific paths
- **Workflow Automation**: Consider using Snakemake or Nextflow for reproducibility and HPC integration
- **Parallelization**: Parallelize compute-intensive steps (MAFFT, MEGAHIT, RGI) using available cluster resources
- **Version Control**: Log versions of all tools used for reproducibility

---

## Suggested Improvements

1. **Workflow Manager**: Implement using Snakemake or Nextflow with step dependencies
2. **Version Logging**: Document versions of MEGAHIT, Prodigal, MAFFT, RGI, FastTree/RAxML-NG, CD-HIT
3. **Testing**: Include a small test dataset with `--quick-test` mode for CI/CD
4. **Reporting**: Generate summary reports (MultiQC-style) and interactive placement visualizations (iTOL, gappa)
5. **Error Handling**: Add logging and automatic retries for failed steps

---

## References

- **CARD Database**: [https://card.mcmaster.ca](https://card.mcmaster.ca)
- **RGI Documentation**: [https://github.com/arpcard/rgi](https://github.com/arpcard/rgi)
- **MEGAHIT**: [https://github.com/voutcn/megahit](https://github.com/voutcn/megahit)
- **Prodigal**: [https://github.com/hyattpd/Prodigal](https://github.com/hyattpd/Prodigal)
- **MAFFT**: [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/)

---


## Citation

If you use this pipeline, please cite:
- The CARD database
- RGI
- All tools used in your analysis

---

## Contact

aravind_iyer@brown.edu
