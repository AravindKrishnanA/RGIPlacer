#!/bin/bash
#SBATCH --job-name=epa_placement
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=epa_placement.%j.out

source ~/rgienv/bin/activate
export PATH=$HOME/data/aiyer51/2180/software/epa-ng/bin:$PATH

cd /oscar/data/kchellap/aiyer51/2180/

epa-ng \
  --ref-msa card_reference_aligned.clean.fasta \
  --tree card_reference_fasttree.nwk \
  --query ./results/wastewater/query_aligned_wastewater_5k_cleaned.fasta \
  --model LG+G \
  --threads 32 \
  --out-dir epa_placements_wastewater \
  --redo
