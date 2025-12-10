#!/bin/bash
#SBATCH --job-name=card_tree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        
#SBATCH --time=24:00:00
#SBATCH --mem=60G
#SBATCH --output=card_tree.%j.out
#SBATCH --error=card_tree.%j.err

module load gcc

cd /oscar/data/kchellap/aiyer51/2180/

raxml-ng --msa  card_reference_aligned.clean.fasta \
         --model LG+G8+F \
         --threads ${SLURM_CPUS_PER_TASK} \
         --seed 12345 \
         --tree pars{1} \
         --prefix card_reference_clean_t${SLURM_CPUS_PER_TASK}
