import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os

rgi_tsv = "../results/wastewater/metagenome_rgi.txt"
out_dir = "../results/wastewater"

os.makedirs(out_dir, exist_ok=True)

print(f"Reading RGI results from {rgi_tsv}")
df = pd.read_csv(rgi_tsv, sep="\t")

df = df[df["Predicted_Protein"].notna()].copy()
print(f"Total hits with protein sequences: {len(df)}")

df = df[df["Cut_Off"].isin(["Loose", "Strict"])].copy()
print(f"Loose+Strict hits: {len(df)}")

records = []
for _, row in df.iterrows():
    seq = str(row["Predicted_Protein"]).strip()
    if not seq:
        continue
    rid = str(row["ORF_ID"])
    desc = f"{row['Best_Hit_ARO']}|{row['AMR Gene Family']}|{row['Drug Class']}|{row['Cut_Off']}"
    records.append(SeqRecord(Seq(seq), id=rid, description=desc))

fasta_all = os.path.join(out_dir, "rgi_hits_all.faa")
SeqIO.write(records, fasta_all, "fasta")
print(f"Wrote {len(records)} sequences to {fasta_all}")

f70 = os.path.join(out_dir, "rgi_hits_70.faa")
print("\n[CD-HIT] 70% identity clustering...")
subprocess.run([
    "cd-hit",
    "-i", fasta_all,
    "-o", f70,
    "-c", "0.70",
    "-n", "5",
    "-M", "16000",
    "-T", "16",
    "-d", "0"
], check=True)

f90 = os.path.join(out_dir, "rgi_hits_90.faa")
print("\n[CD-HIT] 90% identity clustering on 70% reps...")
subprocess.run([
    "cd-hit",
    "-i", f70,
    "-o", f90,
    "-c", "0.90",
    "-n", "5",
    "-M", "16000",
    "-T", "16",
    "-d", "0"
], check=True)

print("\nClustering finished.")
print(f"All sequences: {fasta_all}")
print(f"70% cluster reps: {f70}")
print(f"90% cluster reps (final reps): {f90}")
