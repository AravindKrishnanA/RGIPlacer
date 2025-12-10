from Bio import SeqIO

query_ids = set(rec.id for rec in SeqIO.parse("/users/aiyer51/data/aiyer51/2180/results/wastewater/rgi_hits_70_5k.faa", "fasta"))

query_aligned = [rec for rec in SeqIO.parse("/users/aiyer51/data/aiyer51/2180/results/combined_5k_wastewater.fasta", "fasta") 
                 if rec.id in query_ids]

SeqIO.write(query_aligned, "/users/aiyer51/data/aiyer51/2180/results/wastewater/query_aligned_wastewater_5k.fasta", "fasta")
print(f"Extracted {len(query_aligned)} aligned query sequences")