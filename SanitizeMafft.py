from Bio import SeqIO

in_f = "card_reference_aligned.fasta"
out_f = "card_reference_aligned.clean.fasta"

clean_records = []
for i, rec in enumerate(SeqIO.parse(in_f, "fasta"), start=1):
    safe_id = rec.id
    for ch in [' ', ':', ',', '(', ')', "'", ';', '|']:
        safe_id = safe_id.replace(ch, '_')
    rec.id = safe_id
    rec.name = safe_id
    rec.description = safe_id  
    clean_records.append(rec)

SeqIO.write(clean_records, out_f, "fasta")
print(f"Wrote {len(clean_records)} records to {out_f}")