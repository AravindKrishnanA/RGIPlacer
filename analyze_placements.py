import pandas as pd
import json
import numpy as np

with open('../epa_placements/epa_result.jplace') as f:
    jplace = json.load(f)

placements = []
for pquery in jplace['placements']:
    query_name = pquery['n'][0]
    best = pquery['p'][0]
    
    placements.append({
        'query': query_name,
        'LWR': best[2],
        'pendant_length': best[3] if len(best) > 3 else 0
    })

df = pd.DataFrame(placements)

rgi_df = pd.read_csv('../results/gut_microbiome/metagenome_rgi_1.txt', sep='\t')

print("RGI columns available:")
print(list(rgi_df.columns))
print(f"\nFirst row Best_Identities: {rgi_df['Best_Identities'].iloc[0] if 'Best_Identities' in rgi_df else 'Column not found'}")

df = df.merge(
    rgi_df[['ORF_ID', 'Best_Hit_ARO', 'AMR Gene Family', 'Drug Class', 
            'Cut_Off', 'Best_Identities', 'Best_Hit_Bitscore']],
    left_on='query', right_on='ORF_ID', how='left'
)

def parse_identity(x):
    if pd.isna(x) or x == '':
        return 0
    try:
        if '/' in str(x):
            parts = str(x).split('/')
            return (float(parts[0]) / float(parts[1])) * 100
        else:
            return float(x)
    except:
        return 0

df['percent_identity'] = df['Best_Identities'].apply(parse_identity)

print(f"\nPercent identity after parsing:")
print(f"  Min: {df['percent_identity'].min():.2f}%")
print(f"  Median: {df['percent_identity'].median():.2f}%")
print(f"  Max: {df['percent_identity'].max():.2f}%")
print(f"  Non-zero: {(df['percent_identity'] > 0).sum()}")

median_pendant = df['pendant_length'].median()

df['novel_candidate'] = (
    (df['LWR'] < 0.7) &                      # Poor placement confidence
    (df['pendant_length'] > median_pendant)   # Long evolutionary distance
)

df.to_csv('../results/gut_microbiome/placement_analysis_FIXED.csv', index=False)

novel = df[df['novel_candidate']].copy()
novel = novel.sort_values('LWR')  

print(f"\n{'='*70}")
print("NOVEL MECHANISM CANDIDATES (PHYLOGENETIC CRITERIA)")
print(f"{'='*70}")
print(f"Total sequences: {len(df)}")
print(f"LWR < 0.7: {(df['LWR'] < 0.7).sum()}")
print(f"Pendant > median ({median_pendant:.4f}): {(df['pendant_length'] > median_pendant).sum()}")
print(f"BOTH criteria (novel candidates): {len(novel)}")

if len(novel) > 0:
    novel.to_csv('../results/gut_microbiome/novel_candidates_FIXED.csv', index=False)
    
    print(f"\nTop 30 novel candidates (sorted by worst LWR):")
    print(novel[['query', 'Best_Hit_ARO', 'AMR Gene Family', 'Drug Class',
                 'LWR', 'pendant_length', 'Cut_Off']].head(30).to_string(index=False))

    print(f"\nSaved to: ../results/wastewater/novel_candidates_FIXED.csv")
else:
    print("\nNo candidates with BOTH criteria. Most extreme by LWR only:")
    extreme = df[df['LWR'] < 0.5].sort_values('LWR')
    print(extreme[['query', 'Best_Hit_ARO', 'AMR Gene Family', 'LWR', 
                   'pendant_length']].head(30).to_string(index=False))
