import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


novel_df = pd.read_csv('../results/wastewater/novel_candidates_FIXED.csv')

# Load jplaimport pandas as 
# Load jplace file with all placement data
with open('../epa_placements_wastewater/epa_result.jplace', 'r') as f:
    jplace = json.load(f)

print(f"Loaded {len(novel_df)} novel candidates")
print(f"Loaded jplace file with {len(jplace['placements'])} query sequences")



def analyze_placement_quality(query_name, jplace_data, top_n=10):
    #Analyze if placement distribution is flat (bad) or peaked (good)

    
    # Find this query in jplace data
    placement_data = None
    for pquery in jplace_data['placements']:
        if query_name in pquery['n'][0]:
            placement_data = pquery
            break
    
    if placement_data is None:
        return {
            'top_lwr': np.nan,
            'top_10_lwr_sum': np.nan,
            'lwr_std': np.nan,
            'log_likelihood_range': np.nan,
            'is_flat': True,
            'quality_score': 0,
            'status': 'NOT_FOUND'
        }
    
    # Extract placements (format: [edge_num, log_likelihood, LWR, pendant_length, ...])
    placements = placement_data['p']
    
    placements_sorted = sorted(placements, key=lambda x: x[2], reverse=True)
    
    top_placements = placements_sorted[:min(top_n, len(placements_sorted))]
    
    # Extract metrics
    top_lwr = top_placements[0][2] if len(top_placements) > 0 else 0
    top_10_lwrs = [p[2] for p in top_placements]
    top_10_lwr_sum = sum(top_10_lwrs)
    lwr_std = np.std(top_10_lwrs) if len(top_10_lwrs) > 1 else 0
    
    # Log-likelihood range (should be large if one placement is clearly better)
    top_log_likelihoods = [p[1] for p in top_placements]
    log_likelihood_range = max(top_log_likelihoods) - min(top_log_likelihoods) if len(top_log_likelihoods) > 1 else 0

    is_flat = (
        (top_lwr < 0.15) or
        (lwr_std < 0.01 and top_lwr < 0.3) or
        (log_likelihood_range < 1.0 and top_lwr < 0.3)
    )
    

    quality_score = 0
    
    # Component 1: Top LWR contribution (0-40 points)
    quality_score += min(40, top_lwr * 100)
    
    # Component 2: LWR concentration in top placement (0-30 points)
    lwr_concentration = top_lwr / top_10_lwr_sum if top_10_lwr_sum > 0 else 0
    quality_score += lwr_concentration * 30
    
    # Component 3: Log-likelihood separation (0-30 points)
    ll_score = min(30, log_likelihood_range * 3)
    quality_score += ll_score
    
    return {
        'top_lwr': top_lwr,
        'top_10_lwr_sum': top_10_lwr_sum,
        'lwr_std': lwr_std,
        'log_likelihood_range': log_likelihood_range,
        'is_flat': is_flat,
        'quality_score': quality_score,
        'status': 'ANALYZED'
    }


print("ANALYZING PLACEMENT QUALITY FOR ALL CANDIDATES")

results = []
for idx, row in novel_df.iterrows():
    query = row['query']
    
    if idx % 10 == 0:
        print(f"Processing {idx+1}/{len(novel_df)}...")
    
    analysis = analyze_placement_quality(query, jplace)
    
    results.append({
        'query': query,
        **analysis
    })

results_df = pd.DataFrame(results)

novel_df_enhanced = novel_df.merge(results_df, on='query', how='left')



def categorize_candidate(row):    
    if pd.isna(row['quality_score']):
        return 'UNKNOWN'
    
    if row['is_flat']:
        return 'LIKELY_FALSE_POSITIVE'
    
    if row['quality_score'] >= 50:
        return 'HIGH_CONFIDENCE'
    elif row['quality_score'] >= 30:
        return 'MODERATE_CONFIDENCE'
    else:
        return 'LOW_CONFIDENCE'

novel_df_enhanced['confidence_category'] = novel_df_enhanced.apply(categorize_candidate, axis=1)


print("PLACEMENT QUALITY ANALYSIS RESULTS")

category_counts = novel_df_enhanced['confidence_category'].value_counts()
print("\nConfidence Categories:")
for cat, count in category_counts.items():
    pct = 100 * count / len(novel_df_enhanced)
    print(f"  {cat}: {count} ({pct:.1f}%)")

print("\n" + "="*70)
print("FLAT DISTRIBUTION EXAMPLES (Likely False Positives)")
print("="*70)

flat_candidates = novel_df_enhanced[novel_df_enhanced['is_flat']].sort_values('LWR')
print(f"\nFound {len(flat_candidates)} candidates with flat placement distributions")
print("\nTop 10 most suspicious (lowest quality scores):")
print(flat_candidates[['query', 'Best_Hit_ARO', 'LWR', 
                       'top_lwr', 'lwr_std', 'log_likelihood_range', 
                       'quality_score']].head(10).to_string(index=False))

print("\n" + "="*70)
print("HIGH-CONFIDENCE CANDIDATES (Peaked Distributions)")
print("="*70)

high_conf = novel_df_enhanced[novel_df_enhanced['confidence_category']=='HIGH_CONFIDENCE'].sort_values('quality_score', ascending=False)
print(f"\nFound {len(high_conf)} high-confidence novel candidates")
if len(high_conf) > 0:
    print("\nTop 10 highest quality:")
    print(high_conf[['query', 'Best_Hit_ARO', 'LWR', 
                     'top_lwr', 'log_likelihood_range', 
                     'quality_score']].head(10).to_string(index=False))



fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Quality score vs LWR
ax1 = axes[0, 0]
colors = {'HIGH_CONFIDENCE': '#10b981', 'MODERATE_CONFIDENCE': '#f59e0b', 
          'LOW_CONFIDENCE': '#ef4444', 'LIKELY_FALSE_POSITIVE': '#dc2626'}

for cat in novel_df_enhanced['confidence_category'].unique():
    subset = novel_df_enhanced[novel_df_enhanced['confidence_category'] == cat]
    ax1.scatter(subset['LWR'], subset['quality_score'], 
               label=cat, alpha=0.7, s=50, c=colors.get(cat, 'gray'))

ax1.set_xlabel('Original LWR', fontweight='bold')
ax1.set_ylabel('Placement Quality Score', fontweight='bold')
ax1.set_title('A. Quality Score vs LWR', fontweight='bold')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)
ax1.axhline(y=30, color='red', linestyle='--', alpha=0.5)
ax1.axhline(y=50, color='green', linestyle='--', alpha=0.5)

ax2 = axes[0, 1]
flat = novel_df_enhanced[novel_df_enhanced['is_flat']]
not_flat = novel_df_enhanced[~novel_df_enhanced['is_flat']]

ax2.hist([flat['log_likelihood_range'].dropna(), not_flat['log_likelihood_range'].dropna()],
         bins=20, label=['Flat (False Positive)', 'Peaked (Real Candidate)'],
         color=['#ef4444', '#10b981'], alpha=0.7, edgecolor='black')
ax2.set_xlabel('Log-Likelihood Range (Top 10 Placements)', fontweight='bold')
ax2.set_ylabel('Frequency', fontweight='bold')
ax2.set_title('B. Log-Likelihood Separation', fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)

ax3 = axes[1, 0]
ax3.hist([flat['lwr_std'].dropna(), not_flat['lwr_std'].dropna()],
         bins=20, label=['Flat', 'Peaked'],
         color=['#ef4444', '#10b981'], alpha=0.7, edgecolor='black')
ax3.set_xlabel('LWR Standard Deviation (Top 10)', fontweight='bold')
ax3.set_ylabel('Frequency', fontweight='bold')
ax3.set_title('C. LWR Variability', fontweight='bold')
ax3.legend()
ax3.grid(True, alpha=0.3)

ax4 = axes[1, 1]
category_data = novel_df_enhanced['confidence_category'].value_counts()
colors_pie = [colors.get(cat, 'gray') for cat in category_data.index]
wedges, texts, autotexts = ax4.pie(category_data.values, labels=category_data.index,
                                   autopct='%1.1f%%', colors=colors_pie,
                                   textprops={'fontweight': 'bold'})
ax4.set_title('D. Candidate Classification', fontweight='bold')

plt.tight_layout()
plt.savefig('../results/figures/placement_quality_analysis_wastewater.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved visualization to figures/placement_quality_analysis.png")


novel_df_enhanced.to_csv('../results/wastewater/novel_candidates_quality_filtered.csv', index=False)
print(f"Saved enhanced data to novel_candidates_quality_filtered.csv")

high_quality = novel_df_enhanced[novel_df_enhanced['confidence_category'].isin(['HIGH_CONFIDENCE', 'MODERATE_CONFIDENCE'])]
high_quality.to_csv('../results/wastewater/novel_candidates_HIGH_QUALITY.csv', index=False)
print(f"Saved {len(high_quality)} high-quality candidates to novel_candidates_HIGH_QUALITY.csv")




moderate_conf_count = len(novel_df_enhanced[novel_df_enhanced['confidence_category']=='MODERATE_CONFIDENCE'])


if len(high_quality) > 0:
    print("TOP 5 CANDIDATES FOR EXPERIMENTAL VALIDATION")
    
    top_5 = high_quality.nlargest(5, 'quality_score')
    print(top_5[['query', 'Best_Hit_ARO', 'AMR Gene Family', 'LWR', 
                 'pendant_length', 'quality_score']].to_string(index=False))

print("\nAnalysis complete")
