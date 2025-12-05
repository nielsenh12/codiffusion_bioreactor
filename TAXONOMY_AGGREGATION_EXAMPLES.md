# Taxonomy Aggregation Examples

This guide demonstrates how to use the taxonomy aggregation feature when creating synthetic genomes.

## Overview

When creating synthetic genomes from multiple source genomes, the script can:
1. **Extract** taxonomies from all source genomes
2. **Aggregate** them by taxonomic level (Kingdom, Phylum, Class, Order, Family, Genus, Species)
3. **Select** the most prevalent taxonomy at each level as the consensus
4. **Save** a JSON file showing all taxonomies from the source genomes
5. **Apply** the consensus taxonomy to the synthetic genome

## Output Format

The taxonomy aggregation creates a JSON file in `ASVset_taxonomies/{asv_id}.json`:

```json
{
  "Kingdom": ["Bacteria", "Bacteria", "Bacteria"],
  "Phylum": ["Proteobacteria", "Firmicutes", "Proteobacteria"],
  "Class": ["Gammaproteobacteria", "Bacilli", "Gammaproteobacteria"],
  "Order": ["Enterobacterales", "Lactobacillales", "Enterobacterales"],
  "Family": ["Enterobacteriaceae", "Streptococcaceae", "Enterobacteriaceae"],
  "Genus": ["Escherichia", "Streptococcus", "Escherichia"],
  "Species": ["coli", "thermophilus", "coli"]
}
```

In this example:
- **Consensus**: "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia; coli"
- Most prevalent at each level is selected (e.g., "Proteobacteria" appears 2/3 times)

## Command-Line Usage

### Basic Usage (Auto-aggregation)

**Option 1: From Genome IDs (features/{id}.json and genomes/{id}.fna)**

```bash
# Taxonomy automatically aggregated from source genomes in features/ and genomes/
python bvbrc_to_kbase_genome.py --synthetic ASV_12345 \
  --genome-ids 511145.183,83332.133,1110693.3 \
  --output asv_12345.json

# Creates:
# - asv_12345.json (synthetic genome with consensus taxonomy)
# - ASVset_taxonomies/ASV_12345.json (taxonomy aggregation)
```

**Option 2: From Genome JSON Files**

```bash
# From pre-saved genome JSON files
python bvbrc_to_kbase_genome.py --synthetic ASV_12345 \
  --genomes genome1.json,genome2.json,genome3.json \
  --output asv_12345.json

# Creates:
# - asv_12345.json (synthetic genome with consensus taxonomy)
# - ASVset_taxonomies/ASV_12345.json (taxonomy aggregation)
```

### Custom Taxonomy Directory

```bash
python bvbrc_to_kbase_genome.py --synthetic ASV_12345 \
  --genomes genome1.json,genome2.json,genome3.json \
  --taxonomy-dir my_taxonomies \
  --output asv_12345.json

# Creates taxonomy file at: my_taxonomies/ASV_12345.json
```

### Disable Taxonomy Aggregation

```bash
python bvbrc_to_kbase_genome.py --synthetic ASV_12345 \
  --genomes genome1.json,genome2.json,genome3.json \
  --no-taxonomy \
  --output asv_12345.json

# No taxonomy file created, synthetic genome uses default or provided taxonomy
```

### Explicit Taxonomy Override

```bash
# Provide explicit taxonomy (still saves aggregation for reference)
python bvbrc_to_kbase_genome.py --synthetic ASV_12345 \
  --genomes genome1.json,genome2.json,genome3.json \
  --taxonomy "Bacteria; Firmicutes; Bacilli; Lactobacillales" \
  --output asv_12345.json

# Uses provided taxonomy instead of consensus
# Still creates ASVset_taxonomies/ASV_12345.json showing source taxonomies
```

## Notebook Usage

### Method 1: Using Convenience Function (from Genome IDs - Recommended)

```python
from bvbrc_to_kbase_genome import create_synthetic_genome

# Auto-aggregation from genome IDs (features/{id}.json and genomes/{id}.fna)
genome = create_synthetic_genome(
    asv_id='ASV_12345',
    genome_ids=['511145.183', '83332.133', '1110693.3'],
    features_dir='features',
    genomes_dir='genomes',
    output_file='asv_12345.json'
)

# Taxonomy saved to: ASVset_taxonomies/ASV_12345.json
# Genome uses consensus taxonomy
print(f"Consensus taxonomy: {genome['taxonomy']}")
```

### Method 1b: Using Convenience Function (from JSON Files)

```python
from bvbrc_to_kbase_genome import create_synthetic_genome

# From pre-saved genome JSON files
genome = create_synthetic_genome(
    asv_id='ASV_12345',
    genome_files=['genome1.json', 'genome2.json', 'genome3.json'],
    output_file='asv_12345.json'
)

# Taxonomy saved to: ASVset_taxonomies/ASV_12345.json
print(f"Consensus taxonomy: {genome['taxonomy']}")
```

### Method 2: Manual Taxonomy Aggregation

```python
from bvbrc_to_kbase_genome import (
    load_genome_from_json,
    aggregate_taxonomies,
    create_synthetic_genome
)

# Load source genomes
genomes = [
    load_genome_from_json('genome1.json'),
    load_genome_from_json('genome2.json'),
    load_genome_from_json('genome3.json')
]

# Aggregate taxonomies separately
consensus_taxonomy, taxonomy_dict = aggregate_taxonomies(
    genomes=genomes,
    asv_id='ASV_12345',
    output_dir='my_taxonomies'
)

print(f"Consensus: {consensus_taxonomy}")
print(f"Kingdom distribution: {taxonomy_dict['Kingdom']}")
print(f"Phylum distribution: {taxonomy_dict['Phylum']}")

# Create synthetic genome with consensus
genome = create_synthetic_genome(
    asv_id='ASV_12345',
    genome_files=['genome1.json', 'genome2.json', 'genome3.json'],
    taxonomy=consensus_taxonomy,
    save_taxonomy=False  # Already saved above
)
```

### Method 3: Using Class Directly

```python
from bvbrc_to_kbase_genome import LocalGenomeConverter

converter = LocalGenomeConverter()

# Load genomes
genomes = [
    converter.load_genome_json('genome1.json'),
    converter.load_genome_json('genome2.json'),
    converter.load_genome_json('genome3.json')
]

# Aggregate taxonomies
consensus, tax_dict = converter.aggregate_taxonomies(
    genomes=genomes,
    asv_id='ASV_12345',
    output_dir='ASVset_taxonomies'
)

# Analyze taxonomy distribution
from collections import Counter
for level, values in tax_dict.items():
    counts = Counter(values)
    print(f"\n{level}:")
    for taxon, count in counts.items():
        print(f"  {taxon}: {count}/{len(values)} ({count/len(values):.1%})")
```

## Complete Workflow Example

### Batch Process ASVs from Notebook Data (Recommended Workflow)

```python
from bvbrc_to_kbase_genome import create_synthetic_genome
from json import load, dump
from tqdm import tqdm
import os

# Load ASV mappings
ASV_genomes = load(open('modeling_files/ASV_genomes.json', 'r'))
md5_ID = load(open('model_inputs/md5_ID.json', 'r'))

# Process each ASV
taxonomy_summary = {}

for ASV_md5, genome_md5s in tqdm(list(ASV_genomes.items())[:10]):
    asv_id = ASV_md5  # or use iterativeIDs mapping

    # Get genome IDs from md5 mappings
    genome_ids = [md5_ID[md5] for md5 in genome_md5s if md5 in md5_ID]

    if not genome_ids:
        print(f"  Skipping {asv_id}: No valid genome IDs")
        continue

    # Create synthetic genome directly from genome IDs
    # Loads from features/{genome_id}.json and genomes/{genome_id}.fna
    try:
        os.makedirs('synthetic_genomes', exist_ok=True)

        synthetic = create_synthetic_genome(
            asv_id=asv_id,
            genome_ids=genome_ids,  # Load directly from features/ and genomes/
            features_dir='features',
            genomes_dir='genomes',
            save_taxonomy=True,
            taxonomy_output_dir='ASVset_taxonomies',
            output_file=f'synthetic_genomes/{asv_id}.json'
        )

        # Load and save taxonomy summary
        tax_file = f'ASVset_taxonomies/{asv_id}.json'
        if os.path.exists(tax_file):
            with open(tax_file, 'r') as f:
                tax_data = load(f)
                taxonomy_summary[asv_id] = {
                    'consensus': synthetic['taxonomy'],
                    'num_sources': len(genome_ids),
                    'phylum_diversity': len(set(tax_data.get('Phylum', []))),
                    'genus_diversity': len(set(tax_data.get('Genus', [])))
                }

        print(f"✓ {asv_id}: {len(synthetic['features'])} features, taxonomy: {synthetic['taxonomy']}")

    except Exception as e:
        print(f"✗ {asv_id}: {e}")

# Save summary
with open('taxonomy_summary.json', 'w') as f:
    dump(taxonomy_summary, f, indent=2)

print(f"\nProcessed {len(taxonomy_summary)} ASVs")
```

## Analyzing Taxonomy Aggregations

### Load and Analyze Taxonomy Files

```python
from json import load
from glob import glob
from collections import Counter
import pandas as pd

# Load all taxonomy files
taxonomy_files = glob('ASVset_taxonomies/*.json')

results = []
for tax_file in taxonomy_files:
    asv_id = tax_file.split('/')[-1].replace('.json', '')

    with open(tax_file, 'r') as f:
        tax_data = load(f)

    # Calculate diversity at each level
    row = {'ASV_ID': asv_id}
    for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']:
        if level in tax_data:
            unique = len(set(tax_data[level]))
            total = len(tax_data[level])
            row[f'{level}_diversity'] = unique
            row[f'{level}_total'] = total
            row[f'{level}_consensus'] = Counter(tax_data[level]).most_common(1)[0][0]

    results.append(row)

# Create DataFrame
df = pd.DataFrame(results)
print(df.head())

# Summary statistics
print("\nTaxonomic Diversity Summary:")
print(df[['Phylum_diversity', 'Genus_diversity']].describe())
```

### Visualize Taxonomy Distribution

```python
import matplotlib.pyplot as plt
from json import load

# Load taxonomy for one ASV
asv_id = 'ASV_12345'
with open(f'ASVset_taxonomies/{asv_id}.json', 'r') as f:
    tax_data = load(f)

# Plot phylum distribution
from collections import Counter
phylum_counts = Counter(tax_data['Phylum'])

plt.figure(figsize=(10, 6))
plt.bar(phylum_counts.keys(), phylum_counts.values())
plt.xlabel('Phylum')
plt.ylabel('Count')
plt.title(f'Phylum Distribution in {asv_id}')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f'{asv_id}_phylum_distribution.png')
plt.show()
```

## Output File Structure

After running taxonomy aggregation, you'll have:

```
ASVset_taxonomies/
├── ASV_12345.json
├── ASV_67890.json
└── ...

synthetic_genomes/
├── ASV_12345.json  (genome with consensus taxonomy)
├── ASV_67890.json
└── ...
```

Each taxonomy JSON shows the raw taxonomic data from all source genomes, allowing you to:
- Verify consensus accuracy
- Assess taxonomic diversity within ASV
- Track source genome contributions
- Debug classification issues

## Tips

1. **Check consensus quality**: Load the taxonomy JSON to see how consistent taxonomies are
   ```python
   tax_data = load(open('ASVset_taxonomies/ASV_12345.json'))
   print(f"Phylum consistency: {len(set(tax_data['Phylum']))}/{len(tax_data['Phylum'])}")
   ```

2. **Override when needed**: If consensus is wrong due to outliers, provide explicit taxonomy
   ```python
   genome = create_synthetic_genome(asv_id, files,
                                    taxonomy="Bacteria; Proteobacteria; ...")
   ```

3. **Disable for speed**: If you don't need taxonomy tracking, disable it
   ```python
   genome = create_synthetic_genome(asv_id, files, save_taxonomy=False)
   ```

4. **Batch analysis**: Analyze all taxonomy files together to understand your dataset
   ```python
   tax_files = glob('ASVset_taxonomies/*.json')
   # Analyze diversity, consistency, common phyla, etc.
   ```
