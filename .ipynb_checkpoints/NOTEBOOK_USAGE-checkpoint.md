# Using bvbrc_to_kbase_genome.py in Jupyter Notebooks

This guide shows how to use the genome converter script in Jupyter notebooks for interactive data analysis.

## Quick Start

```python
# Import the module
import sys
sys.path.append('/path/to/GenomeImporter')  # if not in same directory

# Import functions and classes
from bvbrc_to_kbase_genome import (
    fetch_genome_from_api,
    load_genome_from_features,
    load_genome_from_json,
    create_synthetic_genome,
    save_genome_to_json,
    create_fasta_from_genome,
    BVBRCToKBaseConverter,
    LocalGenomeConverter
)
```

## Method 1: Using Convenience Functions

The simplest way to use the script in notebooks is via convenience functions.

### Fetch from BV-BRC API

```python
# Fetch genome from BV-BRC API
genome = fetch_genome_from_api('511145.183')

# Fetch and save to file
genome = fetch_genome_from_api('511145.183', output_file='ecoli.json')

# Access genome data
print(f"Genome: {genome['scientific_name']}")
print(f"Features: {len(genome['features'])}")
print(f"GC content: {genome['gc_content']:.2%}")
```

### Load from Local Features Directory

```python
# Load from features/ and genomes/ directories
genome = load_genome_from_features('511145.183')

# With custom directories and metadata
genome = load_genome_from_features(
    genome_id='511145.183',
    features_dir='path/to/features',
    genomes_dir='path/to/genomes',
    taxonomy='Bacteria; Proteobacteria; Gammaproteobacteria',
    scientific_name='Escherichia coli K-12 MG1655',
    output_file='ecoli.json'
)
```

### Load from Existing JSON

```python
# Load and validate existing genome JSON
genome = load_genome_from_json('my_genome.json')

# Save to new file
save_genome_to_json(genome, 'validated_genome.json')
```

### Create Synthetic Genome

```python
# Create synthetic genome from multiple sources
genome = create_synthetic_genome(
    asv_id='ASV_12345',
    genome_files=['genome1.json', 'genome2.json', 'genome3.json'],
    taxonomy='Bacteria; Firmicutes; Bacilli',
    output_file='synthetic_asv.json'
)

print(f"Synthetic genome has {len(genome['features'])} unique features")
```

### Export to FASTA

```python
# Create FASTA from genome features
create_fasta_from_genome(genome, 'features.fasta')
```

## Method 2: Using Classes Directly

For more control, use the classes directly.

### BVBRCToKBaseConverter Class

```python
from bvbrc_to_kbase_genome import BVBRCToKBaseConverter

# Create converter instance
converter = BVBRCToKBaseConverter('511145.183')

# Fetch metadata
metadata = converter.fetch_genome_metadata()
print(f"Genome name: {metadata['genome_name']}")
print(f"GC content: {metadata['gc_content']}")

# Fetch sequences
sequences = converter.fetch_genome_sequences()
print(f"Number of contigs: {len(sequences)}")

# Fetch features
features = converter.fetch_genome_features()
print(f"Number of features: {len(features)}")

# Build complete genome object
genome = converter.build_kbase_genome()

# Save to file
converter.save_genome(genome, 'ecoli_genome.json')
```

### LocalGenomeConverter Class

```python
from bvbrc_to_kbase_genome import LocalGenomeConverter

# Create converter instance
converter = LocalGenomeConverter()

# Load from features directory
genome = converter.load_genome_from_features_dir(
    genome_id='511145.183',
    features_dir='features',
    genomes_dir='genomes',
    taxonomy='Bacteria; Proteobacteria',
    scientific_name='Escherichia coli'
)

# Create synthetic genome
synthetic = converter.create_synthetic_genome(
    asv_id='ASV_001',
    genome_files=['g1.json', 'g2.json', 'g3.json'],
    taxonomy='Bacteria; Firmicutes'
)

# Validate existing genome
genome = converter.load_genome_json('my_genome.json')
validated = converter.validate_kbase_genome(genome)

# Save genome
converter.save_genome(genome, 'output.json')

# Create FASTA
converter.create_fasta_from_genome(genome, 'features.fasta')
```

## Method 3: Using %run Magic

You can also run the script directly from notebooks using the `%run` magic command.

```python
# Run as script with command-line arguments
%run bvbrc_to_kbase_genome.py --features 511145.183 --output ecoli.json

# Fetch from API
%run bvbrc_to_kbase_genome.py --api 511145.183 --output ecoli.json --fasta ecoli.fasta

# Create synthetic genome
%run bvbrc_to_kbase_genome.py --synthetic ASV_001 --genomes g1.json,g2.json --output synthetic.json
```

## Complete Workflow Examples

### Example 1: Process Multiple Genomes from Features Directory

```python
from bvbrc_to_kbase_genome import load_genome_from_features
from glob import glob
import json

# Get all feature files
feature_files = glob('features/*.json')
genome_ids = [f.split('/')[-1].replace('.json', '') for f in feature_files]

# Process each genome
genomes = {}
for genome_id in genome_ids[:10]:  # Process first 10
    try:
        genome = load_genome_from_features(genome_id)
        genomes[genome_id] = genome
        print(f"✓ Processed {genome_id}: {len(genome['features'])} features")
    except Exception as e:
        print(f"✗ Failed {genome_id}: {e}")

# Analyze results
total_features = sum(len(g['features']) for g in genomes.values())
print(f"\nTotal genomes: {len(genomes)}")
print(f"Total features: {total_features}")
```

### Example 2: Create Synthetic Genomes from ASV Groups

```python
from bvbrc_to_kbase_genome import create_synthetic_genome
from json import load

# Load ASV genome mappings
ASV_genomes = load(open('modeling_files/ASV_genomes.json', 'r'))

# Process each ASV
for asv_id, genome_ids in list(ASV_genomes.items())[:5]:  # First 5 ASVs
    # Prepare genome files
    genome_files = [f'genomes/{gid}.json' for gid in genome_ids]

    # Create synthetic genome
    try:
        synthetic = create_synthetic_genome(
            asv_id=asv_id,
            genome_files=genome_files,
            taxonomy='Bacteria',  # Update with actual taxonomy
            output_file=f'synthetic/{asv_id}.json'
        )
        print(f"✓ Created {asv_id}: {len(synthetic['features'])} unique features")
    except Exception as e:
        print(f"✗ Failed {asv_id}: {e}")
```

### Example 3: Analyze Genome Statistics

```python
from bvbrc_to_kbase_genome import load_genome_from_features
import pandas as pd

# Load genome
genome = load_genome_from_features('511145.183')

# Extract feature statistics
feature_types = {}
for feature in genome['features']:
    ftype = feature['type']
    feature_types[ftype] = feature_types.get(ftype, 0) + 1

# Create DataFrame
df = pd.DataFrame([
    {'Type': k, 'Count': v}
    for k, v in feature_types.items()
])

print(df)

# Genome summary
print(f"\nGenome Summary:")
print(f"  ID: {genome['id']}")
print(f"  Name: {genome['scientific_name']}")
print(f"  DNA size: {genome['dna_size']:,} bp")
print(f"  GC content: {genome['gc_content']:.2%}")
print(f"  Contigs: {genome['num_contigs']}")
print(f"  Features: {len(genome['features'])}")
print(f"  CDS: {len(genome['cdss'])}")
```

### Example 4: Extract Feature Functions

```python
from bvbrc_to_kbase_genome import load_genome_from_features
from collections import Counter

# Load genome
genome = load_genome_from_features('511145.183')

# Extract all functions
all_functions = []
for feature in genome['features']:
    all_functions.extend(feature.get('functions', []))

# Count occurrences
function_counts = Counter(all_functions)

# Top 10 most common functions
print("Top 10 most common functions:")
for func, count in function_counts.most_common(10):
    print(f"  {count:4d} - {func}")
```

### Example 5: Batch Processing with Progress Tracking

```python
from bvbrc_to_kbase_genome import load_genome_from_features
from tqdm import tqdm
import json

# List of genome IDs to process
genome_ids = ['511145.183', '83332.133', '1110693.3']

# Process with progress bar
results = {}
for genome_id in tqdm(genome_ids, desc="Processing genomes"):
    try:
        genome = load_genome_from_features(
            genome_id,
            features_dir='features',
            genomes_dir='genomes'
        )

        results[genome_id] = {
            'name': genome['scientific_name'],
            'features': len(genome['features']),
            'dna_size': genome['dna_size'],
            'gc_content': genome['gc_content']
        }

        # Save individual genome
        with open(f'output/{genome_id}.json', 'w') as f:
            json.dump(genome, f, indent=2)

    except Exception as e:
        results[genome_id] = {'error': str(e)}

# Summary
import pandas as pd
df = pd.DataFrame.from_dict(results, orient='index')
print(df)
```

## Tips for Notebook Usage

1. **Import at the top of your notebook:**
   ```python
   from bvbrc_to_kbase_genome import *
   ```

2. **Use convenience functions for simple tasks:**
   ```python
   genome = load_genome_from_features('511145.183')
   ```

3. **Use classes for complex workflows:**
   ```python
   converter = LocalGenomeConverter()
   genome1 = converter.load_genome_from_features_dir('id1')
   genome2 = converter.load_genome_from_features_dir('id2')
   ```

4. **Access genome data directly:**
   ```python
   # Genome is just a Python dictionary
   print(genome.keys())
   print(genome['scientific_name'])
   for feature in genome['features'][:5]:
       print(feature['id'], feature['type'])
   ```

5. **Combine with pandas for analysis:**
   ```python
   import pandas as pd
   features_df = pd.DataFrame(genome['features'])
   print(features_df.head())
   ```

6. **Handle errors gracefully:**
   ```python
   try:
       genome = load_genome_from_features(genome_id)
   except FileNotFoundError as e:
       print(f"Feature file not found: {e}")
   except Exception as e:
       print(f"Error: {e}")
   ```

## Available Functions Summary

| Function | Purpose |
|----------|---------|
| `fetch_genome_from_api(genome_id)` | Fetch from BV-BRC API |
| `load_genome_from_features(genome_id)` | Load from local features/ dir |
| `load_genome_from_json(json_file)` | Load from genome JSON file |
| `create_synthetic_genome(asv_id, genome_files)` | Merge multiple genomes |
| `save_genome_to_json(genome, output_file)` | Save genome to JSON |
| `create_fasta_from_genome(genome, output_file)` | Export features to FASTA |

## Available Classes Summary

| Class | Purpose |
|-------|---------|
| `BVBRCToKBaseConverter` | Fetch and convert from BV-BRC API |
| `LocalGenomeConverter` | Process local genome files |
