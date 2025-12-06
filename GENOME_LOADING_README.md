# Loading KBase Genome Objects into ModelSEEDpy

This guide explains how to load the genome object JSON files from the `genome_objects` directory and use them with ModelSEEDpy for metabolic modeling.

## Overview

The `genome_objects` directory contains 46 KBase-formatted genome objects in JSON format. These genomes have been created from BV-BRC (formerly PATRIC) data and are ready to be used with ModelSEEDpy for metabolic modeling and analysis.

## Files Created

1. **`load_genomes_to_modelseedpy.py`** - Main library with loading utilities
2. **`example_load_genomes.py`** - Executable examples demonstrating usage
3. **`genome_loading_tutorial.ipynb`** - Interactive Jupyter notebook tutorial
4. **`GENOME_LOADING_README.md`** - This documentation file

## Prerequisites

### Required
- Python 3.8+
- The genome object JSON files in `genome_objects/` directory

### Optional (for metabolic modeling)
- ModelSEEDpy: `pip install modelseedpy`
- COBRApy: `pip install cobra`
- pandas: `pip install pandas`
- matplotlib: `pip install matplotlib`

## Quick Start

### 1. List Available Genomes

```python
from load_genomes_to_modelseedpy import list_available_genomes

genomes = list_available_genomes()
print(f"Found {len(genomes)} genomes")
for genome_id in genomes[:10]:
    print(f"  - {genome_id}")
```

### 2. Load a Single Genome

```python
from load_genomes_to_modelseedpy import load_genome

genome = load_genome('Acetobacterium.1')
print(f"Scientific name: {genome['scientific_name']}")
print(f"Features: {len(genome['features'])}")
print(f"DNA size: {genome['dna_size']:,} bp")
```

### 3. Load All Genomes

```python
from load_genomes_to_modelseedpy import load_all_genomes

all_genomes = load_all_genomes()
print(f"Loaded {len(all_genomes)} genomes")
```

### 4. Export Genome Summary

```python
from load_genomes_to_modelseedpy import export_genome_summary

export_genome_summary('genome_summary.tsv')
# Creates a TSV file with summary statistics for all genomes
```

### 5. Create MSGenome Object (Requires ModelSEEDpy)

```python
from load_genomes_to_modelseedpy import get_msgenome

msgenome = get_msgenome('Acetobacterium.1')
# Now you can use msgenome for metabolic modeling
```

### 6. Build a Metabolic Model (Requires ModelSEEDpy)

```python
from load_genomes_to_modelseedpy import build_model

model = build_model('Acetobacterium.1')
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")

# Run FBA
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")
```

## Command-Line Usage

### List all genomes
```bash
python load_genomes_to_modelseedpy.py --list
```

### Load specific genomes
```bash
python load_genomes_to_modelseedpy.py Acetobacterium.1 Bacteroides.1
```

### Export summary table
```bash
python load_genomes_to_modelseedpy.py --summary genome_summary.tsv
```

### Build models for all genomes
```bash
python load_genomes_to_modelseedpy.py --build-models --save-models models/
```

## Running Examples

The `example_load_genomes.py` script contains 6 complete examples:

```bash
# Run all examples interactively
python example_load_genomes.py

# Run a specific example
python example_load_genomes.py 1  # List genomes
python example_load_genomes.py 2  # Load single genome
python example_load_genomes.py 3  # Load all genomes
python example_load_genomes.py 4  # Export summary
python example_load_genomes.py 5  # Create MSGenome
python example_load_genomes.py 6  # Build model
```

## Using the Jupyter Notebook

The interactive tutorial notebook provides a comprehensive walkthrough:

```bash
jupyter notebook genome_loading_tutorial.ipynb
```

The notebook covers:
- Loading and inspecting genomes
- Creating summary tables and visualizations
- Building metabolic models
- Comparing genomes
- Batch processing

## Genome Object Structure

Each genome JSON file contains a KBase Genome object with the following structure:

```json
{
  "id": "Acetobacterium.1",
  "scientific_name": "Bacteria",
  "taxonomy": "Bacteria; Firmicutes; ...",
  "domain": "Bacteria",
  "genetic_code": 11,
  "dna_size": 4200000,
  "num_contigs": 150,
  "contig_ids": ["contig_1", "contig_2", ...],
  "contig_lengths": [28000, 32000, ...],
  "gc_content": 0.45,
  "features": [
    {
      "id": "Acetobacterium.1_0",
      "type": "CDS",
      "functions": ["DNA polymerase III"],
      "dna_sequence": "ATGC...",
      "protein_translation": "MKTL...",
      "location": [[contig_id, start, strand, length]]
    }
  ],
  "cdss": [...],
  "non_coding_features": [...]
}
```

## API Reference

### GenomeLoader Class

```python
from load_genomes_to_modelseedpy import GenomeLoader

loader = GenomeLoader(genome_dir='genome_objects')

# List available genomes
genomes = loader.list_available_genomes()

# Load a genome
genome = loader.load_genome('Acetobacterium.1')

# Load all genomes
all_genomes = loader.load_all_genomes()

# Get MSGenome object (requires ModelSEEDpy)
msgenome = loader.get_msgenome('Acetobacterium.1')

# Build metabolic model (requires ModelSEEDpy)
model = loader.build_model('Acetobacterium.1', template='GramPositive')

# Export summary
loader.export_genome_summary('summary.tsv')
```

### Convenience Functions

```python
# All convenience functions work the same as GenomeLoader methods
# but use a global loader instance

from load_genomes_to_modelseedpy import (
    list_available_genomes,
    load_genome,
    load_all_genomes,
    get_msgenome,
    build_model,
    export_genome_summary
)
```

## Working with ModelSEEDpy

Once you have an MSGenome object, you can use all ModelSEEDpy features:

```python
from load_genomes_to_modelseedpy import get_msgenome
from modelseedpy import MSBuilder, MSGapfill
from modelseedpy.helpers import get_template

# Get genome
msgenome = get_msgenome('Acetobacterium.1')

# Build model
template = get_template('GramPositive')
builder = MSBuilder(msgenome, template)
model = builder.build('acetobacterium_model')

# Gap-fill model
gapfiller = MSGapfill(model)
gapfill_solution = gapfiller.run_gapfilling()

# Test growth
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")
```

## Genome Statistics

To get statistics about all genomes:

```python
from load_genomes_to_modelseedpy import load_all_genomes

genomes = load_all_genomes()

# Calculate statistics
total_features = sum(len(g['features']) for g in genomes.values())
total_dna = sum(g['dna_size'] for g in genomes.values())
avg_gc = sum(g['gc_content'] for g in genomes.values()) / len(genomes)

print(f"Total genomes: {len(genomes)}")
print(f"Total features: {total_features:,}")
print(f"Total DNA: {total_dna/1e9:.2f} Gbp")
print(f"Average GC: {avg_gc:.2%}")

# Group by domain
from collections import Counter
domains = Counter(g['domain'] for g in genomes.values())
print(f"By domain: {dict(domains)}")
```

## Troubleshooting

### ModuleNotFoundError: No module named 'modelseedpy'

Install ModelSEEDpy:
```bash
pip install modelseedpy
```

### FileNotFoundError: Genome file not found

Make sure you're running the script from the repository root directory where the `genome_objects` folder is located.

### Model optimization fails (infeasible)

The model may need gap-filling:
```python
from modelseedpy import MSGapfill

model = build_model('genome_id')
gapfiller = MSGapfill(model)
solution = gapfiller.run_gapfilling()
```

## Advanced Usage

### Custom Genome Directory

```python
from load_genomes_to_modelseedpy import GenomeLoader

loader = GenomeLoader(genome_dir='/path/to/genomes')
genomes = loader.load_all_genomes()
```

### Batch Model Building

```python
from load_genomes_to_modelseedpy import GenomeLoader
from cobra.io import write_sbml_model

loader = GenomeLoader()
genome_ids = loader.list_available_genomes()

for genome_id in genome_ids:
    try:
        model = loader.build_model(genome_id)
        write_sbml_model(model, f'models/{genome_id}_model.xml')
        print(f"✓ {genome_id}")
    except Exception as e:
        print(f"✗ {genome_id}: {e}")
```

### Genome Comparison

```python
from load_genomes_to_modelseedpy import load_genome

# Load genomes to compare
genome1 = load_genome('Acetobacterium.1')
genome2 = load_genome('Bacteroides.1')

# Compare features
functions1 = set()
for feature in genome1['features']:
    functions1.update(feature.get('functions', []))

functions2 = set()
for feature in genome2['features']:
    functions2.update(feature.get('functions', []))

# Find shared and unique functions
shared = functions1 & functions2
unique1 = functions1 - functions2
unique2 = functions2 - functions1

print(f"Shared functions: {len(shared)}")
print(f"Unique to genome1: {len(unique1)}")
print(f"Unique to genome2: {len(unique2)}")
```

## Related Files

- **`bvbrc_to_kbase_genome.py`** - Script used to create these genome objects from BV-BRC data
- **`genome_objects/`** - Directory containing 46 genome JSON files
- **`ASVset_taxonomies/`** - Taxonomic information for synthetic genomes

## Support

For issues or questions:
1. Check this README
2. Review the example scripts
3. Try the Jupyter notebook tutorial
4. See ModelSEEDpy documentation: https://github.com/ModelSEED/ModelSEEDpy

## Citation

If you use these tools in your research, please cite:
- ModelSEEDpy: https://github.com/ModelSEED/ModelSEEDpy
- BV-BRC: https://www.bv-brc.org/
- KBase: https://www.kbase.us/
