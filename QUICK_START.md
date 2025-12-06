# Quick Start Guide: Loading Genomes into ModelSEEDpy

## What Was Created

Three main files to help you work with your genome objects:

1. **`load_genomes_to_modelseedpy.py`** - Main library (can import or run directly)
2. **`example_load_genomes.py`** - 6 runnable examples
3. **`genome_loading_tutorial.ipynb`** - Interactive Jupyter notebook
4. **`GENOME_LOADING_README.md`** - Complete documentation

## Instant Start (Copy & Paste)

### In Python Script or Notebook

```python
from load_genomes_to_modelseedpy import load_genome, load_all_genomes

# Load one genome
genome = load_genome('Acetobacterium.1')
print(f"Loaded: {genome['scientific_name']}")
print(f"Features: {len(genome['features'])}")

# Load all genomes
all_genomes = load_all_genomes()
print(f"Total: {len(all_genomes)} genomes")
```

### From Command Line

```bash
# List all genomes
python3 load_genomes_to_modelseedpy.py --list

# Export summary table
python3 load_genomes_to_modelseedpy.py --summary genomes.tsv

# Run examples
python3 example_load_genomes.py
```

## With ModelSEEDpy (Optional)

First install: `pip install modelseedpy`

```python
from load_genomes_to_modelseedpy import get_msgenome, build_model

# Create MSGenome object
msgenome = get_msgenome('Acetobacterium.1')

# Build metabolic model
model = build_model('Acetobacterium.1')

# Test growth
solution = model.optimize()
print(f"Growth: {solution.objective_value}")
```

## Your Genome Collection

**Location:** `genome_objects/` directory
**Count:** 71 genome objects in KBase JSON format
**Ready to use:** No additional setup needed!

Each genome includes:
- ✓ Full feature annotations
- ✓ DNA and protein sequences
- ✓ Functional assignments
- ✓ Taxonomy information
- ✓ Contig information

## Next Steps

1. **Try the examples:** `python3 example_load_genomes.py`
2. **Read full docs:** Open `GENOME_LOADING_README.md`
3. **Interactive tutorial:** `jupyter notebook genome_loading_tutorial.ipynb`
4. **Build models:** Install ModelSEEDpy and start modeling!

## Available Genomes (First 20)

```
Acetoanaerobium.1
Acetobacterium.1
Acinetobacter.2
Acinetobacter.3
Alcaligenes.14
Aminivibrio.1
Anaerovorax.1
Aquamicrobium.7
Azoarcus.4
Bacillus.10
Bacteroides.1
Brevundimonas.2
Burkholderiales.1
Christensenellaceae_R-7_group.1
Clostridium_sensu_stricto_13.1
Colidextribacter.5
Desulfovibrio.1
Desulfovibrio.2
Dysgonomonas.1
EBM-39.1
... and 51 more
```

Run `python3 load_genomes_to_modelseedpy.py --list` to see all genomes.

## Questions?

- Check `GENOME_LOADING_README.md` for detailed documentation
- Run examples: `python3 example_load_genomes.py`
- ModelSEEDpy docs: https://github.com/ModelSEED/ModelSEEDpy
