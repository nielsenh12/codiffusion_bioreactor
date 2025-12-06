# Printing Parameter Guide

The `printing` parameter has been added to all genome loading and model building functions to give you control over output verbosity.

## Important: Errors Always Print

**Regardless of the `printing` parameter value, errors are always displayed.** This ensures you're never silently missing critical issues.

```python
# Even with printing=False, errors still print
genome = load_genome('InvalidGenome', printing=False)
# ✓ Error will still be displayed

genomes = load_all_genomes(verbose=False, printing=False)
# ✓ Any errors will still be displayed
```

## Quick Reference

All functions now support `printing=True/False`:

```python
from load_genomes_to_modelseedpy import (
    load_genome,
    load_all_genomes,
    get_msgenome,
    build_model
)
```

## Usage Examples

### 1. Load Genome

**With printing (default):**
```python
genome = load_genome('Acetobacterium.1')
# Output:
# Loading genome: Acetobacterium.1
#   Scientific name: Bacteria
#   Taxonomy: Bacteria...
#   Features: 3320
#   DNA size: 0 bp
#   GC content: 42.58%
```

**Without printing:**
```python
genome = load_genome('Acetobacterium.1', printing=False)
# No output - silent operation
```

### 2. Load All Genomes

**Full verbosity (default):**
```python
genomes = load_all_genomes()
# Shows:
#   - Overall progress
#   - Individual genome stats for each
```

**Progress only, no individual stats:**
```python
genomes = load_all_genomes(printing=False)
# Output:
# Loading 459 genomes from genome_objects
# Successfully loaded 459 genomes
```

**Completely silent:**
```python
genomes = load_all_genomes(verbose=False, printing=False)
# No output at all
```

### 3. Get MSGenome

**With printing (default):**
```python
msgenome = get_msgenome('Acetobacterium.1')
# Output:
# Loading genome: Acetobacterium.1
#   Scientific name: Bacteria
#   ...
# Creating MSGenome object for Acetobacterium.1
```

**Without printing:**
```python
msgenome = get_msgenome('Acetobacterium.1', printing=False)
# No output
```

### 4. Build Model

**With printing (default):**
```python
model = build_model('Acetobacterium.1')
# Output:
# Loading genome: Acetobacterium.1
#   Scientific name: Bacteria
#   ...
# Creating MSGenome object for Acetobacterium.1
# Auto-selected template: GramPositive
# Loading template: GramPositive
# Building model: Acetobacterium.1_model
# Model built successfully!
#   Reactions: 1234
#   Metabolites: 1056
#   Genes: 890
```

**Without printing:**
```python
model = build_model('Acetobacterium.1', printing=False)
# No output - returns model silently
```

## Batch Processing Example

When processing many genomes, use `printing=False` to reduce clutter:

```python
from load_genomes_to_modelseedpy import GenomeLoader

loader = GenomeLoader()
genome_ids = loader.list_available_genomes()

print(f"Building models for {len(genome_ids)} genomes...")

models = {}
for i, genome_id in enumerate(genome_ids, 1):
    try:
        # Build without printing - we'll show our own progress
        model = loader.build_model(genome_id, printing=False)
        models[genome_id] = model

        # Custom progress output
        if i % 10 == 0:
            print(f"  Progress: {i}/{len(genome_ids)} models built")

    except Exception as e:
        print(f"  Error with {genome_id}: {e}")

print(f"\n✓ Successfully built {len(models)} models")
```

## Parameter Combinations

### load_genome()
```python
# Default - show everything
genome = load_genome('id')

# Silent
genome = load_genome('id', printing=False)
```

### load_all_genomes()
```python
# Default - show everything
genomes = load_all_genomes()

# Show progress, hide individual stats
genomes = load_all_genomes(printing=False)

# Hide progress, show individual stats
genomes = load_all_genomes(verbose=False)

# Completely silent
genomes = load_all_genomes(verbose=False, printing=False)
```

### get_msgenome()
```python
# Default - show everything
msgenome = get_msgenome('id')

# Silent
msgenome = get_msgenome('id', printing=False)
```

### build_model()
```python
# Default - show everything
model = build_model('id')

# Silent model building
model = build_model('id', printing=False)

# Custom template, silent
model = build_model('id', template='GramNegative', printing=False)
```

## Use Cases

### Use `printing=True` (default) when:
- Working interactively in a notebook
- Debugging genome loading issues
- Exploring data for the first time
- Need to see what's happening

### Use `printing=False` when:
- Batch processing many genomes
- Running automated scripts
- Building pipelines
- Want clean logs with custom messages
- Performance monitoring (reduce I/O overhead)

## Testing

Test the printing parameter:

```bash
python3 test_printing_option.py
```

This will demonstrate all combinations of the printing parameter.

## Quick Comparison

```python
# VERBOSE - Default behavior
genome = load_genome('Acetobacterium.1')
# Loading genome: Acetobacterium.1
#   Scientific name: Bacteria
#   Taxonomy: Bacteria...
#   Features: 3320
#   DNA size: 0 bp
#   GC content: 42.58%

# SILENT - With printing=False
genome = load_genome('Acetobacterium.1', printing=False)
# (no output)
```

## API Summary

All functions now accept `printing` parameter:

| Function | Parameters | Default | What Prints |
|----------|------------|---------|-------------|
| `load_genome()` | `printing=True` | Print stats | Genome stats + errors |
| `load_all_genomes()` | `verbose=True, printing=True` | Print all | Progress, stats + errors |
| `get_msgenome()` | `printing=True` | Print status | Status msgs + errors |
| `build_model()` | `printing=True` | Print progress | Build msgs + errors |
| `export_genome_summary()` | `printing=True` | Print stats | Export msgs + errors |

**Note:** Errors always print regardless of parameter settings.

The `GenomeLoader` class methods have the same parameters.

## Backward Compatibility

All existing code will continue to work as before - the default is `printing=True`, which maintains the original verbose behavior.

```python
# These are equivalent:
genome = load_genome('Acetobacterium.1')
genome = load_genome('Acetobacterium.1', printing=True)
```
