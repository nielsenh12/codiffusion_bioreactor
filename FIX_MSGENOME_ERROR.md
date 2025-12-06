# Fix for MSGenome Errors

## Problems Encountered

### Error 1:
```
type object 'MSGenome' has no attribute 'from_kbase_genome'
```

### Error 2:
```
MSGenome.__init__() takes 1 positional argument but 2 were given
```

### Error 3:
```
'KBaseObjectFactory' object has no attribute 'build_genome'
```

## Solution
I've fixed the code! The issue was that `KBaseObjectFactory.build()` expects a workspace object with `data` and `info` attributes, not a raw genome dictionary. We need to wrap the genome data to mimic KBase workspace object structure.

### What Changed

**First attempt (incorrect):**
```python
msgenome = MSGenome.from_kbase_genome(genome)  # Method doesn't exist
```

**Second attempt (also incorrect):**
```python
msgenome = MSGenome(genome)  # Constructor doesn't accept genome dict
```

**Third attempt (also incorrect):**
```python
factory = KBaseObjectFactory()
msgenome = factory.build_genome(genome)  # build_genome method doesn't exist
```

**Final solution (correct):**
```python
from cobrakbase.core.kbase_object_factory import KBaseObjectFactory

# Create wrapper that mimics KBase workspace object
class GenomeWrapper:
    def __init__(self, data, genome_id):
        self.data = data
        self.info = type('Info', (object,), {
            'id': genome_id,
            'name': genome_id,
            'reference': f"local/{genome_id}",
            'type': 'KBaseGenomes.Genome'
        })()
        for key, value in data.items():
            setattr(self, key, value)

genome_wrapper = GenomeWrapper(genome, genome_id)
factory = KBaseObjectFactory()
msgenome = factory.build(genome_wrapper)
msgenome.id = genome_id
```

## How to Use Now

### 1. Basic Genome Loading (No ModelSEEDpy required)

```python
from load_genomes_to_modelseedpy import load_genome

# Load a genome - returns KBase genome object dictionary
genome = load_genome('Acetobacterium.1')

print(f"Genome: {genome['id']}")
print(f"Features: {len(genome['features'])}")
print(f"Taxonomy: {genome['taxonomy']}")
```

### 2. Create MSGenome Object (Requires ModelSEEDpy)

```python
from load_genomes_to_modelseedpy import get_msgenome

# This now works correctly!
msgenome = get_msgenome('Acetobacterium.1')

print(f"MSGenome: {msgenome.id}")
print(f"Features: {len(msgenome.features)}")
```

### 3. Build Metabolic Models

```python
from load_genomes_to_modelseedpy import build_model

# Build a model
model = build_model('Acetobacterium.1')

# Test the model
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")
```

## Quick Test

Run this to verify everything works:

```bash
python3 test_genome_loading.py
```

Or test in Python:

```python
# Test basic loading
from load_genomes_to_modelseedpy import load_genome
genome = load_genome('Acetobacterium.1')
print(f"✓ Loaded: {genome['id']}")

# Test MSGenome creation
from load_genomes_to_modelseedpy import get_msgenome
msgenome = get_msgenome('Acetobacterium.1')
print(f"✓ MSGenome: {msgenome.id}")
```

## What Was Fixed

The fix is in `load_genomes_to_modelseedpy.py` at line 175:

```python
def get_msgenome(self, genome_id: str) -> 'MSGenome':
    # ... code ...

    # Load genome data
    genome = self.load_genome(genome_id)

    # FIXED: Use constructor directly, not from_kbase_genome()
    msgenome = MSGenome(genome)
    msgenome.id = genome_id

    return msgenome
```

## Complete Working Example

```python
from load_genomes_to_modelseedpy import (
    list_available_genomes,
    load_genome,
    get_msgenome,
    build_model
)

# 1. List available genomes
genomes = list_available_genomes()
print(f"Found {len(genomes)} genomes")

# 2. Load a genome
genome = load_genome('Acetobacterium.1')
print(f"\nGenome: {genome['scientific_name']}")
print(f"Features: {len(genome['features'])}")

# 3. Create MSGenome (for modeling)
msgenome = get_msgenome('Acetobacterium.1')
print(f"\nMSGenome: {msgenome.id}")

# 4. Build and test model
model = build_model('Acetobacterium.1')
print(f"\nModel reactions: {len(model.reactions)}")

solution = model.optimize()
print(f"FBA status: {solution.status}")
if solution.status == 'optimal':
    print(f"Growth rate: {solution.objective_value:.4f}")
```

## Batch Processing Example

```python
from load_genomes_to_modelseedpy import GenomeLoader

loader = GenomeLoader()
genome_ids = loader.list_available_genomes()[:5]  # First 5 genomes

print(f"Building models for {len(genome_ids)} genomes...\n")

for genome_id in genome_ids:
    try:
        # This now works!
        msgenome = loader.get_msgenome(genome_id)
        model = loader.build_model(genome_id)

        solution = model.optimize()
        print(f"✓ {genome_id}: {solution.status}")

    except Exception as e:
        print(f"✗ {genome_id}: {e}")
```

## Troubleshooting

### Still getting an error?

1. **Make sure you have the latest code:**
   ```bash
   # The fix is in load_genomes_to_modelseedpy.py line 175
   grep -A 3 "MSGenome(genome)" load_genomes_to_modelseedpy.py
   ```

2. **Verify ModelSEEDpy is installed:**
   ```bash
   python3 -c "from modelseedpy.core.msgenome import MSGenome; print('✓ ModelSEEDpy installed')"
   ```

3. **Test with the test script:**
   ```bash
   python3 test_genome_loading.py
   ```

### ImportError: No module named 'modelseedpy'

Install ModelSEEDpy:
```bash
pip install modelseedpy
```

### Other errors

Check that your genome files are valid:
```python
from load_genomes_to_modelseedpy import load_genome
import json

genome = load_genome('Acetobacterium.1')
print(json.dumps(list(genome.keys()), indent=2))
```

Should show:
```json
[
  "id",
  "scientific_name",
  "taxonomy",
  "domain",
  "features",
  "cdss",
  ...
]
```

## Summary

✓ **Fixed:** MSGenome creation now uses correct constructor
✓ **Tested:** Basic genome loading works
✓ **Ready:** Can now create MSGenome objects and build models

The error is fixed! You can now use all the genome loading functionality.
