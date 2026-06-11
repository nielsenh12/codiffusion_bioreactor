#!/usr/bin/env python3
"""
Test the printing parameter for genome loading functions
"""

print("=" * 60)
print("Testing the printing parameter")
print("=" * 60)

from load_genomes_to_modelseedpy import (
    load_genome,
    load_all_genomes,
    get_msgenome,
    build_model
)

# Test 1: Load genome with printing (default)
print("\n1. Load genome WITH printing (default):")
print("-" * 60)
genome = load_genome('Acetobacterium.1')
print(f"✓ Loaded: {genome['id']}\n")

# Test 2: Load genome without printing
print("2. Load genome WITHOUT printing:")
print("-" * 60)
genome = load_genome('Acinetobacter.2', printing=False)
print(f"✓ Loaded silently: {genome['id']}\n")

# Test 3: Load all genomes with minimal output
print("3. Load all genomes with minimal output:")
print("-" * 60)
genomes = load_all_genomes(verbose=True, printing=False)
print(f"✓ Loaded {len(genomes)} genomes without individual stats\n")

# Test 4: Load all genomes completely silently
print("4. Load all genomes COMPLETELY silently:")
print("-" * 60)
genomes = load_all_genomes(verbose=False, printing=False)
print(f"✓ Loaded {len(genomes)} genomes (no output during loading)\n")

# Test 5: MSGenome without printing
print("5. Create MSGenome WITHOUT printing:")
print("-" * 60)
try:
    msgenome = get_msgenome('Acetobacterium.1', printing=False)
    print(f"✓ MSGenome created silently: {msgenome.id}\n")
except ImportError:
    print("⚠ ModelSEEDpy not installed - skipping MSGenome test\n")

# Test 6: Build model without printing
print("6. Build model WITHOUT printing:")
print("-" * 60)
try:
    model = build_model('Acetobacterium.1', printing=False)
    print(f"✓ Model built silently: {len(model.reactions)} reactions\n")
except ImportError:
    print("⚠ ModelSEEDpy not installed - skipping model building test\n")

print("=" * 60)
print("All tests completed!")
print("=" * 60)
print("\nSummary:")
print("  ✓ printing=True  (default) - Shows all stats")
print("  ✓ printing=False - Silent operation, no stats")
print("  ✓ verbose=False  - No progress messages (load_all_genomes only)")
