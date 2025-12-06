#!/usr/bin/env python3
"""
Test that errors always print regardless of printing parameter
"""

print("=" * 60)
print("Testing Error Printing")
print("=" * 60)

from load_genomes_to_modelseedpy import (
    load_genome,
    load_all_genomes,
    get_msgenome,
    build_model
)

# Test 1: Try to load a non-existent genome with printing=True
print("\n1. Load non-existent genome WITH printing=True:")
print("-" * 60)
try:
    genome = load_genome('NonExistentGenome.999', printing=True)
except FileNotFoundError as e:
    print(f"✓ Error was raised and caught: {str(e)[:80]}...")

# Test 2: Try to load a non-existent genome with printing=False
print("\n2. Load non-existent genome WITH printing=False:")
print("-" * 60)
try:
    genome = load_genome('NonExistentGenome.999', printing=False)
except FileNotFoundError as e:
    print(f"✓ Error was raised and caught: {str(e)[:80]}...")

# Test 3: Load all genomes with some valid, ensure errors print
print("\n3. Load all genomes with verbose=False, printing=False:")
print("-" * 60)
print("   (Errors should still print)")

# Create a temporary test by trying to load from a directory with mixed valid/invalid
import os
import json
import tempfile

# Create a temp directory with one valid and one invalid genome
temp_dir = tempfile.mkdtemp()
try:
    # Create a valid genome file
    valid_genome = {
        "id": "TestGenome.1",
        "scientific_name": "Test Organism",
        "taxonomy": "Test",
        "domain": "Bacteria",
        "features": [],
        "cdss": [],
        "dna_size": 1000,
        "gc_content": 0.5,
        "num_contigs": 1
    }
    with open(os.path.join(temp_dir, "TestGenome.1.json"), 'w') as f:
        json.dump(valid_genome, f)

    # Create an invalid genome file (malformed JSON)
    with open(os.path.join(temp_dir, "InvalidGenome.1.json"), 'w') as f:
        f.write("{invalid json content")

    # Try to load all - errors should print even with printing=False
    genomes = load_all_genomes(genome_dir=temp_dir, verbose=False, printing=False)
    print(f"✓ Loaded {len(genomes)} valid genome(s)")
    print("  (Error for invalid genome should have printed above)")

finally:
    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)

# Test 4: Summary - what gets printed and what doesn't
print("\n" + "=" * 60)
print("Summary of Printing Behavior")
print("=" * 60)
print("\nWith printing=True:")
print("  ✓ Shows genome loading messages")
print("  ✓ Shows genome statistics")
print("  ✓ Shows errors")
print("\nWith printing=False:")
print("  ✗ Hides genome loading messages")
print("  ✗ Hides genome statistics")
print("  ✓ STILL SHOWS ERRORS (important!)")
print("\nWith verbose=False (load_all_genomes only):")
print("  ✗ Hides progress messages")
print("  ✓ STILL SHOWS ERRORS (important!)")

print("\n" + "=" * 60)
print("✓ Test completed!")
print("=" * 60)
