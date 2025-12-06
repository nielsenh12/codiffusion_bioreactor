#!/usr/bin/env python3
"""
Quick test script to verify genome loading works correctly
"""

print("Testing genome loading...")

# Test 1: Load a genome (doesn't require ModelSEEDpy)
print("\n1. Loading genome...")
from load_genomes_to_modelseedpy import load_genome

genome = load_genome('Acetobacterium.1')
print(f"✓ Loaded genome: {genome['id']}")
print(f"  Features: {len(genome['features'])}")

# Test 2: Try to create MSGenome (requires ModelSEEDpy)
print("\n2. Creating MSGenome object...")
try:
    from load_genomes_to_modelseedpy import get_msgenome

    msgenome = get_msgenome('Acetobacterium.1')
    print(f"✓ MSGenome created: {msgenome.id}")
    print(f"  Features: {len(msgenome.features)}")
    print("\n✓ All tests passed!")

except ImportError as e:
    print(f"⚠ ModelSEEDpy not installed: {e}")
    print("  Install with: pip install modelseedpy")
    print("\n✓ Basic genome loading works (without ModelSEEDpy)")
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
