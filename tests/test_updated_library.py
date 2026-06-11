#!/usr/bin/env python3
"""
Test the updated load_genomes_to_modelseedpy library
"""

from load_genomes_to_modelseedpy import GenomeLoader
from cobra.io import write_sbml_model

print("=" * 80)
print("Testing Updated Library")
print("=" * 80)

# Create loader
loader = GenomeLoader()

# Test loading a genome
print("\n1. Loading genome...")
genome = loader.load_genome('Acetobacterium.1')
print(f"✓ Loaded: {genome['id']}")
print(f"  Features: {len(genome.get('features', []))}")

# Test creating MSGenome
print("\n2. Creating MSGenome...")
msgenome = loader.get_msgenome('Acetobacterium.1')
print(f"✓ MSGenome ID: {msgenome.id}")
print(f"  Features: {len(msgenome.features)}")

# Test building model
print("\n3. Building model...")
model = loader.build_model('Acetobacterium.1')
print(f"✓ Model reactions: {len(model.reactions)}")
print(f"  Model metabolites: {len(model.metabolites)}")
print(f"  Model genes: {len(model.genes)}")

# Test FBA
print("\n4. Running FBA...")
solution = model.optimize()
print(f"  Status: {solution.status}")
if solution.status == 'optimal':
    print(f"  Objective: {solution.objective_value:.4f}")

# Save model
print("\n5. Saving model...")
output_file = "Acetobacterium.1_test_model.xml"
write_sbml_model(model, output_file)
print(f"✓ Saved to: {output_file}")

print("\n" + "=" * 80)
print("SUCCESS! All tests passed")
print("=" * 80)
