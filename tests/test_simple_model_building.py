#!/usr/bin/env python3
"""
Test simple model building approach without MSBuilder
"""

import json
from pathlib import Path

print("=" * 80)
print("Testing Simple Model Building")
print("=" * 80)

# Load genome
genome_file = "genome_objects/AAP99.1.json"
with open(genome_file, 'r') as f:
    genome_data = json.load(f)

genome_id = genome_data['id']
print(f"\n1. Loaded genome: {genome_id}")
print(f"   Features: {len(genome_data.get('features', []))}")

# Try using FBAHelper or cobrakbase to build model
print(f"\n2. Trying different model building approaches...")

# Approach 1: Use cobrakbase KBaseAPI
print("\n  Approach 1: cobrakbase KBaseAPI.fba_tools.genome_to_fbamodel()")
try:
    from cobrakbase.kbase_api import KBaseAPI

    # Create mock KBase object
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

    genome_wrapper = GenomeWrapper(genome_data, genome_id)

    # Use KBaseAPI to build model
    kbase_api = KBaseAPI()
    model = kbase_api.fba_tools.genome_to_fbamodel(genome_wrapper)

    print(f"    ✓ Model created!")
    print(f"    - Reactions: {len(model.reactions)}")
    print(f"    - Metabolites: {len(model.metabolites)}")
    print(f"    - Genes: {len(model.genes)}")

    # Save model
    from cobra.io import write_sbml_model
    output_file = f"{genome_id}_model.xml"
    write_sbml_model(model, output_file)
    print(f"    ✓ Saved to: {output_file}")

    print("\n" + "=" * 80)
    print("SUCCESS! Model built and saved")
    print("=" * 80)

except AttributeError as e:
    print(f"    ✗ Failed: {e}")
    print("    KBaseAPI doesn't have fba_tools.genome_to_fbamodel()")
except Exception as e:
    print(f"    ✗ Failed: {e}")
    import traceback
    traceback.print_exc()

# Approach 2: Create a basic COBRA model manually from genome data
print("\n  Approach 2: Create basic COBRA model from genome features")
try:
    from cobra.core import Model, Reaction, Metabolite, Gene

    model = Model(f"{genome_id}_model")
    model.name = genome_data.get('scientific_name', genome_id)

    # Add genes from features
    for feat in genome_data.get('features', [])[:100]:  # Start with first 100
        gene_id = feat.get('id', '')
        if gene_id:
            gene = Gene(gene_id)
            gene.name = feat.get('function', '')
            model.genes.add(gene)

    print(f"    ✓ Basic model created!")
    print(f"    - Genes: {len(model.genes)}")

    # Note: This is just a skeleton model with no reactions
    # Real model building requires mapping genes to reactions via annotations

except Exception as e:
    print(f"    ✗ Failed: {e}")
    import traceback
    traceback.print_exc()

# Approach 3: Try using ModelSEEDpy FBAHelper
print("\n  Approach 3: ModelSEEDpy FBAHelper")
try:
    from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
    from modelseedpy.core.fbahelper import FBAHelper
    from cobra.core import Model

    # Create empty model
    model = Model(f"{genome_id}_model")

    # Try using FBAHelper
    helper = FBAHelper(model)
    print(f"    FBAHelper created")

    # This approach doesn't work without a properly built model

except Exception as e:
    print(f"    ✗ Failed: {e}")

print("\n" + "=" * 80)
print("Summary:")
print("  - MSBuilder has issues with template biomasses")
print("  - Need alternative approach for model building")
print("  - May need to use different tools or fix template parsing")
print("=" * 80)
