#!/usr/bin/env python3
"""
Iteratively debug MSGenome creation and model building

This script tests different approaches to create MSGenome objects from
genome JSON files and build COBRA models.
"""

import json
import sys
from pathlib import Path

print("=" * 80)
print("DEBUG: MSGenome Creation and Model Building")
print("=" * 80)

# Load the test genome
genome_file = "genome_objects/AAP99.1.json"
print(f"\n1. Loading test genome: {genome_file}")

if not Path(genome_file).exists():
    print(f"✗ Error: {genome_file} not found")
    print("Available genome files:")
    for f in list(Path("genome_objects").glob("*.json"))[:5]:
        print(f"  - {f}")
    sys.exit(1)

with open(genome_file, 'r') as f:
    genome_data = json.load(f)

genome_id = genome_data['id']
print(f"✓ Loaded genome: {genome_id}")
print(f"  Features: {len(genome_data.get('features', []))}")
print(f"  Scientific name: {genome_data.get('scientific_name', 'N/A')}")

# Check what's available
print("\n2. Checking available packages...")
packages_available = {}

try:
    import modelseedpy
    packages_available['modelseedpy'] = True
    print(f"✓ modelseedpy version: {modelseedpy.__version__}")
except ImportError as e:
    packages_available['modelseedpy'] = False
    print(f"✗ modelseedpy not available: {e}")

try:
    import cobrakbase
    packages_available['cobrakbase'] = True
    print(f"✓ cobrakbase version: {cobrakbase.__version__}")
except ImportError as e:
    packages_available['cobrakbase'] = False
    print(f"✗ cobrakbase not available: {e}")

try:
    import cobra
    packages_available['cobra'] = True
    print(f"✓ cobra version: {cobra.__version__}")
except ImportError as e:
    packages_available['cobra'] = False
    print(f"✗ cobra not available: {e}")

if not all([packages_available['modelseedpy'], packages_available['cobrakbase']]):
    print("\n✗ Required packages not available. Install with:")
    print("  pip install modelseedpy cobrakbase")
    sys.exit(1)

# Approach 1: Check KBaseObjectFactory methods
print("\n3. Exploring KBaseObjectFactory...")
try:
    from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
    factory = KBaseObjectFactory()

    methods = [m for m in dir(factory) if not m.startswith('_')]
    print(f"✓ KBaseObjectFactory methods: {methods[:10]}")

    # Check if any method looks genome-related
    genome_methods = [m for m in methods if 'genome' in m.lower() or 'build' in m.lower()]
    print(f"  Genome-related methods: {genome_methods}")

except Exception as e:
    print(f"✗ Error with KBaseObjectFactory: {e}")

# Approach 2: Try MSGenome directly with different signatures
print("\n4. Testing MSGenome creation approaches...")

from modelseedpy.core.msgenome import MSGenome
import inspect

print(f"  MSGenome.__init__ signature: {inspect.signature(MSGenome.__init__)}")
print(f"  MSGenome class methods:")
class_methods = [m for m in dir(MSGenome) if not m.startswith('_') and callable(getattr(MSGenome, m))]
print(f"    {class_methods[:15]}")

# Approach 3: Check if there's a from_dict or similar
print("\n5. Looking for MSGenome factory methods...")
factory_methods = [m for m in dir(MSGenome) if 'from' in m.lower() or 'create' in m.lower()]
print(f"  Factory-like methods: {factory_methods}")

# Approach 4: Try to create MSGenome by examining actual structure
print("\n6. Attempting different MSGenome creation strategies...")

strategies = []

# Strategy 1: Direct instantiation with no args and set attributes
print("\n  Strategy 1: Empty MSGenome + set attributes")
try:
    msgenome = MSGenome()
    # Try to set genome data attributes
    for key, value in genome_data.items():
        try:
            setattr(msgenome, key, value)
        except:
            pass
    msgenome.id = genome_id
    print(f"    ✓ Created MSGenome, testing attributes...")
    print(f"      - msgenome.id: {msgenome.id}")
    print(f"      - hasattr features: {hasattr(msgenome, 'features')}")
    if hasattr(msgenome, 'features'):
        print(f"      - len(features): {len(msgenome.features) if msgenome.features else 0}")
    strategies.append(("Strategy 1", msgenome))
except Exception as e:
    print(f"    ✗ Failed: {e}")

# Strategy 2: Try using cobrakbase.KBaseAPI
print("\n  Strategy 2: Using cobrakbase components")
try:
    # Maybe we can create a mock workspace object and use it
    from cobrakbase.core.kbasefba.fbamodel import FBAModel

    # Create a simple object with data and info
    class MockWorkspaceObject:
        def __init__(self, data, obj_id):
            self.data = data
            class Info:
                id = obj_id
                name = obj_id
                type = "KBaseGenomes.Genome"
                reference = f"local/{obj_id}"
            self.info = Info()

    mock_obj = MockWorkspaceObject(genome_data, genome_id)

    # Try different factory approaches
    from cobrakbase.core.kbase_object_factory import KBaseObjectFactory
    factory = KBaseObjectFactory()

    # Check what methods are available
    print(f"    Trying factory methods...")

    # Try the most likely candidates
    if hasattr(factory, 'build'):
        try:
            result = factory.build(mock_obj)
            print(f"    ✓ factory.build() worked: {type(result)}")
            strategies.append(("Strategy 2 - factory.build", result))
        except Exception as e:
            print(f"    ✗ factory.build() failed: {e}")

    if hasattr(factory, 'create'):
        try:
            result = factory.create(mock_obj)
            print(f"    ✓ factory.create() worked: {type(result)}")
            strategies.append(("Strategy 2 - factory.create", result))
        except Exception as e:
            print(f"    ✗ factory.create() failed: {e}")

except Exception as e:
    print(f"    ✗ Failed: {e}")
    import traceback
    traceback.print_exc()

# Strategy 3: Bypass MSGenome and go directly to MSBuilder with genome dict
print("\n  Strategy 3: Direct to MSBuilder with genome dict")
try:
    from modelseedpy import MSBuilder
    from modelseedpy.helpers import get_template

    template = get_template('GramNegative')

    # Try passing genome_data directly
    builder = MSBuilder(genome_data, template, model_id=f"{genome_id}_test")
    model = builder.build(f"{genome_id}_test")

    print(f"    ✓ MSBuilder accepted genome dict directly!")
    print(f"      - Reactions: {len(model.reactions)}")
    print(f"      - Metabolites: {len(model.metabolites)}")
    strategies.append(("Strategy 3 - Direct MSBuilder", model))

except Exception as e:
    print(f"    ✗ Failed: {e}")

# Strategy 4: Check MSGenome source for proper initialization
print("\n  Strategy 4: Examine MSGenome initialization requirements")
try:
    # Get the actual signature and parameters
    sig = inspect.signature(MSGenome.__init__)
    params = sig.parameters
    print(f"    MSGenome.__init__ parameters: {list(params.keys())}")

    # Try with different parameter combinations
    if 'data' in params:
        try:
            msgenome = MSGenome(data=genome_data)
            msgenome.id = genome_id
            print(f"    ✓ MSGenome(data=...) worked!")
            strategies.append(("Strategy 4 - data param", msgenome))
        except Exception as e:
            print(f"    ✗ MSGenome(data=...) failed: {e}")

except Exception as e:
    print(f"    ✗ Failed: {e}")

# Test the best strategy
print("\n" + "=" * 80)
print("7. Testing successful strategies...")
print("=" * 80)

if not strategies:
    print("✗ No strategies succeeded. Cannot proceed.")
    sys.exit(1)

for strategy_name, obj in strategies:
    print(f"\nTesting: {strategy_name}")
    try:
        # If it's already a model, just save it
        if hasattr(obj, 'reactions'):
            model = obj
            print(f"  ✓ Already a model")
        else:
            # It's an MSGenome, build a model
            from modelseedpy import MSBuilder
            from modelseedpy.helpers import get_template

            template = get_template('GramNegative')
            builder = MSBuilder(obj, template, model_id=f"{genome_id}_model")
            model = builder.build(f"{genome_id}_model")
            print(f"  ✓ Built model from MSGenome")

        print(f"  Model stats:")
        print(f"    - Reactions: {len(model.reactions)}")
        print(f"    - Metabolites: {len(model.metabolites)}")
        print(f"    - Genes: {len(model.genes)}")

        # Try to optimize
        solution = model.optimize()
        print(f"  FBA optimization:")
        print(f"    - Status: {solution.status}")
        if solution.status == 'optimal':
            print(f"    - Objective: {solution.objective_value:.4f}")

        # Save model
        from cobra.io import write_sbml_model
        output_file = f"{genome_id}_model.xml"
        write_sbml_model(model, output_file)
        print(f"  ✓ Saved to: {output_file}")

        print(f"\n{'=' * 80}")
        print(f"SUCCESS! {strategy_name} works!")
        print(f"{'=' * 80}")
        break

    except Exception as e:
        print(f"  ✗ Failed to build/save model: {e}")
        import traceback
        traceback.print_exc()
        continue
