#!/usr/bin/env python3
"""
Test the MSGenome creation fix
"""

print("=" * 60)
print("Testing MSGenome Creation")
print("=" * 60)

try:
    from load_genomes_to_modelseedpy import load_genome, get_msgenome, build_model

    # Test 1: Load genome (should work without ModelSEEDpy)
    print("\n1. Loading genome...")
    genome = load_genome('Acetobacterium.1', printing=False)
    print(f"✓ Loaded genome: {genome['id']}")
    print(f"  Features: {len(genome['features'])}")

    # Test 2: Create MSGenome (requires ModelSEEDpy + cobrakbase)
    print("\n2. Creating MSGenome object...")
    try:
        msgenome = get_msgenome('Acetobacterium.1', printing=True)
        print(f"✓ MSGenome created successfully!")
        print(f"  ID: {msgenome.id}")
        print(f"  Features: {len(msgenome.features)}")
        print(f"  Type: {type(msgenome)}")

        # Test 3: Build model (requires ModelSEEDpy)
        print("\n3. Building metabolic model...")
        model = build_model('Acetobacterium.1', printing=True)
        print(f"✓ Model built successfully!")
        print(f"  Reactions: {len(model.reactions)}")
        print(f"  Metabolites: {len(model.metabolites)}")
        print(f"  Genes: {len(model.genes)}")

        # Test 4: Optimize model
        print("\n4. Testing FBA optimization...")
        solution = model.optimize()
        print(f"  Status: {solution.status}")
        if solution.status == 'optimal':
            print(f"  Growth rate: {solution.objective_value:.4f}")

        print("\n" + "=" * 60)
        print("✓ ALL TESTS PASSED!")
        print("=" * 60)

    except ImportError as e:
        print(f"✗ Import Error: {e}")
        print("\nRequired packages:")
        print("  pip install modelseedpy cobrakbase")

    except AttributeError as e:
        print(f"✗ Attribute Error: {e}")
        print("\nThis might be a version compatibility issue.")
        print("Try: pip install --upgrade modelseedpy cobrakbase")

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()

except Exception as e:
    print(f"✗ Fatal Error: {e}")
    import traceback
    traceback.print_exc()
