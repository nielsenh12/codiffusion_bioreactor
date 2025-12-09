#!/usr/bin/env python3
"""
Integration Example: Using Genome Loader with Existing Code

This example shows how to integrate the genome loader with your existing
bioreactor modeling workflow.
"""

import sys
import os

# Import existing utilities
# (Adjust these imports based on your actual setup)
try:
    from util import CommUtil
    UTIL_AVAILABLE = True
except ImportError:
    print("Note: util.py not available, skipping CommUtil examples")
    UTIL_AVAILABLE = False

# Import genome loader
from load_genomes_to_modelseedpy import (
    GenomeLoader,
    load_genome,
    load_all_genomes,
    get_msgenome,
    build_model
)


def example_1_basic_integration():
    """
    Example 1: Load genomes and use with ModelSEEDpy
    """
    print("=" * 60)
    print("Example 1: Basic Integration with ModelSEEDpy")
    print("=" * 60)

    # Create loader
    loader = GenomeLoader(genome_dir='genome_objects')

    # Get available genomes
    genome_ids = loader.list_available_genomes()
    print(f"\nFound {len(genome_ids)} genomes")

    # Load a few genomes
    print("\nLoading first 3 genomes...")
    for genome_id in genome_ids[:3]:
        genome = loader.load_genome(genome_id)
        print(f"  {genome_id}: {len(genome['features'])} features")

    print("\nGenomes loaded successfully!")
    print()


def example_2_build_models_batch():
    """
    Example 2: Build metabolic models for multiple genomes
    """
    print("=" * 60)
    print("Example 2: Batch Model Building")
    print("=" * 60)

    try:
        # Import ModelSEEDpy components
        from modelseedpy import MSBuilder
        from cobra.io import write_sbml_model

        loader = GenomeLoader()
        genome_ids = loader.list_available_genomes()[:5]  # First 5 for demo

        print(f"\nBuilding models for {len(genome_ids)} genomes...\n")

        models = {}
        for genome_id in genome_ids:
            try:
                # Build model
                model = loader.build_model(genome_id)
                models[genome_id] = model

                # Test optimization
                solution = model.optimize()
                print(f"✓ {genome_id}")
                print(f"  Reactions: {len(model.reactions)}")
                print(f"  Status: {solution.status}")
                if solution.status == 'optimal':
                    print(f"  Growth: {solution.objective_value:.4f}")
                print()

            except Exception as e:
                print(f"✗ {genome_id}: {e}\n")
                continue

        print(f"Successfully built {len(models)} models")
        print()

    except ImportError:
        print("\nModelSEEDpy not installed. Install with:")
        print("  pip install modelseedpy\n")


def example_3_integrate_with_util():
    """
    Example 3: Integrate with existing CommUtil workflow
    """
    print("=" * 60)
    print("Example 3: Integration with CommUtil")
    print("=" * 60)

    if not UTIL_AVAILABLE:
        print("\nutil.py not available, skipping this example\n")
        return

    try:
        # Load genomes
        loader = GenomeLoader()
        genome_ids = loader.list_available_genomes()[:3]

        print(f"\nProcessing {len(genome_ids)} genomes with CommUtil...\n")

        # Initialize CommUtil (your existing utility)
        # util = CommUtil()

        for genome_id in genome_ids:
            # Load genome
            genome = loader.load_genome(genome_id)

            # Create MSGenome
            msgenome = loader.get_msgenome(genome_id)

            print(f"✓ {genome_id} ready for modeling")
            # Now you can use msgenome with your existing workflow
            # Example: build models, run simulations, etc.

        print("\nGenomes integrated successfully!")
        print()

    except Exception as e:
        print(f"\nError: {e}\n")
        import traceback
        traceback.print_exc()


def example_4_create_genome_collection():
    """
    Example 4: Create a collection of genomes for community modeling
    """
    print("=" * 60)
    print("Example 4: Genome Collection for Community Modeling")
    print("=" * 60)

    loader = GenomeLoader()

    # Load all genomes
    print("\nLoading all genomes...")
    all_genomes = loader.load_all_genomes(verbose=False)

    # Filter by domain or taxonomy
    bacteria = {gid: g for gid, g in all_genomes.items()
                if g.get('domain') == 'Bacteria'}
    archaea = {gid: g for gid, g in all_genomes.items()
               if g.get('domain') == 'Archaea'}

    print(f"\nGenome collection:")
    print(f"  Total: {len(all_genomes)} genomes")
    print(f"  Bacteria: {len(bacteria)}")
    print(f"  Archaea: {len(archaea)}")

    # Group by taxonomy (e.g., by phylum)
    from collections import defaultdict
    by_phylum = defaultdict(list)

    for gid, genome in all_genomes.items():
        taxonomy = genome.get('taxonomy', '')
        parts = taxonomy.split(';')
        phylum = parts[1].strip() if len(parts) > 1 else 'Unknown'
        by_phylum[phylum].append(gid)

    print(f"\nGenomes by phylum:")
    for phylum, genome_list in sorted(by_phylum.items(),
                                     key=lambda x: len(x[1]),
                                     reverse=True)[:10]:
        print(f"  {phylum}: {len(genome_list)} genomes")

    print()


def example_5_community_model_workflow():
    """
    Example 5: Complete workflow for community metabolic modeling
    """
    print("=" * 60)
    print("Example 5: Community Modeling Workflow")
    print("=" * 60)

    try:
        from modelseedpy import MSBuilder
        from cobra.io import write_sbml_model

        loader = GenomeLoader()

        # Select genomes for community model
        # (e.g., based on 16S data or abundance)
        community_genomes = [
            'Acetobacterium.1',
            'Bacteroides.1',
            'Desulfovibrio.1'
        ]

        print(f"\nBuilding community with {len(community_genomes)} species...\n")

        # Step 1: Build individual models
        models = {}
        for genome_id in community_genomes:
            try:
                print(f"Building model for {genome_id}...")
                model = loader.build_model(genome_id)
                models[genome_id] = model

                # Save individual model
                os.makedirs('models', exist_ok=True)
                output_file = f"models/{genome_id}_model.xml"
                write_sbml_model(model, output_file)
                print(f"  Saved to {output_file}")

            except Exception as e:
                print(f"  Error: {e}")
                continue

        print(f"\n✓ Built {len(models)} individual models")

        # Step 2: Community model construction
        # (This would use MSCommunity or similar tools)
        print("\nReady for community model construction!")
        print("Next steps:")
        print("  1. Use MSCommunity to combine models")
        print("  2. Define shared metabolite pool")
        print("  3. Set abundance constraints")
        print("  4. Simulate community dynamics")
        print()

    except ImportError:
        print("\nModelSEEDpy not installed. Install with:")
        print("  pip install modelseedpy\n")


def example_6_export_for_analysis():
    """
    Example 6: Export genome data for external analysis
    """
    print("=" * 60)
    print("Example 6: Export Data for Analysis")
    print("=" * 60)

    loader = GenomeLoader()

    # Export summary table
    print("\nExporting genome summary...")
    loader.export_genome_summary('genome_summary.tsv')
    print("✓ Created genome_summary.tsv")

    # Export feature functions
    print("\nExporting feature functions...")
    all_genomes = loader.load_all_genomes(verbose=False)

    with open('genome_functions.tsv', 'w') as f:
        f.write("genome_id\tfeature_id\tfunction\n")
        for genome_id, genome in all_genomes.items():
            for feature in genome.get('features', []):
                for function in feature.get('functions', []):
                    f.write(f"{genome_id}\t{feature['id']}\t{function}\n")

    print("✓ Created genome_functions.tsv")

    # Export taxonomy
    print("\nExporting taxonomy...")
    with open('genome_taxonomy.tsv', 'w') as f:
        f.write("genome_id\tdomain\ttaxonomy\n")
        for genome_id, genome in all_genomes.items():
            domain = genome.get('domain', 'Unknown')
            taxonomy = genome.get('taxonomy', '')
            f.write(f"{genome_id}\t{domain}\t{taxonomy}\n")

    print("✓ Created genome_taxonomy.tsv")

    print("\nExport complete! Files ready for analysis in R, Python, or Excel")
    print()


def main():
    """Run all integration examples"""
    examples = [
        ("Basic Integration", example_1_basic_integration),
        ("Batch Model Building", example_2_build_models_batch),
        ("CommUtil Integration", example_3_integrate_with_util),
        ("Genome Collection", example_4_create_genome_collection),
        ("Community Modeling", example_5_community_model_workflow),
        ("Export for Analysis", example_6_export_for_analysis),
    ]

    if len(sys.argv) > 1:
        # Run specific example
        try:
            example_num = int(sys.argv[1])
            if 1 <= example_num <= len(examples):
                name, func = examples[example_num - 1]
                print(f"\nRunning: {name}\n")
                func()
            else:
                print(f"Example number must be between 1 and {len(examples)}")
        except ValueError:
            print("Usage: python integrate_with_existing_code.py [example_number]")
            print("\nAvailable examples:")
            for i, (name, _) in enumerate(examples, 1):
                print(f"  {i}. {name}")
    else:
        # Show menu
        print("\n" + "=" * 60)
        print("Integration Examples Menu")
        print("=" * 60)
        print("\nAvailable examples:")
        for i, (name, _) in enumerate(examples, 1):
            print(f"  {i}. {name}")
        print("\nUsage: python integrate_with_existing_code.py [example_number]")
        print("Example: python integrate_with_existing_code.py 1")
        print()


if __name__ == '__main__':
    main()
