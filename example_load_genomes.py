#!/usr/bin/env python3
"""
Example: Loading Genome Objects into ModelSEEDpy

This script demonstrates how to load KBase genome objects from the genome_objects
directory and use them with ModelSEEDpy for metabolic modeling.
"""

import sys
from load_genomes_to_modelseedpy import (
    list_available_genomes,
    load_genome,
    load_all_genomes,
    get_msgenome,
    build_model,
    export_genome_summary
)

def example_1_list_genomes():
    """Example 1: List all available genomes"""
    print("=" * 60)
    print("Example 1: List all available genomes")
    print("=" * 60)

    genomes = list_available_genomes()
    print(f"\nFound {len(genomes)} genomes in genome_objects/")
    print("\nFirst 10 genomes:")
    for genome_id in genomes[:10]:
        print(f"  - {genome_id}")
    print()


def example_2_load_single_genome():
    """Example 2: Load a single genome and inspect it"""
    print("=" * 60)
    print("Example 2: Load and inspect a single genome")
    print("=" * 60)

    # Load one genome
    genome_id = 'Acetobacterium.1'
    print(f"\nLoading genome: {genome_id}")

    genome = load_genome(genome_id)

    # Inspect the genome
    print("\nGenome Information:")
    print(f"  ID: {genome['id']}")
    print(f"  Scientific name: {genome['scientific_name']}")
    print(f"  Domain: {genome['domain']}")
    print(f"  Taxonomy: {genome['taxonomy'][:100]}...")
    print(f"  Number of features: {len(genome['features'])}")
    print(f"  Number of CDSs: {len(genome['cdss'])}")
    print(f"  DNA size: {genome['dna_size']:,} bp")
    print(f"  GC content: {genome['gc_content']:.2%}")
    print(f"  Number of contigs: {genome['num_contigs']}")

    # Inspect a feature
    if genome['features']:
        feature = genome['features'][0]
        print("\nFirst feature:")
        print(f"  ID: {feature['id']}")
        print(f"  Type: {feature['type']}")
        print(f"  Functions: {feature.get('functions', [])}")
        print(f"  DNA length: {feature['dna_sequence_length']} bp")
        if 'protein_translation' in feature:
            print(f"  Protein length: {feature['protein_translation_length']} aa")
    print()


def example_3_load_all_genomes():
    """Example 3: Load all genomes"""
    print("=" * 60)
    print("Example 3: Load all genomes")
    print("=" * 60)

    genomes = load_all_genomes()

    print(f"\nLoaded {len(genomes)} genomes")

    # Summary statistics
    total_features = sum(len(g['features']) for g in genomes.values())
    total_dna = sum(g['dna_size'] for g in genomes.values())

    print(f"Total features across all genomes: {total_features:,}")
    print(f"Total DNA: {total_dna:,} bp ({total_dna/1e9:.2f} Gbp)")
    print()


def example_4_export_summary():
    """Example 4: Export genome summary to TSV"""
    print("=" * 60)
    print("Example 4: Export genome summary")
    print("=" * 60)

    output_file = "genome_summary.tsv"
    print(f"\nExporting genome summary to {output_file}")

    export_genome_summary(output_file)

    print(f"\nSummary file created: {output_file}")
    print("You can open this in Excel, Google Sheets, or pandas")
    print()


def example_5_create_msgenome():
    """Example 5: Create MSGenome object for modeling"""
    print("=" * 60)
    print("Example 5: Create MSGenome object")
    print("=" * 60)

    try:
        genome_id = 'Acetobacterium.1'
        print(f"\nCreating MSGenome object for {genome_id}")

        msgenome = get_msgenome(genome_id)

        print(f"\nMSGenome created successfully!")
        print(f"  Genome ID: {msgenome.id}")
        print(f"  Number of features: {len(msgenome.features)}")
        print("\nThis MSGenome object can now be used for:")
        print("  - Building metabolic models")
        print("  - Gap-filling analysis")
        print("  - Metabolic pathway analysis")
        print()

    except ImportError:
        print("\nModelSEEDpy not installed. Install with:")
        print("  pip install modelseedpy")
        print()


def example_6_build_model():
    """Example 6: Build a metabolic model"""
    print("=" * 60)
    print("Example 6: Build a metabolic model")
    print("=" * 60)

    try:
        genome_id = 'Acetobacterium.1'
        print(f"\nBuilding metabolic model for {genome_id}")

        model = build_model(genome_id)

        print(f"\nModel statistics:")
        print(f"  Reactions: {len(model.reactions)}")
        print(f"  Metabolites: {len(model.metabolites)}")
        print(f"  Genes: {len(model.genes)}")

        # Try to optimize
        print("\nAttempting FBA optimization...")
        solution = model.optimize()

        if solution.status == 'optimal':
            print(f"  Status: {solution.status}")
            print(f"  Objective value: {solution.objective_value:.4f}")
        else:
            print(f"  Status: {solution.status}")
            print("  Model may need gap-filling")

        # Save model
        output_file = f"{genome_id}_model.xml"
        print(f"\nSaving model to {output_file}")

        from cobra.io import write_sbml_model
        write_sbml_model(model, output_file)
        print(f"  Model saved successfully!")
        print()

    except ImportError:
        print("\nModelSEEDpy not installed. Install with:")
        print("  pip install modelseedpy")
        print()
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        print()


def main():
    """Run all examples"""
    examples = [
        ("List genomes", example_1_list_genomes),
        ("Load single genome", example_2_load_single_genome),
        ("Load all genomes", example_3_load_all_genomes),
        ("Export summary", example_4_export_summary),
        ("Create MSGenome", example_5_create_msgenome),
        ("Build model", example_6_build_model),
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
            print("Usage: python example_load_genomes.py [example_number]")
            print("\nAvailable examples:")
            for i, (name, _) in enumerate(examples, 1):
                print(f"  {i}. {name}")
    else:
        # Run all examples
        print("\n" + "=" * 60)
        print("Running all examples")
        print("=" * 60 + "\n")

        for i, (name, func) in enumerate(examples, 1):
            try:
                func()
                input("Press Enter to continue to next example...")
                print("\n")
            except KeyboardInterrupt:
                print("\n\nExiting...")
                break
            except Exception as e:
                print(f"\nError in example {i}: {e}")
                import traceback
                traceback.print_exc()
                print()


if __name__ == '__main__':
    main()
