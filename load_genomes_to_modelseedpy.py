#!/usr/bin/env python3
"""
Load KBase Genome Objects into ModelSEEDpy

This script loads genome object JSON files from the genome_objects directory
and creates ModelSEEDpy MSGenome objects that can be used for metabolic modeling.

Usage:
    # Load all genomes from genome_objects directory
    python load_genomes_to_modelseedpy.py

    # Load specific genomes
    python load_genomes_to_modelseedpy.py Acetobacterium.1 Bacteroides.1

    # Load genomes and build models
    python load_genomes_to_modelseedpy.py --build-models

Examples:
    >>> from load_genomes_to_modelseedpy import load_genome, load_all_genomes
    >>>
    >>> # Load a single genome
    >>> genome = load_genome('Acetobacterium.1')
    >>>
    >>> # Load all genomes
    >>> genomes = load_all_genomes()
    >>>
    >>> # Get MSGenome object for modeling
    >>> ms_genome = get_msgenome('Acetobacterium.1')
"""

import os
import sys
import json
import argparse
from typing import Dict, List, Optional, Any
from glob import glob

# Try to import ModelSEEDpy components
try:
    from modelseedpy.core.msgenome import MSGenome
    MODELSEEDPY_AVAILABLE = True
except ImportError:
    print("Warning: ModelSEEDpy not installed. Install with:")
    print("  pip install modelseedpy")
    MODELSEEDPY_AVAILABLE = False


class GenomeLoader:
    """
    Loader for KBase genome objects to ModelSEEDpy format
    """

    def __init__(self, genome_dir: str = "genome_objects"):
        """
        Initialize genome loader

        Args:
            genome_dir: Directory containing genome JSON files
        """
        self.genome_dir = genome_dir
        self._genomes = {}  # Cache for loaded genomes
        self._msgenomes = {}  # Cache for MSGenome objects

    def list_available_genomes(self) -> List[str]:
        """
        List all available genome IDs in the genome_objects directory

        Returns:
            List of genome IDs (without .json extension)
        """
        pattern = os.path.join(self.genome_dir, "*.json")
        files = glob(pattern)
        return sorted([os.path.splitext(os.path.basename(f))[0] for f in files])

    def load_genome(self, genome_id: str, printing: bool = True) -> Dict[str, Any]:
        """
        Load a genome object from JSON file

        Args:
            genome_id: Genome identifier (e.g., 'Acetobacterium.1')
            printing: Print genome statistics (default: True)

        Returns:
            KBase Genome object dictionary
        """
        # Check cache first
        if genome_id in self._genomes:
            return self._genomes[genome_id]

        # Load from file
        genome_file = os.path.join(self.genome_dir, f"{genome_id}.json")

        if not os.path.exists(genome_file):
            raise FileNotFoundError(
                f"Genome file not found: {genome_file}\n"
                f"Available genomes: {', '.join(self.list_available_genomes()[:5])}..."
            )

        if printing:
            print(f"Loading genome: {genome_id}")
        with open(genome_file, 'r') as f:
            genome = json.load(f)

        # Cache the genome
        self._genomes[genome_id] = genome

        # Print summary
        if printing:
            print(f"  Scientific name: {genome.get('scientific_name', 'N/A')}")
            print(f"  Taxonomy: {genome.get('taxonomy', 'N/A')[:80]}...")
            print(f"  Features: {len(genome.get('features', []))}")
            print(f"  DNA size: {genome.get('dna_size', 0):,} bp")
            print(f"  GC content: {genome.get('gc_content', 0):.2%}")

        return genome

    def load_all_genomes(self, verbose: bool = True, printing: bool = True) -> Dict[str, Dict[str, Any]]:
        """
        Load all genomes from the genome_objects directory

        Args:
            verbose: Print progress information (default: True)
            printing: Print individual genome statistics (default: True)

        Returns:
            Dictionary mapping genome_id -> genome object
        """
        genome_ids = self.list_available_genomes()

        if verbose:
            print(f"Loading {len(genome_ids)} genomes from {self.genome_dir}")

        genomes = {}
        for genome_id in genome_ids:
            try:
                genome = self.load_genome(genome_id, printing=printing)
                genomes[genome_id] = genome
            except Exception as e:
                # Always print errors regardless of verbose/printing settings
                print(f"  Error loading {genome_id}: {e}")
                continue

        if verbose:
            print(f"\nSuccessfully loaded {len(genomes)} genomes")

        return genomes

    def get_msgenome(self, genome_id: str, printing: bool = True) -> 'MSGenome':
        """
        Get MSGenome object for a genome

        Args:
            genome_id: Genome identifier
            printing: Print status messages (default: True)

        Returns:
            MSGenome object ready for metabolic modeling

        Raises:
            ImportError: If ModelSEEDpy is not installed
        """
        if not MODELSEEDPY_AVAILABLE:
            raise ImportError(
                "ModelSEEDpy is not installed. Install with: pip install modelseedpy"
            )

        # Check cache first
        if genome_id in self._msgenomes:
            return self._msgenomes[genome_id]

        # Load genome if not already loaded
        genome = self.load_genome(genome_id, printing=printing)

        # Create MSGenome object
        if printing:
            print(f"Creating MSGenome object for {genome_id}")
            print(f"  Converting {len(genome.get('features', []))} features...")

        from modelseedpy.core.msgenome import MSGenome, MSFeature

        # Create empty MSGenome
        msgenome = MSGenome()

        # Set basic attributes (excluding features which need special handling)
        for key, value in genome.items():
            if key != 'features':  # Skip features, handle them specially
                try:
                    setattr(msgenome, key, value)
                except:
                    pass
        msgenome.id = genome_id

        # Convert feature dicts to MSFeature objects
        if 'features' in genome and genome['features']:
            msfeatures = []
            for feat_dict in genome['features']:
                try:
                    # Create MSFeature with required fields
                    feature_id = feat_dict.get('id', '')
                    sequence = feat_dict.get('dna_sequence', feat_dict.get('protein_translation', ''))
                    description = feat_dict.get('function', feat_dict.get('type', ''))
                    aliases = feat_dict.get('aliases', [])

                    msfeature = MSFeature(feature_id, sequence, description=description, aliases=aliases)

                    # Add ontology terms if present
                    if 'ontology_terms' in feat_dict:
                        for term_type, terms in feat_dict['ontology_terms'].items():
                            if isinstance(terms, list):
                                for term in terms:
                                    msfeature.add_ontology_term(term_type, term)
                            else:
                                msfeature.add_ontology_term(term_type, terms)

                    msfeatures.append(msfeature)
                except Exception as e:
                    # Skip features that can't be converted
                    pass

            # Add all features to MSGenome
            msgenome.features = msfeatures
            if printing:
                print(f"  ✓ Converted {len(msfeatures)} features to MSFeature objects")

        # Cache the MSGenome
        self._msgenomes[genome_id] = msgenome

        if printing:
            print(f"✓ MSGenome created: {msgenome.id}")

        return msgenome

    def build_model(
        self,
        genome_id: str,
        template: Optional[str] = None,
        model_id: Optional[str] = None,
        printing: bool = True
    ):
        """
        Build a metabolic model from a genome using ModelSEEDpy

        Args:
            genome_id: Genome identifier
            template: Template to use (default: GramNegative or GramPositive based on taxonomy)
            model_id: Model ID (default: genome_id with _model suffix)
            printing: Print status messages and model statistics (default: True)

        Returns:
            COBRApy model object

        Raises:
            ImportError: If ModelSEEDpy is not installed
        """
        if not MODELSEEDPY_AVAILABLE:
            raise ImportError(
                "ModelSEEDpy is not installed. Install with: pip install modelseedpy"
            )

        # Get MSGenome object
        msgenome = self.get_msgenome(genome_id, printing=printing)

        # Load genome to determine template
        genome = self.load_genome(genome_id, printing=False)  # Already loaded, don't print again

        # Load template - use template_core.json
        if printing:
            print(f"Loading template...")

        import json
        from pathlib import Path as PathlibPath
        from os import path
        from modelseedpy.core.mstemplate import MSTemplateBuilder

        # Use template_core.json from ModelSEEDpy repo
        template_path = path.expanduser('~/repos/ModelSEEDpy/modelseedpy/data/templates/template_core.json')

        try:
            with open(template_path, 'r') as f:
                template_data = json.load(f)
            template_obj = MSTemplateBuilder.from_dict(template_data).build()
            if printing:
                print(f"  ✓ Loaded template: {template_obj.id} ({len(template_obj.reactions)} reactions)")
        except Exception as e:
            error_msg = f"Failed to load template: {e}"
            print(f"✗ {error_msg}")
            raise RuntimeError(error_msg)

        # Build model
        model_id = model_id or f"{genome_id}_model"
        if printing:
            print(f"Building model: {model_id}")

        from modelseedpy.core.msbuilder import MSBuilder

        builder = MSBuilder(msgenome, template_obj, name=model_id)
        model = builder.build(model_id, allow_all_non_grp_reactions=True, annotate_with_rast=False)

        if printing:
            print(f"✓ Model built successfully!")
            print(f"  Reactions: {len(model.reactions)}")
            print(f"  Metabolites: {len(model.metabolites)}")
            print(f"  Genes: {len(model.genes)}")

        return model

    def export_genome_summary(
        self,
        output_file: str = "genome_summary.tsv",
        genome_ids: Optional[List[str]] = None,
        printing: bool = True
    ):
        from pandas import read_csv
        from io import StringIO
        """
        Export a summary table of all genomes

        Args:
            output_file: Path to output TSV file
            genome_ids: List of genome IDs to include (default: all)
            printing: Print individual genome loading messages (default: True)
        """
        if genome_ids is None:
            genome_ids = self.list_available_genomes()

        print(f"Exporting summary for {len(genome_ids)} genomes to {output_file}")

        with open(output_file, 'w') as f:
            # Header
            headers = ["genome_id","scientific_name","domain","taxonomy",
                      "num_features","num_cdss","dna_size","gc_content","num_contigs"]
            header = "\t".join(headers)
            f.write(header)

            # Data rows
            rows = []
            for i, genome_id in enumerate(genome_ids):
                try:
                    genome = self.load_genome(genome_id, printing=printing)
                    rows.append([
                        f"{genome_id}", f"{genome.get('scientific_name', 'N/A')}",
                        f"{genome.get('domain', 'N/A')}", f"{genome.get('taxonomy', 'N/A')}",
                        f"{len(genome.get('features', []))}", f"{len(genome.get('cdss', []))}",
                        f"{genome.get('dna_size', 0)}", f"{genome.get('gc_content', 0):.4f}",
                        f"{genome.get('num_contigs', 0)}"
                    ])
                    rows[-1] = "\t".join(rows[-1])
                    
                    # f.write(f"{genome_id}\t")
                    # f.write(f"{genome.get('scientific_name', 'N/A')}\t")
                    # f.write(f"{genome.get('domain', 'N/A')}\t")
                    # f.write(f"{genome.get('taxonomy', 'N/A')}\t")
                    # f.write(f"{len(genome.get('features', []))}\t")
                    # f.write(f"{len(genome.get('cdss', []))}\t")
                    # f.write(f"{genome.get('dna_size', 0)}\t")
                    # f.write(f"{genome.get('gc_content', 0):.4f}\t")
                    # f.write(f"{genome.get('num_contigs', 0)}\n")

                except Exception as e:
                    print(f"  Error processing {genome_id}: {e}")
                    continue
            total_string = "\n".join(rows)
            f.write(total_string)

        df = read_csv(StringIO(total_string), sep="\t", header=None, names=headers).set_index("genome_id")
        print(f"Summary exported to {output_file}")
        return df


# ============================================================================
# Convenience Functions
# ============================================================================

# Global loader instance
_loader = None

def get_loader(genome_dir: str = "genome_objects") -> GenomeLoader:
    """Get or create the global GenomeLoader instance"""
    global _loader
    if _loader is None or _loader.genome_dir != genome_dir:
        _loader = GenomeLoader(genome_dir)
    return _loader


def list_available_genomes(genome_dir: str = "genome_objects") -> List[str]:
    """
    List all available genome IDs

    Args:
        genome_dir: Directory containing genome JSON files

    Returns:
        List of genome IDs

    Example:
        >>> genomes = list_available_genomes()
        >>> print(f"Found {len(genomes)} genomes")
    """
    loader = get_loader(genome_dir)
    return loader.list_available_genomes()


def load_genome(genome_id: str, genome_dir: str = "genome_objects", printing: bool = True) -> Dict[str, Any]:
    """
    Load a genome object from JSON file

    Args:
        genome_id: Genome identifier
        genome_dir: Directory containing genome JSON files
        printing: Print genome statistics (default: True)

    Returns:
        KBase Genome object dictionary

    Example:
        >>> genome = load_genome('Acetobacterium.1')
        >>> print(genome['scientific_name'])

        >>> # Load without printing
        >>> genome = load_genome('Acetobacterium.1', printing=False)
    """
    loader = get_loader(genome_dir)
    return loader.load_genome(genome_id, printing=printing)


def load_all_genomes(
    genome_dir: str = "genome_objects",
    verbose: bool = True,
    printing: bool = True
) -> Dict[str, Dict[str, Any]]:
    """
    Load all genomes from directory

    Args:
        genome_dir: Directory containing genome JSON files
        verbose: Print progress information (default: True)
        printing: Print individual genome statistics (default: True)

    Returns:
        Dictionary mapping genome_id -> genome object

    Example:
        >>> genomes = load_all_genomes()
        >>> print(f"Loaded {len(genomes)} genomes")

        >>> # Load quietly (only show overall progress, not individual genomes)
        >>> genomes = load_all_genomes(printing=False)
    """
    loader = get_loader(genome_dir)
    return loader.load_all_genomes(verbose=verbose, printing=printing)


def get_msgenome(genome_id: str, genome_dir: str = "genome_objects", printing: bool = True) -> 'MSGenome':
    """
    Get MSGenome object for modeling

    Args:
        genome_id: Genome identifier
        genome_dir: Directory containing genome JSON files
        printing: Print status messages (default: True)

    Returns:
        MSGenome object

    Example:
        >>> msgenome = get_msgenome('Acetobacterium.1')
        >>> # Use msgenome for model building

        >>> # Create MSGenome without printing
        >>> msgenome = get_msgenome('Acetobacterium.1', printing=False)
    """
    loader = get_loader(genome_dir)
    return loader.get_msgenome(genome_id, printing=printing)


def build_model(
    genome_id: str,
    template: Optional[str] = None,
    model_id: Optional[str] = None,
    genome_dir: str = "genome_objects",
    printing: bool = True
):
    """
    Build a metabolic model from a genome

    Args:
        genome_id: Genome identifier
        template: Template to use (default: auto-select based on taxonomy)
        model_id: Model ID (default: genome_id with _model suffix)
        genome_dir: Directory containing genome JSON files
        printing: Print status messages and model statistics (default: True)

    Returns:
        COBRApy model object

    Example:
        >>> model = build_model('Acetobacterium.1')
        >>> solution = model.optimize()
        >>> print(f"Growth rate: {solution.objective_value}")

        >>> # Build model quietly
        >>> model = build_model('Acetobacterium.1', printing=False)
    """
    loader = get_loader(genome_dir)
    return loader.build_model(genome_id, template=template, model_id=model_id, printing=printing)


def export_genome_summary(
    output_file: str = "genome_summary.tsv",
    genome_ids: Optional[List[str]] = None,
    genome_dir: str = "genome_objects",
    printing: bool = True
):
    """
    Export genome summary table

    Args:
        output_file: Path to output TSV file
        genome_ids: List of genome IDs to include (default: all)
        genome_dir: Directory containing genome JSON files
        printing: Print individual genome loading messages (default: True)

    Example:
        >>> export_genome_summary('my_genomes.tsv')

        >>> # Export quietly
        >>> export_genome_summary('my_genomes.tsv', printing=False)
    """
    loader = get_loader(genome_dir)
    return loader.export_genome_summary(output_file=output_file, genome_ids=genome_ids, printing=printing)


# ============================================================================
# Command-Line Interface
# ============================================================================

def main():
    """Main entry point for command-line usage"""
    parser = argparse.ArgumentParser(
        description='Load KBase genome objects into ModelSEEDpy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('genome_ids', nargs='*',
                       help='Genome IDs to load (default: all genomes)')
    parser.add_argument('--genome-dir', default='genome_objects',
                       help='Directory containing genome JSON files (default: genome_objects)')
    parser.add_argument('--list', action='store_true',
                       help='List all available genomes')
    parser.add_argument('--summary', metavar='FILE',
                       help='Export genome summary to TSV file')
    parser.add_argument('--build-models', action='store_true',
                       help='Build metabolic models for the genomes')
    parser.add_argument('--template',
                       help='Template to use for model building (default: auto-select)')
    parser.add_argument('--save-models', metavar='DIR',
                       help='Directory to save SBML models')

    args = parser.parse_args()

    # Create loader
    loader = GenomeLoader(genome_dir=args.genome_dir)

    # List genomes
    if args.list:
        genomes = loader.list_available_genomes()
        print(f"Available genomes ({len(genomes)}):")
        for genome_id in genomes:
            print(f"  {genome_id}")
        return

    # Determine which genomes to process
    if args.genome_ids:
        genome_ids = args.genome_ids
    else:
        genome_ids = loader.list_available_genomes()

    print(f"Processing {len(genome_ids)} genomes")
    print("=" * 50)

    # Load genomes
    genomes = {}
    for genome_id in genome_ids:
        try:
            genome = loader.load_genome(genome_id)
            genomes[genome_id] = genome
            print()
        except Exception as e:
            print(f"Error loading {genome_id}: {e}\n")
            continue

    # Export summary if requested
    if args.summary:
        loader.export_genome_summary(
            output_file=args.summary,
            genome_ids=list(genomes.keys())
        )

    # Build models if requested
    if args.build_models:
        if not MODELSEEDPY_AVAILABLE:
            print("\nError: ModelSEEDpy not installed. Install with: pip install modelseedpy")
            return

        print("\n" + "=" * 50)
        print("Building metabolic models")
        print("=" * 50)

        models = {}
        for genome_id in genomes.keys():
            try:
                model = loader.build_model(
                    genome_id,
                    template=args.template
                )
                models[genome_id] = model

                # Save model if requested
                if args.save_models:
                    os.makedirs(args.save_models, exist_ok=True)
                    model_file = os.path.join(args.save_models, f"{genome_id}_model.xml")

                    from cobra.io import write_sbml_model
                    write_sbml_model(model, model_file)
                    print(f"  Saved to {model_file}")

                print()

            except Exception as e:
                print(f"Error building model for {genome_id}: {e}\n")
                import traceback
                traceback.print_exc()
                continue

        print(f"\nSuccessfully built {len(models)} models")

    print("\n" + "=" * 50)
    print("Done!")


if __name__ == '__main__':
    main()
