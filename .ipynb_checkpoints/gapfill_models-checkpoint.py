#!/usr/bin/env python3
"""
Gapfill SBML models using ModelSEEDpy with auxotrophic and anaerobic minimal media.

This script adapts KBase's ModelSEEDReconstruction gapfilling workflow for local use
with ModelSEEDpy. It processes all SBML models in the models/ directory.

Usage:
    # Gapfill all models
    python gapfill_models.py

    # Gapfill specific models
    python gapfill_models.py Acetobacterium.1 Bacteroides.1

    # Gapfill with custom output directory
    python gapfill_models.py --output-dir gapfilled_models/

    # Dry run (list models without processing)
    python gapfill_models.py --dry-run
"""

import os
import sys
import argparse
import logging
from glob import glob
from typing import List, Optional, Dict, Any

# COBRApy for model I/O
from cobra.io import read_sbml_model, write_sbml_model

# ModelSEEDpy for gapfilling
from modelseedpy import MSMedia, MSGapfill, MSModelUtil
from modelseedpy.helpers import get_template
from modelseedpy.core.mstemplate import MSTemplateBuilder

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_auxo_media() -> MSMedia:
    """
    Create auxotrophic media (AuxoMedia).

    This is a minimal media supplemented with amino acids and vitamins
    to support auxotrophic organisms that cannot synthesize certain compounds.

    Based on KBase's AuxoMedia definition.
    """
    # Core minimal media compounds with ModelSEED IDs
    # Format: cpd_id -> uptake rate (positive value)
    compounds = {
        # Inorganic nutrients
        "cpd00009": 1000,   # Phosphate
        "cpd00013": 1000,   # NH3 (ammonia)
        "cpd00048": 1000,   # Sulfate
        "cpd00034": 1000,   # Zn2+
        "cpd00058": 1000,   # Cu2+
        "cpd00063": 1000,   # Ca2+
        "cpd00099": 1000,   # Cl-
        "cpd00149": 1000,   # Co2+
        "cpd00030": 1000,   # Mn2+
        "cpd00254": 1000,   # Mg
        "cpd00971": 1000,   # Na+
        "cpd00205": 1000,   # K+
        "cpd00001": 1000,   # H2O
        "cpd00067": 1000,   # H+
        "cpd10515": 1000,   # Fe2+
        "cpd10516": 1000,   # Fe3+

        # Amino acids (for auxotrophs)
        "cpd00035": 10,    # L-Alanine
        "cpd00051": 10,    # L-Arginine
        "cpd00132": 10,    # L-Asparagine
        "cpd00041": 10,    # L-Aspartate
        "cpd00084": 10,    # L-Cysteine
        "cpd00053": 10,    # L-Glutamine
        "cpd00023": 10,    # L-Glutamate
        "cpd00033": 10,    # Glycine
        "cpd00119": 10,    # L-Histidine
        "cpd00322": 10,    # L-Isoleucine
        "cpd00107": 10,    # L-Leucine
        "cpd00039": 10,    # L-Lysine
        "cpd00060": 10,    # L-Methionine
        "cpd00066": 10,    # L-Phenylalanine
        "cpd00129": 10,    # L-Proline
        "cpd00054": 10,    # L-Serine
        "cpd00161": 10,    # L-Threonine
        "cpd00065": 10,    # L-Tryptophan
        "cpd00069": 10,    # L-Tyrosine
        "cpd00156": 10,    # L-Valine

        # Vitamins and cofactors
        "cpd00104": 10,    # Biotin
        "cpd00393": 10,    # Folate
        "cpd00263": 10,    # Pantothenate
        "cpd00218": 10,    # Niacin
        "cpd00220": 10,    # Riboflavin
        "cpd00305": 10,    # Thiamine
        "cpd00215": 10,    # Pyridoxine
        "cpd00166": 10,    # Cobalamin (B12)

        # Carbon source - glucose
        "cpd00027": 10,    # D-Glucose
    }

    # MSMedia.from_dict expects: {'cpd_id': uptake_rate} or {'cpd_id': (lower, upper)}
    media = MSMedia.from_dict(compounds)
    media.id = "AuxoMedia"
    media.name = "Auxotrophic Minimal Media"

    return media


def create_pyruvate_minimal_anaerobic() -> MSMedia:
    """
    Create pyruvate minimal anaerobic media (PyruateMinimalAnaerobic).

    This is a minimal anaerobic media with pyruvate as the carbon source,
    suitable for fermentative and anaerobic organisms.

    Based on KBase's PyruateMinimalAnaerobic definition.
    """
    # Anaerobic minimal media with pyruvate
    # Format: cpd_id -> uptake rate (positive value)
    compounds = {
        # Inorganic nutrients
        "cpd00009": 1000,   # Phosphate
        "cpd00013": 1000,   # NH3 (ammonia)
        "cpd00048": 1000,   # Sulfate
        "cpd00034": 1000,   # Zn2+
        "cpd00058": 1000,   # Cu2+
        "cpd00063": 1000,   # Ca2+
        "cpd00099": 1000,   # Cl-
        "cpd00149": 1000,   # Co2+
        "cpd00030": 1000,   # Mn2+
        "cpd00254": 1000,   # Mg
        "cpd00971": 1000,   # Na+
        "cpd00205": 1000,   # K+
        "cpd00001": 1000,   # H2O
        "cpd00067": 1000,   # H+
        "cpd10515": 1000,   # Fe2+
        "cpd10516": 1000,   # Fe3+

        # Carbon source - pyruvate (anaerobic-friendly)
        "cpd00020": 10,    # Pyruvate

        # CO2 for carboxylation reactions
        "cpd00011": 1000,   # CO2

        # Bicarbonate buffer
        "cpd00242": 1000,   # HCO3-

        # H2 for hydrogenotrophic reactions
        "cpd11640": 1000,  # H2

        # Acetate (common anaerobic metabolite/substrate)
        "cpd00029": 10,    # Acetate
    }

    # MSMedia.from_dict expects: {'cpd_id': uptake_rate} or {'cpd_id': (lower, upper)}
    media = MSMedia.from_dict(compounds)
    media.id = "PyruateMinimalAnaerobic"
    media.name = "Pyruvate Minimal Anaerobic Media"

    return media


def get_template_for_model(model) -> Any:
    """
    Determine and load appropriate template based on model characteristics.

    Args:
        model: COBRApy model object

    Returns:
        MSTemplate object
    """
    # Check model notes/annotations for taxonomy hints
    taxonomy = ""
    if hasattr(model, 'notes') and model.notes:
        taxonomy = str(model.notes).lower()

    # Gram-positive indicators
    gram_positive_indicators = [
        'firmicutes', 'actinobacteria', 'bacill', 'clostrid',
        'lactobacill', 'streptococc', 'staphylococc', 'enterococc'
    ]

    # Check model name for hints
    model_name = model.id.lower() if model.id else ""

    if any(indicator in taxonomy or indicator in model_name
           for indicator in gram_positive_indicators):
        template_name = 'template_gram_pos'
    else:
        template_name = 'template_gram_neg'

    template_data = get_template(template_name)
    template = MSTemplateBuilder.from_dict(template_data).build()

    return template


def gapfill_model(
    model_path: str,
    media_list: List[MSMedia],
    output_dir: str,
    gapfilling_mode: str = "Sequential",
    atp_safe: bool = True
) -> Optional[str]:
    """
    Gapfill a single model with the specified media.

    Args:
        model_path: Path to SBML model file
        media_list: List of MSMedia objects for gapfilling
        output_dir: Directory to save gapfilled model
        gapfilling_mode: "Sequential" or "Simultaneous"
        atp_safe: Whether to ensure ATP production safety

    Returns:
        Path to gapfilled model if successful, None otherwise
    """
    model_name = os.path.splitext(os.path.basename(model_path))[0]

    try:
        # Load model
        logger.info(f"Loading model: {model_name}")
        model = read_sbml_model(model_path)

        # Get appropriate template
        logger.info(f"  Loading template...")
        template = get_template_for_model(model)
        logger.info(f"  Using template: {template.id}")

        # Create model utility wrapper
        model_util = MSModelUtil.get(model)

        # Track gapfilled reactions
        total_gapfilled = []

        # Perform gapfilling with each media
        for media in media_list:
            logger.info(f"  Gapfilling with {media.id}...")

            try:
                # Create gapfiller
                gapfiller = MSGapfill(
                    model,
                    default_gapfill_templates=[template],
                    default_gapfill_models=[],
                    test_conditions=[],
                    reaction_scores={},
                    blacklist=[]
                )

                # Run gapfilling
                gapfill_solution = gapfiller.run_gapfilling(
                    media=media,
                    target="bio1",  # Standard biomass reaction ID
                    minimum_obj=0.01,
                    binary_check=False,
                    prefilter=True
                )

                if gapfill_solution:
                    # Integrate gapfill solution
                    # integrate_gapfill_solution is an instance method - model is self.model
                    integrated = gapfiller.integrate_gapfill_solution(
                        gapfill_solution,
                        cumulative_solution=total_gapfilled,
                        gapfilling_mode=gapfilling_mode
                    )
                    if isinstance(integrated, dict):
                        new_rxns = integrated.get('new', [])
                    else:
                        new_rxns = integrated if integrated else []
                    total_gapfilled.extend(new_rxns)
                    logger.info(f"    Added {len(new_rxns)} reactions")
                else:
                    logger.warning(f"    No gapfill solution found for {media.id}")

            except Exception as e:
                logger.warning(f"    Gapfilling with {media.id} failed: {e}")
                continue

        # ATP safety check if requested
        if atp_safe and total_gapfilled:
            logger.info(f"  Running ATP safety check...")
            try:
                from modelseedpy import MSATPCorrection
                atp_checker = MSATPCorrection(model, template)
                atp_checker.run_atp_correction()
            except Exception as e:
                logger.warning(f"    ATP correction failed: {e}")

        # Save gapfilled model
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{model_name}_gf.xml")
        write_sbml_model(model, output_path)

        logger.info(f"  Saved gapfilled model: {output_path}")
        logger.info(f"  Total reactions added: {len(total_gapfilled)}")

        return output_path

    except Exception as e:
        logger.error(f"Error processing {model_name}: {e}")
        import traceback
        traceback.print_exc()
        return None


def list_models(models_dir: str = "models") -> List[str]:
    """
    List all SBML model files in the models directory.

    Args:
        models_dir: Directory containing model files

    Returns:
        List of model file paths
    """
    pattern = os.path.join(models_dir, "*.xml")
    models = sorted(glob(pattern))

    # Filter out community models (those ending with _comm.xml)
    models = [m for m in models if not m.endswith("_comm.xml")]

    return models


def gapfill_all_models(
    models_dir: str = "models",
    output_dir: str = "gapfilled_models",
    model_ids: Optional[List[str]] = None,
    gapfilling_mode: str = "Sequential",
    atp_safe: bool = True,
    dry_run: bool = False
) -> Dict[str, Optional[str]]:
    """
    Gapfill all models in the models directory.

    Args:
        models_dir: Directory containing SBML models
        output_dir: Directory to save gapfilled models
        model_ids: Optional list of specific model IDs to process
        gapfilling_mode: "Sequential" or "Simultaneous"
        atp_safe: Whether to ensure ATP production safety
        dry_run: If True, list models without processing

    Returns:
        Dictionary mapping model names to output paths (or None if failed)
    """
    # Create media
    logger.info("Creating media definitions...")
    auxo_media = create_auxo_media()
    anaerobic_media = create_pyruvate_minimal_anaerobic()
    media_list = [auxo_media, anaerobic_media]
    logger.info(f"  Created {len(media_list)} media: {[m.id for m in media_list]}")

    # Get list of models
    all_models = list_models(models_dir)
    logger.info(f"Found {len(all_models)} models in {models_dir}")

    # Filter to specific models if requested
    if model_ids:
        model_paths = []
        for model_id in model_ids:
            # Try with and without .xml extension
            path = os.path.join(models_dir, f"{model_id}.xml")
            if os.path.exists(path):
                model_paths.append(path)
            elif os.path.exists(os.path.join(models_dir, model_id)):
                model_paths.append(os.path.join(models_dir, model_id))
            else:
                logger.warning(f"Model not found: {model_id}")
        all_models = model_paths

    if dry_run:
        logger.info("Dry run - models that would be processed:")
        for model_path in all_models:
            print(f"  {os.path.basename(model_path)}")
        return {}

    # Process each model
    results = {}
    total = len(all_models)

    for i, model_path in enumerate(all_models, 1):
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        logger.info(f"\n[{i}/{total}] Processing: {model_name}")

        output_path = gapfill_model(
            model_path=model_path,
            media_list=media_list,
            output_dir=output_dir,
            gapfilling_mode=gapfilling_mode,
            atp_safe=atp_safe
        )

        results[model_name] = output_path

    # Summary
    successful = sum(1 for v in results.values() if v is not None)
    logger.info(f"\n{'='*50}")
    logger.info(f"Gapfilling complete!")
    logger.info(f"  Successful: {successful}/{total}")
    logger.info(f"  Failed: {total - successful}/{total}")
    logger.info(f"  Output directory: {output_dir}")

    return results


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Gapfill SBML models using ModelSEEDpy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'model_ids', nargs='*',
        help='Specific model IDs to process (default: all models)'
    )
    parser.add_argument(
        '--models-dir', default='models',
        help='Directory containing SBML models (default: models)'
    )
    parser.add_argument(
        '--output-dir', default='gapfilled_models',
        help='Directory to save gapfilled models (default: gapfilled_models)'
    )
    parser.add_argument(
        '--gapfilling-mode', default='Sequential',
        choices=['Sequential', 'Simultaneous'],
        help='Gapfilling mode (default: Sequential)'
    )
    parser.add_argument(
        '--no-atp-safe', action='store_true',
        help='Disable ATP safety check'
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='List models without processing'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enable verbose output'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    results = gapfill_all_models(
        models_dir=args.models_dir,
        output_dir=args.output_dir,
        model_ids=args.model_ids if args.model_ids else None,
        gapfilling_mode=args.gapfilling_mode,
        atp_safe=not args.no_atp_safe,
        dry_run=args.dry_run
    )

    # Return exit code based on success
    if results:
        failed = sum(1 for v in results.values() if v is None)
        sys.exit(1 if failed > 0 else 0)


if __name__ == '__main__':
    main()
