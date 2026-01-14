#!/usr/bin/env python3
"""
Gapfill SBML models using ModelSEEDpy with auxotrophic and anaerobic minimal media.

This script adapts KBase's ModelSEEDReconstruction gapfilling workflow for local use
with ModelSEEDpy. It processes all SBML models in the models/ directory.

Genome Classification:
    The script classifies each genome in the genomes directory to select the appropriate
    gapfilling template. Classification uses:
    1. MSPredict classifier (when annotated genomes with RAST terms are available)
    2. GTDB taxonomy lineage from genome_hash.json
    3. Model-based heuristics as fallback

    Classifications (string labels for JSON export):
    - "Gram-P": Bacillota/Firmicutes, Actinobacteriota
    - "Gram-N": Proteobacteria, Bacteroidota, etc.
    - "Cyanobacteria": Cyanobacteriota
    - "Archaea": All archaeal phyla

Usage:
    # Gapfill all models with genome classification
    python gapfill_models.py

    # Gapfill specific models
    python gapfill_models.py Acetobacterium.1 Bacteroides.1

    # Gapfill with custom output directory
    python gapfill_models.py --output-dir gapfilled_models/

    # Dry run (list models with classifications)
    python gapfill_models.py --dry-run

    # Disable genome classification (use model heuristics only)
    python gapfill_models.py --no-genome-classification

    # Custom genome metadata location
    python gapfill_models.py --genome-hash path/to/genome_hash.json
"""

import os
import sys
import argparse
import logging
import json
from glob import glob
from typing import List, Optional, Dict, Any, Tuple, Union
from multiprocessing import Pool, cpu_count

# COBRApy for model I/O
from cobra.io import read_sbml_model, write_sbml_model

# ModelSEEDpy for gapfilling
from modelseedpy import MSMedia, MSGapfill, MSModelUtil
from modelseedpy.helpers import get_template, get_classifier
from modelseedpy.core.mstemplate import MSTemplateBuilder, MSTemplate
from modelseedpy.core.msgenome import MSGenome

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Genome classification string constants
GENOME_CLASS_GRAM_P = "Gram Positive"
GENOME_CLASS_GRAM_N = "Gram Negative"
GENOME_CLASS_ARCHAEA = "Archaea"
GENOME_CLASS_CYANOBACTERIA = "Cyano"

# Valid genome classifications for validation
VALID_GENOME_CLASSES = {GENOME_CLASS_GRAM_P, GENOME_CLASS_GRAM_N, GENOME_CLASS_ARCHAEA, GENOME_CLASS_CYANOBACTERIA}


# Mapping from GTDB phylum names to genome classification
PHYLUM_TO_CLASS = {
    # Gram-positive phyla
    "bacillota": GENOME_CLASS_GRAM_P,        # Firmicutes
    "firmicutes": GENOME_CLASS_GRAM_P,
    "actinobacteriota": GENOME_CLASS_GRAM_P,  # Actinobacteria
    "actinobacteria": GENOME_CLASS_GRAM_P,
    "coriobacteriota": GENOME_CLASS_GRAM_P,
    "chloroflexi": GENOME_CLASS_GRAM_P,       # Some are gram-variable but closer to pos
    # Cyanobacteria
    "cyanobacteria": GENOME_CLASS_CYANOBACTERIA,
    "cyanobacteriota": GENOME_CLASS_CYANOBACTERIA,
    # Archaea domain
    "euryarchaeota": GENOME_CLASS_ARCHAEA,
    "crenarchaeota": GENOME_CLASS_ARCHAEA,
    "thaumarchaeota": GENOME_CLASS_ARCHAEA,
    "nanoarchaeota": GENOME_CLASS_ARCHAEA,
    "korarchaeota": GENOME_CLASS_ARCHAEA,
    "asgardarchaeota": GENOME_CLASS_ARCHAEA,
    "halobacteriota": GENOME_CLASS_ARCHAEA,
    "thermoplasmatota": GENOME_CLASS_ARCHAEA,
    "methanobacteriota": GENOME_CLASS_ARCHAEA,
    # Everything else defaults to Gram-negative
}


def classify_genome_from_taxonomy(
    taxonomy: str,
    domain: Optional[str] = None
) -> str:
    """
    Classify a genome based on GTDB taxonomy lineage.

    Args:
        taxonomy: GTDB lineage string (e.g., "d__Bacteria;p__Bacillota;c__Bacilli...")
        domain: Optional domain code ('A' for Archaea, 'B' for Bacteria)

    Returns:
        Genome classification string: "Gram-P", "Gram-N", "Archaea", or "Cyanobacteria"
    """
    taxonomy_lower = taxonomy.lower()

    # Check domain first
    if domain == "A" or "d__archaea" in taxonomy_lower:
        return GENOME_CLASS_ARCHAEA

    # Parse phylum from GTDB lineage format
    phylum = None
    for part in taxonomy.split(";"):
        if part.startswith("p__"):
            phylum = part[3:].lower().strip()
            break

    if phylum:
        # Check against known phylum mappings
        for phylum_key, genome_class in PHYLUM_TO_CLASS.items():
            if phylum_key in phylum:
                return genome_class

    # Default to Gram-negative for bacteria/unknown
    return GENOME_CLASS_GRAM_N


def classify_genome(
    genome: Optional[MSGenome] = None,
    taxonomy: Optional[str] = None,
    domain: Optional[str] = None,
    classifier_name: str = "knn_ACNP_RAST_filter_01_17_2023"
) -> str:
    """
    Classify a genome using MSPredict classifier or taxonomy data.

    This function uses a hybrid approach:
    1. If an MSGenome with RAST annotations is provided, use the MSPredict classifier
    2. Otherwise, fall back to taxonomy-based classification

    Args:
        genome: Optional MSGenome object with RAST ontology terms
        taxonomy: Optional GTDB taxonomy lineage string
        domain: Optional domain code ('A' for Archaea, 'B' for Bacteria)
        classifier_name: Name of the classifier to use for MSPredict

    Returns:
        Genome classification string: "Gram Positive", "Gram Negative", "Archaea", or "Cyano"
    """
    # Try MSPredict classifier if genome with annotations is provided
    if genome is not None:
        # Check if genome has RAST annotations
        has_rast = any(
            "RAST" in feature.ontology_terms
            for feature in genome.features
            if hasattr(feature, 'ontology_terms')
        )

        if has_rast:
            try:
                logger.debug(f"Classifying genome using MSPredict classifier")
                genome_classifier = get_classifier(classifier_name)
                genome_class_str = genome_classifier.classify(genome)

                # Map classifier result to string constants
                class_map = {
                    "P": GENOME_CLASS_GRAM_P,
                    "N": GENOME_CLASS_GRAM_N,
                    "C": GENOME_CLASS_CYANOBACTERIA,
                    "A": GENOME_CLASS_ARCHAEA
                }
                if genome_class_str in class_map:
                    return class_map[genome_class_str]
            except Exception as e:
                logger.warning(f"MSPredict classification failed: {e}, falling back to taxonomy")

    # Fall back to taxonomy-based classification
    if taxonomy:
        return classify_genome_from_taxonomy(taxonomy, domain)

    # Default to Gram-negative if no classification data available
    logger.warning("No classification data available, defaulting to Gram-negative")
    return GENOME_CLASS_GRAM_N


def get_template_for_genome_class(genome_class: str) -> Tuple[MSTemplate, MSTemplate]:
    """
    Get the appropriate templates for a genome class.

    Args:
        genome_class: Genome classification string ("Gram-P", "Gram-N", "Archaea", or "Cyanobacteria")

    Returns:
        Tuple of (core_template, genome_scale_template)
    """
    if genome_class not in VALID_GENOME_CLASSES:
        raise ValueError(f"Invalid genome class: {genome_class}. Must be one of {VALID_GENOME_CLASSES}")

    template_genome_scale_map = {
        GENOME_CLASS_GRAM_N: "template_gram_neg",
        GENOME_CLASS_CYANOBACTERIA: "template_gram_neg",  # Cyano uses gram-neg base
        GENOME_CLASS_ARCHAEA: "template_gram_neg",  # Archaea uses gram-neg base
        GENOME_CLASS_GRAM_P: "template_gram_pos",
    }
    template_core_map = {
        GENOME_CLASS_ARCHAEA: "template_core",
        GENOME_CLASS_CYANOBACTERIA: "template_core",
        GENOME_CLASS_GRAM_N: "template_core",
        GENOME_CLASS_GRAM_P: "template_core",
    }

    d_template_core = get_template(template_core_map[genome_class])
    d_template_genome_scale = get_template(template_genome_scale_map[genome_class])

    template_core = MSTemplateBuilder.from_dict(d_template_core).build()
    template_genome_scale = MSTemplateBuilder.from_dict(d_template_genome_scale).build()

    return template_core, template_genome_scale


def load_genome_metadata(
    genome_hash_path: str = "Sludge/datacache/genome_hash.json"
) -> Dict[str, Dict[str, Any]]:
    """
    Load genome metadata from genome_hash.json.

    Args:
        genome_hash_path: Path to genome_hash.json file

    Returns:
        Dictionary mapping genome IDs to their metadata
    """
    if not os.path.exists(genome_hash_path):
        logger.warning(f"Genome hash file not found: {genome_hash_path}")
        return {}

    with open(genome_hash_path, 'r') as f:
        genome_hash = json.load(f)

    # Parse genome hash into ID -> metadata mapping
    metadata = {}
    for key, value in genome_hash.items():
        # Extract genome ID from key (format: "1352.8344__sludge_genomes.RAST")
        genome_id = key.split("__")[0] if "__" in key else key

        # Value is a list with metadata dict at index 10
        if isinstance(value, list) and len(value) > 10 and isinstance(value[10], dict):
            meta = value[10]
            metadata[genome_id] = {
                "taxonomy": meta.get("GTDB_lineage", meta.get("Taxonomy", "")),
                "domain": meta.get("Domain", ""),
                "name": meta.get("Name", ""),
                "source_id": meta.get("Source ID", "")
            }

    return metadata


def classify_genomes_in_directory(
    genomes_dir: str = "genomes",
    genome_hash_path: Optional[str] = "Sludge/datacache/genome_hash.json"
) -> Dict[str, str]:
    """
    Classify all genomes in a directory based on available metadata.

    Args:
        genomes_dir: Directory containing genome files (.fna)
        genome_hash_path: Path to genome_hash.json with taxonomy data

    Returns:
        Dictionary mapping genome IDs to their classification string
    """
    classifications = {}

    # Load genome metadata if available
    metadata = {}
    if genome_hash_path:
        metadata = load_genome_metadata(genome_hash_path)
        logger.info(f"Loaded metadata for {len(metadata)} genomes")

    # Get all genome files
    genome_files = glob(os.path.join(genomes_dir, "*.fna"))
    logger.info(f"Found {len(genome_files)} genome files in {genomes_dir}")

    for genome_path in genome_files:
        # Extract genome ID from filename (e.g., "100174.3.fna" -> "100174.3")
        genome_id = os.path.splitext(os.path.basename(genome_path))[0]

        # Try to classify using metadata
        if genome_id in metadata:
            meta = metadata[genome_id]
            genome_class = classify_genome_from_taxonomy(
                meta.get("taxonomy", ""),
                meta.get("domain", "")
            )
            logger.debug(f"Classified {genome_id} as {genome_class} from taxonomy")
        else:
            # Default classification
            genome_class = GENOME_CLASS_GRAM_N
            logger.debug(f"No metadata for {genome_id}, defaulting to {genome_class}")

        classifications[genome_id] = genome_class

    # Summary
    class_counts = {}
    for gc in classifications.values():
        class_counts[gc] = class_counts.get(gc, 0) + 1
    logger.info(f"Genome classification summary: {class_counts}")

    return classifications


def get_genome_id_for_model(model_id: str, genome_model_mapping: Optional[Dict[str, str]] = None) -> Optional[str]:
    """
    Get the genome ID associated with a model.

    Args:
        model_id: Model identifier (e.g., "Acetobacterium.1")
        genome_model_mapping: Optional explicit mapping of model IDs to genome IDs

    Returns:
        Genome ID if found, None otherwise
    """
    if genome_model_mapping and model_id in genome_model_mapping:
        return genome_model_mapping[model_id]
    return None


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


def get_template_for_model(
    model,
    genome_class: Optional[str] = None,
    genome_classifications: Optional[Dict[str, str]] = None,
    genome_model_mapping: Optional[Dict[str, str]] = None
) -> MSTemplate:
    """
    Determine and load appropriate template based on genome classification or model characteristics.

    This function uses a priority-based approach:
    1. If genome_class is provided directly, use it
    2. If genome_classifications is provided, look up the genome for this model
    3. Fall back to model-based heuristics (taxonomy hints in model name/notes)

    Args:
        model: COBRApy model object
        genome_class: Optional pre-computed genome classification string
        genome_classifications: Optional dict of genome_id -> classification string
        genome_model_mapping: Optional dict of model_id -> genome_id

    Returns:
        MSTemplate object
    """
    determined_class = None

    # Priority 1: Direct genome class provided
    if genome_class is not None:
        determined_class = genome_class
        logger.debug(f"Using provided genome class: {genome_class}")

    # Priority 2: Look up from genome classifications
    elif genome_classifications is not None:
        model_id = model.id if model.id else ""
        genome_id = get_genome_id_for_model(model_id, genome_model_mapping)
        if genome_id and genome_id in genome_classifications:
            determined_class = genome_classifications[genome_id]
            logger.debug(f"Found genome class for {model_id}: {determined_class}")

    # Priority 3: Model-based heuristics
    if determined_class is None:
        # Check model notes/annotations for taxonomy hints
        taxonomy = ""
        if hasattr(model, 'notes') and model.notes:
            taxonomy = str(model.notes).lower()

        # Gram-positive indicators
        gram_positive_indicators = [
            'firmicutes', 'actinobacteria', 'bacill', 'clostrid',
            'lactobacill', 'streptococc', 'staphylococc', 'enterococc',
            'bacillota', 'actinobacteriota'
        ]

        # Check model name for hints
        model_name = model.id.lower() if model.id else ""

        if any(indicator in taxonomy or indicator in model_name
               for indicator in gram_positive_indicators):
            determined_class = GENOME_CLASS_GRAM_P
        else:
            determined_class = GENOME_CLASS_GRAM_N

        logger.debug(f"Classified {model.id} from model heuristics: {determined_class}")

    # Get appropriate template
    template_name = 'template_gram_pos' if determined_class == GENOME_CLASS_GRAM_P else 'template_gram_neg'
    template_data = get_template(template_name)
    template = MSTemplateBuilder.from_dict(template_data).build()

    return template


def gapfill_model(
    model_path: str,
    media_list: List[MSMedia],
    output_dir: str,
    gapfilling_mode: str = "Sequential",
    atp_safe: bool = True,
    genome_class: Optional[str] = None,
    genome_classifications: Optional[Dict[str, str]] = None,
    genome_model_mapping: Optional[Dict[str, str]] = None
) -> Optional[str]:
    """
    Gapfill a single model with the specified media.

    Args:
        model_path: Path to SBML model file
        media_list: List of MSMedia objects for gapfilling
        output_dir: Directory to save gapfilled model
        gapfilling_mode: "Sequential" or "Simultaneous"
        atp_safe: Whether to ensure ATP production safety
        genome_class: Optional pre-computed genome classification string
        genome_classifications: Optional dict of genome_id -> classification string
        genome_model_mapping: Optional dict of model_id -> genome_id

    Returns:
        Path to gapfilled model if successful, None otherwise
    """
    model_name = os.path.splitext(os.path.basename(model_path))[0]

    try:
        # Load model
        logger.info(f"Loading model: {model_name}")
        model = read_sbml_model(model_path)

        # Get appropriate template based on genome classification
        logger.info(f"  Loading template...")
        template = get_template_for_model(
            model,
            genome_class=genome_class,
            genome_classifications=genome_classifications,
            genome_model_mapping=genome_model_mapping
        )
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


def _gapfill_worker(args: Tuple) -> Tuple[str, Optional[str]]:
    """
    Worker function for parallel gapfilling.

    Args:
        args: Tuple of (model_path, output_dir, gapfilling_mode, atp_safe, genome_class)
              genome_class is the string classification (e.g., "Gram-P", "Gram-N")

    Returns:
        Tuple of (model_name, output_path or None if failed)
    """
    model_path, output_dir, gapfilling_mode, atp_safe, genome_class = args
    model_name = os.path.splitext(os.path.basename(model_path))[0]

    # Create media in each worker (not picklable across processes)
    auxo_media = create_auxo_media()
    anaerobic_media = create_pyruvate_minimal_anaerobic()
    media_list = [auxo_media, anaerobic_media]

    try:
        output_path = gapfill_model(
            model_path=model_path,
            media_list=media_list,
            output_dir=output_dir,
            gapfilling_mode=gapfilling_mode,
            atp_safe=atp_safe,
            genome_class=genome_class
        )
        return (model_name, output_path)
    except Exception as e:
        logger.error(f"Worker error for {model_name}: {e}")
        return (model_name, None)


def gapfill_all_models(
    models_dir: str = "models",
    output_dir: str = "gapfilled_models",
    model_ids: Optional[List[str]] = None,
    gapfilling_mode: str = "Sequential",
    atp_safe: bool = True,
    dry_run: bool = False,
    n_workers: Optional[int] = None,
    genomes_dir: Optional[str] = "genomes",
    genome_hash_path: Optional[str] = "Sludge/datacache/genome_hash.json",
    genome_model_mapping: Optional[Dict[str, str]] = None
) -> Dict[str, Optional[str]]:
    """
    Gapfill all models in the models directory using parallel processing.

    This function classifies genomes before gapfilling to select the appropriate
    template for each model.

    Args:
        models_dir: Directory containing SBML models
        output_dir: Directory to save gapfilled models
        model_ids: Optional list of specific model IDs to process
        gapfilling_mode: "Sequential" or "Simultaneous"
        atp_safe: Whether to ensure ATP production safety
        dry_run: If True, list models without processing
        n_workers: Number of parallel workers (default: 1/4 of CPU cores)
        genomes_dir: Directory containing genome files for classification
        genome_hash_path: Path to genome_hash.json with taxonomy data
        genome_model_mapping: Optional dict mapping model_id -> genome_id

    Returns:
        Dictionary mapping model names to output paths (or None if failed)
    """
    # Determine number of workers
    if n_workers is None:
        n_workers = max(1, cpu_count() // 4)
    logger.info(f"Using {n_workers} parallel workers (out of {cpu_count()} cores)")

    # Classify genomes if genomes_dir is provided
    genome_classifications = {}
    if genomes_dir and os.path.isdir(genomes_dir):
        logger.info(f"Classifying genomes in {genomes_dir}...")
        genome_classifications = classify_genomes_in_directory(
            genomes_dir=genomes_dir,
            genome_hash_path=genome_hash_path
        )

    # Create media (just for info logging)
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
            model_name = os.path.splitext(os.path.basename(model_path))[0]
            genome_id = get_genome_id_for_model(model_name, genome_model_mapping)
            genome_class = genome_classifications.get(genome_id) if genome_id else None
            class_str = genome_class if genome_class else "unknown"
            print(f"  {os.path.basename(model_path)} -> {class_str}")
        return {}

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Prepare arguments for workers
    # Include genome classification for each model
    total = len(all_models)
    worker_args = []
    for model_path in all_models:
        model_name = os.path.splitext(os.path.basename(model_path))[0]
        genome_id = get_genome_id_for_model(model_name, genome_model_mapping)
        genome_class = genome_classifications.get(genome_id) if genome_id else None
        worker_args.append(
            (model_path, output_dir, gapfilling_mode, atp_safe, genome_class)
        )

    # Process models in parallel
    logger.info(f"\nStarting parallel gapfilling of {total} models...")
    results = {}

    with Pool(processes=n_workers) as pool:
        for i, (model_name, output_path) in enumerate(
            pool.imap_unordered(_gapfill_worker, worker_args), 1
        ):
            results[model_name] = output_path
            status = "OK" if output_path else "FAILED"
            logger.info(f"[{i}/{total}] {model_name}: {status}")

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
    parser.add_argument(
        '--workers', type=int, default=None,
        help='Number of parallel workers (default: 1/4 of CPU cores)'
    )
    parser.add_argument(
        '--genomes-dir', default='genomes',
        help='Directory containing genome files for classification (default: genomes)'
    )
    parser.add_argument(
        '--genome-hash', default='Sludge/datacache/genome_hash.json',
        help='Path to genome_hash.json with taxonomy data'
    )
    parser.add_argument(
        '--no-genome-classification', action='store_true',
        help='Disable genome-based template selection'
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
        dry_run=args.dry_run,
        n_workers=args.workers,
        genomes_dir=None if args.no_genome_classification else args.genomes_dir,
        genome_hash_path=None if args.no_genome_classification else args.genome_hash
    )

    # Return exit code based on success
    if results:
        failed = sum(1 for v in results.values() if v is None)
        sys.exit(1 if failed > 0 else 0)


if __name__ == '__main__':
    main()
