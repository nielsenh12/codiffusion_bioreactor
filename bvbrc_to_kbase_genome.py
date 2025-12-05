#!/usr/bin/env python3
"""
BV-BRC to KBase Genome Converter

This script can:
1. Fetch genome data from the BV-BRC (formerly PATRIC) API
2. Create KBase genome objects from local genome JSON files
3. Merge multiple genomes into synthetic genomes

Converts to KBase Genome object format. The resulting JSON is saved locally.

Usage:
    # From BV-BRC API:
    python bvbrc_to_kbase_genome.py --api <genome_id> [--output output_file]

    # From local genome file:
    python bvbrc_to_kbase_genome.py --local <genome_json_file> [--output output_file]

    # From local features directory (features/{genome_id}.json):
    python bvbrc_to_kbase_genome.py --features <genome_id> [--features-dir features] [--genomes-dir genomes]

    # Create synthetic genome from multiple local genomes:
    python bvbrc_to_kbase_genome.py --synthetic <asv_id> --genomes genome1.json,genome2.json [--output output_file]

Example:
    python bvbrc_to_kbase_genome.py --api 1110693.3 --output genome_output.json
    python bvbrc_to_kbase_genome.py --local my_genome.json
    python bvbrc_to_kbase_genome.py --features 511145.183 --features-dir features --genomes-dir genomes
    python bvbrc_to_kbase_genome.py --synthetic ASV_001 --genomes g1.json,g2.json,g3.json
"""

import sys
import os
import json
import requests
import hashlib
from datetime import datetime
from typing import Dict, List, Any, Optional
from collections import defaultdict
import argparse


class BVBRCToKBaseConverter:
    """Converter for BV-BRC genome data to KBase format"""

    def __init__(self, genome_id: str, base_url: str = "https://www.patricbrc.org/api"):
        self.genome_id = genome_id
        self.base_url = base_url
        self.session = requests.Session()
        self.session.verify = False  # Disable SSL verification as in Perl code

        # Suppress SSL warnings
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    def fetch_genome_metadata(self) -> Dict[str, Any]:
        """Fetch genome metadata from BV-BRC API"""
        url = f"{self.base_url}/genome/?eq(genome_id,{self.genome_id})&http_accept=application/json"

        print(f"Fetching genome metadata for {self.genome_id}...")
        response = self.session.get(url)
        response.raise_for_status()

        data = response.json()
        if not data:
            raise ValueError(f"No genome found with ID {self.genome_id}")

        return data[0]

    def fetch_genome_sequences(self) -> List[Dict[str, Any]]:
        """Fetch genome sequences (contigs) from BV-BRC API"""
        url = f"{self.base_url}/genome_sequence/?eq(genome_id,{self.genome_id})&http_accept=application/json"

        print(f"Fetching genome sequences...")
        response = self.session.get(url)
        response.raise_for_status()

        return response.json()

    def fetch_genome_features(self) -> List[Dict[str, Any]]:
        """Fetch all genome features from BV-BRC API (paginated)"""
        features = []
        start = 0
        limit = 10000

        print(f"Fetching genome features...")
        while True:
            url = f"{self.base_url}/genome_feature/?eq(genome_id,{self.genome_id})&http_accept=application/json&limit({limit},{start})"

            response = self.session.get(url)
            response.raise_for_status()

            batch = response.json()
            if not batch:
                break

            features.extend(batch)
            print(f"  Retrieved {len(features)} features so far...")
            start += limit

            if len(batch) < limit:
                break

        print(f"Total features retrieved: {len(features)}")
        return features

    def fetch_feature_sequences(self, md5_hashes: List[str]) -> Dict[str, Dict[str, Any]]:
        """Fetch feature sequences by MD5 hash (batched)"""
        sequences = {}
        batch_size = 100

        print(f"Fetching feature sequences for {len(md5_hashes)} unique sequences...")

        for i in range(0, len(md5_hashes), batch_size):
            batch = md5_hashes[i:i+batch_size]
            md5_list = ",".join(batch)
            url = f"{self.base_url}/feature_sequence/?in(md5,({md5_list}))&http_accept=application/json"

            try:
                response = self.session.get(url)
                response.raise_for_status()

                for seq_data in response.json():
                    md5 = seq_data.get('md5')
                    seq_type = seq_data.get('sequence_type', '')
                    sequence = seq_data.get('sequence', '')

                    if md5:
                        if md5 not in sequences:
                            sequences[md5] = {}
                        sequences[md5][seq_type] = sequence

                print(f"  Retrieved {len(sequences)} sequences so far...")
            except Exception as e:
                print(f"  Warning: Failed to fetch batch: {e}")
                continue

        return sequences

    def calculate_md5(self, sequence: str) -> str:
        """Calculate MD5 hash of a sequence"""
        return hashlib.md5(sequence.encode()).hexdigest()

    def build_kbase_genome(self) -> Dict[str, Any]:
        """Build complete KBase genome object"""

        # Fetch all data from BV-BRC
        genome_meta = self.fetch_genome_metadata()
        contigs = self.fetch_genome_sequences()
        features_data = self.fetch_genome_features()

        # Parse taxonomy
        taxonomy = genome_meta.get('taxon_lineage_names', [])
        taxonomy_str = "; ".join(taxonomy) if taxonomy else genome_meta.get('genome_name', '')

        # Determine domain
        domain = "Bacteria"
        if taxonomy:
            first_level = taxonomy[0].lower()
            if 'archaea' in first_level:
                domain = "Archaea"
            elif 'eukaryota' in first_level or 'eukarya' in first_level:
                domain = "Eukaryota"

        # Process contigs
        sorted_contigs = sorted(contigs, key=lambda x: x.get('accession', ''))
        contig_ids = []
        contig_lengths = []
        contig_sequences = []
        total_dna_size = 0

        for contig in sorted_contigs:
            contig_id = contig.get('accession', contig.get('sequence_id', ''))
            sequence = contig.get('sequence', '')
            length = len(sequence)

            contig_ids.append(contig_id)
            contig_lengths.append(length)
            contig_sequences.append(sequence)
            total_dna_size += length

        # Calculate genome MD5 (from sorted contig sequences)
        genome_md5 = self.calculate_md5("".join(contig_sequences))

        # Create contig ID mapping for features
        contig_map = {c.get('sequence_id', ''): c.get('accession', c.get('sequence_id', ''))
                      for c in contigs}

        # Collect MD5 hashes for feature sequences
        md5_hashes = set()
        for feature in features_data:
            if feature.get('na_sequence_md5'):
                md5_hashes.add(feature['na_sequence_md5'])
            if feature.get('aa_sequence_md5'):
                md5_hashes.add(feature['aa_sequence_md5'])

        # Fetch feature sequences
        sequences = self.fetch_feature_sequences(list(md5_hashes))

        # Process features
        kbase_features = []
        non_coding_features = []
        feature_counts = defaultdict(int)
        ontologies = {
            'SSO': {},
            'RefSeq': {},
            'FIGFAM': {},
            'PGFAM': {},
            'PLFAM': {},
            'GO': {}
        }

        print("Processing features...")
        for idx, feature in enumerate(features_data):
            kbase_feature = self._convert_feature(
                feature, idx, contig_map, sequences, ontologies
            )

            if kbase_feature:
                feature_type = kbase_feature['type']
                feature_counts[feature_type] += 1

                # Categorize feature
                if feature_type in ['CDS', 'gene', 'protein_encoding_gene']:
                    kbase_features.append(kbase_feature)
                    if feature_type in ['CDS', 'protein_encoding_gene']:
                        feature_counts['protein_encoding_gene'] += 1
                else:
                    non_coding_features.append(kbase_feature)
                    if feature_type not in ['CDS', 'gene']:
                        feature_counts['non-protein_encoding_gene'] += 1

        # Create CDS features for all genes
        cdss = self._create_cds_features(kbase_features)

        # Build final genome object
        genome = {
            'id': self.genome_id,
            'scientific_name': genome_meta.get('genome_name', ''),
            'domain': domain,
            'taxonomy': taxonomy_str,
            'genetic_code': int(genome_meta.get('genetic_code', 11)),
            'dna_size': total_dna_size,
            'num_contigs': len(contig_ids),
            'contig_ids': contig_ids,
            'contig_lengths': contig_lengths,
            'gc_content': float(genome_meta.get('gc_content', 0.5)),
            'md5': genome_md5,
            'molecule_type': 'DNA',
            'source': 'PATRIC',
            'source_id': self.genome_id,
            'assembly_ref': '',  # Blank as requested
            'external_source_origination_date': genome_meta.get('completion_date',
                                                                datetime.now().isoformat()),
            'notes': f'Imported from BV-BRC on {datetime.now().isoformat()}',
            'features': kbase_features,
            'non_coding_features': non_coding_features,
            'cdss': cdss,
            'mrnas': [],
            'feature_counts': dict(feature_counts),
            'publications': [],
            'genome_tiers': ['ExternalDB', 'User'],
            'warnings': [],
            'taxon_ref': '',  # Would need to be looked up in KBase
        }

        print(f"Genome object created successfully!")
        print(f"  - Features: {len(kbase_features)}")
        print(f"  - Non-coding features: {len(non_coding_features)}")
        print(f"  - CDS features: {len(cdss)}")
        print(f"  - Contigs: {len(contig_ids)}")
        print(f"  - DNA size: {total_dna_size:,} bp")

        return genome

    def _convert_feature(self, feature: Dict[str, Any], index: int,
                        contig_map: Dict[str, str], sequences: Dict[str, Dict[str, Any]],
                        ontologies: Dict[str, Dict]) -> Optional[Dict[str, Any]]:
        """Convert a BV-BRC feature to KBase format"""

        feature_type = feature.get('feature_type', 'gene')
        patric_id = feature.get('patric_id', '')

        # Get sequences
        na_md5 = feature.get('na_sequence_md5', '')
        aa_md5 = feature.get('aa_sequence_md5', '')

        na_sequence = sequences.get(na_md5, {}).get('dna', '') if na_md5 else ''
        aa_sequence = sequences.get(aa_md5, {}).get('protein', '') if aa_md5 else ''

        # Get contig ID
        sequence_id = feature.get('sequence_id', '')
        contig_id = contig_map.get(sequence_id, sequence_id)

        # Build location
        start = feature.get('start', 0)
        strand = feature.get('strand', '+')
        length = feature.get('na_length', len(na_sequence))
        location = [[contig_id, start, strand, length]]

        # Build functions list
        functions = []
        product = feature.get('product', '')
        if product:
            functions.append(product)

        # Build aliases
        aliases = [['PATRIC_id', patric_id]]

        # Add RefSeq data
        refseq_locus_tag = feature.get('refseq_locus_tag', '')
        gene_name = feature.get('gene', '')
        if refseq_locus_tag:
            aliases.append(['RefSeq_locus_tag', refseq_locus_tag])
        if gene_name:
            aliases.append(['gene_name', gene_name])

        # Build feature object
        feature_id = f"{self.genome_id}_{index}"
        kbase_feature = {
            'id': feature_id,
            'type': feature_type,
            'location': location,
            'functions': functions,
            'aliases': aliases,
            'dna_sequence': na_sequence,
            'dna_sequence_length': len(na_sequence),
            'md5': self.calculate_md5(na_sequence) if na_sequence else '',
        }

        # Add protein data if available
        if aa_sequence:
            kbase_feature['protein_translation'] = aa_sequence
            kbase_feature['protein_translation_length'] = len(aa_sequence)
            kbase_feature['protein_md5'] = self.calculate_md5(aa_sequence)

        # Collect ontology terms
        if product:
            if feature_id not in ontologies['SSO']:
                ontologies['SSO'][feature_id] = []
            ontologies['SSO'][feature_id].append(product)

            if feature_id not in ontologies['RefSeq']:
                ontologies['RefSeq'][feature_id] = []
            ontologies['RefSeq'][feature_id].append(product)

        # Add family IDs
        for family_type, ont_key in [('figfam_id', 'FIGFAM'),
                                      ('pgfam_id', 'PGFAM'),
                                      ('plfam_id', 'PLFAM')]:
            family_id = feature.get(family_type, '')
            if family_id:
                if feature_id not in ontologies[ont_key]:
                    ontologies[ont_key][feature_id] = []
                ontologies[ont_key][feature_id].append(family_id)

        # Add GO terms
        go_terms = feature.get('go', '').split(',') if feature.get('go') else []
        if go_terms:
            if feature_id not in ontologies['GO']:
                ontologies['GO'][feature_id] = []
            ontologies['GO'][feature_id].extend([g.strip() for g in go_terms if g.strip()])

        return kbase_feature

    def _create_cds_features(self, features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Create CDS features for all protein-coding genes"""
        cdss = []

        for feature in features:
            if feature.get('protein_translation'):  # Only for protein-coding genes
                cds = feature.copy()
                cds['id'] = f"{feature['id']}_CDS_1"
                cds['type'] = 'CDS'
                cds['parent_gene'] = feature['id']

                # Add CDS reference to parent gene
                feature['cdss'] = [cds['id']]

                cdss.append(cds)

        return cdss

    def save_genome(self, genome: Dict[str, Any], output_file: str):
        """Save genome object to JSON file"""
        print(f"\nSaving genome to {output_file}...")

        with open(output_file, 'w') as f:
            json.dump(genome, f, indent=2)

        print(f"Genome saved successfully!")
        print(f"File size: {len(json.dumps(genome)) / 1024 / 1024:.2f} MB")


class LocalGenomeConverter:
    """
    Converter for creating KBase genome objects from local feature JSON files.
    Adapted from the BV_BRC-Copy1.ipynb notebook code.
    """

    def __init__(self):
        pass

    def load_genome_json(self, genome_file: str) -> Dict[str, Any]:
        """Load a genome from a local JSON file"""
        # print(f"Loading genome from {genome_file}...")
        with open(genome_file, 'r') as f:
            genome = json.load(f)
        return genome

    def load_template_genome(self, template_file: str) -> Dict[str, Any]:
        """Load a template genome JSON file"""
        # print(f"Loading template genome from {template_file}...")
        with open(template_file, 'r') as f:
            template = json.load(f)
        return template

    def validate_kbase_genome(self, genome: Dict[str, Any]) -> Dict[str, Any]:
        """
        Ensure genome has all required KBase fields with proper defaults
        """
        # Ensure all required fields exist
        defaults = {
            'id': genome.get('id', 'unknown'),
            'scientific_name': genome.get('scientific_name', 'Unknown organism'),
            'domain': genome.get('domain', 'Bacteria'),
            'taxonomy': genome.get('taxonomy', ''),
            'genetic_code': int(genome.get('genetic_code', 11)),
            'dna_size': int(genome.get('dna_size', 0)),
            'num_contigs': int(genome.get('num_contigs', 0)),
            'contig_ids': genome.get('contig_ids', []),
            'contig_lengths': genome.get('contig_lengths', []),
            'gc_content': float(genome.get('gc_content', 0.5)),
            'md5': genome.get('md5', ''),
            'molecule_type': genome.get('molecule_type', 'DNA'),
            'source': genome.get('source', 'User'),
            'source_id': genome.get('source_id', ''),
            'assembly_ref': genome.get('assembly_ref', ''),
            'features': genome.get('features', []),
            'non_coding_features': genome.get('non_coding_features', []),
            'cdss': genome.get('cdss', []),
            'mrnas': genome.get('mrnas', []),
            'feature_counts': genome.get('feature_counts', {}),
            'publications': genome.get('publications', []),
            'genome_tiers': genome.get('genome_tiers', ['User']),
            'warnings': genome.get('warnings', []),
            'taxon_ref': genome.get('taxon_ref', ''),
        }

        # Add optional fields if present
        for key in ['external_source_origination_date', 'notes', 'ontologies_present',
                    'ontology_events', 'genome_type']:
            if key in genome:
                defaults[key] = genome[key]

        return defaults

    def aggregate_taxonomies(
        self,
        genomes: List[Dict[str, Any]],
        asv_id: str,
        output_dir: str = "ASVset_taxonomies"
    ) -> tuple[str, Dict[str, List[str]]]:
        """
        Aggregate taxonomies from multiple genomes and select most prevalent taxonomy.

        Args:
            genomes: List of genome dictionaries
            asv_id: Identifier for the ASV/synthetic genome
            output_dir: Directory to save taxonomy JSON (default: ASVset_taxonomies)

        Returns:
            Tuple of (consensus_taxonomy_string, taxonomy_dict)

        The taxonomy_dict format:
        {
            "Kingdom": ["Bacteria", "Bacteria", "Bacteria"],
            "Phylum": ["Proteobacteria", "Firmicutes", "Proteobacteria"],
            ...
        }
        """
        from collections import Counter

        # Standard taxonomic levels
        tax_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

        # Collect all taxonomies
        all_taxonomies = []
        for genome in genomes:
            taxonomy = genome.get('taxonomy', '')
            if taxonomy:
                all_taxonomies.append(taxonomy)

        if not all_taxonomies:
            print("  Warning: No taxonomies found in source genomes")
            return "Unknown", {}

        # Parse taxonomies into levels
        taxonomy_by_level = {level: [] for level in tax_levels}

        for taxonomy_str in all_taxonomies:
            # Split by semicolon or other common delimiters
            parts = [p.strip() for p in taxonomy_str.replace(';', '|').split('|')]

            # Assign to levels (assuming order matches standard levels)
            for i, part in enumerate(parts):
                if i < len(tax_levels) and part:
                    taxonomy_by_level[tax_levels[i]].append(part)

        # Find most common taxonomy at each level
        consensus_taxonomy = []
        for level in tax_levels:
            if taxonomy_by_level[level]:
                # Count occurrences
                counts = Counter(taxonomy_by_level[level])
                # Get most common
                most_common = counts.most_common(1)[0][0]
                consensus_taxonomy.append(most_common)
            else:
                # No data at this level, stop here
                break

        # Build consensus taxonomy string
        consensus_str = "; ".join(consensus_taxonomy)

        # Build output dictionary (only include levels with data)
        output_dict = {
            level: taxonomy_by_level[level]
            for level in tax_levels
            if taxonomy_by_level[level]
        }

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Save to JSON
        output_file = os.path.join(output_dir, f"{asv_id}.json")
        with open(output_file, 'w') as f:
            json.dump(output_dict, f, indent=2)

        print(f"  Taxonomy saved to {output_file}")
        # print(f"  Consensus taxonomy: {consensus_str}")

        return consensus_str, output_dict

    def create_synthetic_genome(
        self,
        asv_id: str,
        genome_files: Optional[List[str]] = None,
        genome_ids: Optional[List[str]] = None,
        genomes: Optional[List[Dict[str, Any]]] = None,
        features_dir: str = "features",
        genomes_dir: str = "genomes",
        metadata_dir: str = "genome_metadata",
        taxonomy: Optional[str] = None,
        template_file: Optional[str] = "model_inputs/TemplateGenomes.json",
        save_taxonomy: bool = True,
        taxonomy_output_dir: str = "ASVset_taxonomies",
        genome_output_dir: str = "genome_objects"
    ) -> Dict[str, Any]:
        """
        Create a synthetic genome from multiple source genomes.
        Adapted from the notebook code under "# Create Synthetic genomes; ## Local version"

        Args:
            asv_id: Identifier for the synthetic genome (e.g., ASV or genus name)
            genome_files: List of paths to genome JSON files to merge (optional)
            genome_ids: List of genome IDs to load from features/{id}.json and genomes/{id}.fna (optional)
            genomes: List of pre-loaded genome dictionaries (optional)
            features_dir: Directory containing features JSON files (default: "features")
            genomes_dir: Directory containing genome FASTA files (default: "genomes")
            metadata_dir: Directory containing genome metadata JSON files (default: "genome_metadata")
            taxonomy: Optional taxonomy string for the synthetic genome.
                     If not provided and save_taxonomy=True, will use consensus taxonomy.
            template_file: Optional path to template genome JSON
            save_taxonomy: If True, aggregates and saves taxonomy from source genomes
                          (default: True)
            taxonomy_output_dir: Directory to save taxonomy JSON
                                (default: "ASVset_taxonomies")

        Returns:
            KBase Genome object dictionary

        Note:
            Provide exactly one of: genome_files, genome_ids, or genomes

        Side Effects:
            If save_taxonomy=True, creates {taxonomy_output_dir}/{asv_id}.json with format:
            {
                "Kingdom": ["Bacteria", "Bacteria", "Bacteria"],
                "Phylum": ["Proteobacteria", "Firmicutes", "Proteobacteria"],
                ...
            }
        """
        # Validate input
        inputs_provided = sum([
            genome_files is not None,
            genome_ids is not None,
            genomes is not None
        ])
        if inputs_provided == 0:
            raise ValueError("Must provide one of: genome_files, genome_ids, or genomes")
        if inputs_provided > 1:
            raise ValueError("Provide only one of: genome_files, genome_ids, or genomes")

        print(f"\nCreating synthetic genome: {asv_id}")

        # Load source genomes based on input type
        source_genomes = []
        source_ids = []

        if genome_ids:
            # print(f"Loading {len(genome_ids)} genomes from features directory...")
            for genome_id in genome_ids:
                try:
                    genome = self.load_genome_from_features_dir(
                        genome_id=genome_id,
                        features_dir=features_dir,
                        genomes_dir=genomes_dir,
                        metadata_dir=metadata_dir
                    )
                    source_genomes.append(genome)
                    source_ids.append(genome_id)
                except Exception as e:
                    print(f"  Warning: Failed to load {genome_id}: {e}")
                    continue

        elif genome_files:
            # print(f"Loading {len(genome_files)} genomes from JSON files...")
            for genome_file in genome_files:
                try:
                    genome = self.load_genome_json(genome_file)
                    source_genomes.append(genome)
                    source_ids.append(genome_file)
                except Exception as e:
                    print(f"  Warning: Failed to load {genome_file}: {e}")
                    continue

        elif genomes:
            # print(f"Using {len(genomes)} pre-loaded genomes...")
            source_genomes = genomes
            source_ids = [g.get('id', f'genome_{i}') for i, g in enumerate(genomes)]

        if not source_genomes:
            raise ValueError("No valid source genomes could be loaded")

        print(f"Successfully loaded {len(source_genomes)} source genomes")

        # Load or create template
        if template_file and os.path.exists(template_file):
            template_genome = self.load_template_genome(template_file)
        else:
            # Create minimal template
            template_genome = {
                'id': asv_id,
                'scientific_name': taxonomy or asv_id,
                'taxonomy': taxonomy or '',
                'domain': 'Bacteria',
                'genetic_code': 11,
                'dna_size': 0,
                'num_contigs': 0,
                'contig_ids': [],
                'contig_lengths': [],
                'gc_content': 0.0,
                'md5': '',
                'molecule_type': 'DNA',
                'source': 'Synthetic',
                'source_id': '|'.join(source_ids),
                'assembly_ref': '',
                'features': [],
                'non_coding_features': [],
                'cdss': [],
                'mrnas': [],
                'feature_counts': {},
                'publications': [],
                'genome_tiers': ['User'],
                'warnings': ['Synthetic genome created by merging multiple source genomes'],
                'taxon_ref': '',
            }

        # Calculate average GC content
        gc_contents = [float(g.get('gc_content', 0.5)) for g in source_genomes if 'gc_content' in g]
        if not gc_contents:
            raise ValueError("No valid source genomes could be loaded")

        # Aggregate taxonomies from source genomes
        if save_taxonomy:
            print(f"Aggregating taxonomies from {len(source_genomes)} source genomes...")
            consensus_taxonomy, taxonomy_dict = self.aggregate_taxonomies(
                genomes=source_genomes,
                asv_id=asv_id,
                output_dir=taxonomy_output_dir
            )
            # Use consensus taxonomy if not explicitly provided
            if not taxonomy:
                taxonomy = consensus_taxonomy
                template_genome['taxonomy'] = consensus_taxonomy
                template_genome['scientific_name'] = consensus_taxonomy.split(';')[-1].strip() if consensus_taxonomy else asv_id

                # Set domain from taxonomy
                if consensus_taxonomy:
                    first_level = consensus_taxonomy.split(';')[0].strip().lower()
                    if 'archaea' in first_level:
                        template_genome['domain'] = 'Archaea'
                    elif 'eukaryota' in first_level or 'eukarya' in first_level:
                        template_genome['domain'] = 'Eukaryota'
                    else:
                        template_genome['domain'] = 'Bacteria'
        elif taxonomy:
            # User provided taxonomy
            template_genome['taxonomy'] = taxonomy
            template_genome['scientific_name'] = taxonomy.split(';')[-1].strip() if taxonomy else asv_id

        # Calculate average GC content
        if gc_contents:
            from numpy import mean
            template_genome['gc_content'] = float(mean(gc_contents))

        # Track unique functions and features
        functions = {}  # function -> feature info
        features = {}  # feature_id -> feature dict
        md5_list = []

        print(f"Processing features from {len(source_genomes)} genomes...")

        # Iterate through source genomes and collect unique functions
        for genome_idx, source_genome in enumerate(source_genomes):
            genome_id = source_genome.get('id', f'genome_{genome_idx}')
            genome_functions = {}  # Track functions in this genome

            # Process each feature from this genome
            for source_feature in source_genome.get('features', []):
                if 'functions' not in source_feature or not source_feature['functions']:
                    continue

                for function in source_feature['functions']:
                    # Check if this function is new to the synthetic genome
                    if function not in functions:
                        # Create new feature ID
                        feature_id = f"{asv_id}_{len(template_genome['contig_ids']) + 1}"

                        # Update contig information
                        template_genome['contig_ids'].append(f"{feature_id}.contig")
                        dna_length = len(source_feature.get('dna_sequence', ''))
                        template_genome['contig_lengths'].append(dna_length)
                        template_genome['num_contigs'] += 1
                        template_genome['dna_size'] += dna_length

                        # Calculate protein MD5
                        protein_seq = source_feature.get('protein_translation', '')
                        protein_md5 = hashlib.md5(protein_seq.encode()).hexdigest() if protein_seq else ''
                        if protein_md5:
                            md5_list.append(protein_md5)

                        # Track function with probability (starts at 1)
                        functions[function] = {
                            'feature_id': feature_id,
                            'probability': 1
                        }

                        # Create feature object
                        features[feature_id] = {
                            'id': feature_id,
                            'type': source_feature.get('type', 'gene'),
                            'aliases': source_feature.get('aliases', [])[:],  # Copy aliases
                            'cdss': [f"{feature_id}.CDS"],
                            'functions': [function],
                            'dna_sequence': source_feature.get('dna_sequence', ''),
                            'dna_sequence_length': dna_length,
                            'location': [[f"{feature_id}.contig", 1, "+", dna_length]],
                            'md5': hashlib.md5(source_feature.get('dna_sequence', '').encode()).hexdigest(),
                            'ontology_terms': source_feature.get('ontology_terms', {}),
                            'protein_md5': protein_md5,
                            'protein_translation': protein_seq,
                            'protein_translation_length': len(protein_seq),
                            'warnings': []
                        }

                        # Create CDS feature
                        cds_feature = features[feature_id].copy()
                        del cds_feature['cdss']
                        cds_feature['id'] = f"{feature_id}.CDS"
                        cds_feature['type'] = 'CDS'
                        cds_feature['parent_gene'] = feature_id
                        template_genome['cdss'].append(cds_feature)

                        # Add to features list
                        template_genome['features'].append(features[feature_id])

                    elif function not in genome_functions:
                        # Function exists but not counted in this genome yet
                        # Increment probability counter
                        functions[function]['probability'] += 1

                    # Mark function as seen in this genome
                    genome_functions[function] = True

        # Normalize probabilities
        num_genomes = len(source_genomes)
        for function in functions:
            functions[function]['probability'] /= num_genomes

        # Update feature counts
        template_genome['feature_counts'] = {
            'CDS': len(template_genome['cdss']),
            'gene': len(template_genome['features']),
            'protein_encoding_gene': len(template_genome['features']),
            'non_coding_features': len(template_genome.get('non_coding_features', [])),
        }

        # Calculate genome MD5 from sorted protein MD5s
        md5_list.sort()
        genome_md5 = hashlib.md5(";".join(md5_list).encode()).hexdigest()
        template_genome['md5'] = genome_md5

        print(f"Synthetic genome created:")
        print(f"  - Total features: {len(template_genome['features'])}")
        print(f"  - Total CDS: {len(template_genome['cdss'])}")
        print(f"  - Total contigs: {template_genome['num_contigs']}")
        print(f"  - Total DNA size: {template_genome['dna_size']:,} bp")

        json.dump(template_genome, open(f"genome_objects/{asv_id}.json", 'w'))
        return template_genome

    def parse_fasta(self, fasta_file: str) -> Dict[str, str]:
        """
        Parse a FASTA file and return a dictionary of sequence_id -> sequence
        """
        sequences = {}
        current_id = None
        current_seq = []

        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get first part of header
                    current_seq = []
                else:
                    current_seq.append(line)

            # Save last sequence
            if current_id:
                sequences[current_id] = ''.join(current_seq)

        return sequences

    def load_genome_from_features_dir(
        self,
        genome_id: str,
        features_dir: str = "features",
        genomes_dir: str = "genomes",
        metadata_dir: str = "genome_metadata",
        taxonomy: Optional[str] = None,
        scientific_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Load genome from local BV-BRC files (three-file structure).

        Reads from:
        - genome_metadata/{genome_id}.json - Taxonomy, GC content, genome stats
        - genomes/{genome_id}.fna - Genome sequences in FASTA format
        - features/{genome_id}.json - Feature metadata from BV-BRC API

        Args:
            genome_id: The BV-BRC genome ID
            features_dir: Directory containing feature JSON files (default: "features")
            genomes_dir: Directory containing genome FASTA files (default: "genomes")
            metadata_dir: Directory containing genome metadata files (default: "genome_metadata")
            taxonomy: Optional taxonomy string (overrides metadata)
            scientific_name: Optional scientific name (overrides metadata)

        Returns:
            KBase Genome object dictionary
        """
        # print(f"\nLoading genome {genome_id} from local files...")

        # Construct file paths
        metadata_file = os.path.join(metadata_dir, f"{genome_id}.json")
        features_file = os.path.join(features_dir, f"{genome_id}.json")
        genome_file = os.path.join(genomes_dir, f"{genome_id}.fna")

        # Load genome metadata (preferred source for taxonomy and stats)
        metadata = {}
        if os.path.exists(metadata_file):
            # print(f"Loading metadata from {metadata_file}...")
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            # print(f"  Loaded metadata for {metadata.get('genome_name', genome_id)}")
        # else:
            # print(f"  Warning: Metadata file not found: {metadata_file}")

        # Load features
        if os.path.exists(features_file):
            # print(f"Loading features from {features_file}...")
            with open(features_file, 'r') as f:
                features_data = json.load(f)
            # print(f"  Loaded {len(features_data)} features")
        else:
            # print(f"  Warning: Features file not found: {features_file}")
            features_data = []

        # Load sequences if available
        sequences = {}
        if os.path.exists(genome_file):
            # print(f"Loading sequences from {genome_file}...")
            sequences = self.parse_fasta(genome_file)
            # print(f"  Loaded {len(sequences)} contig sequences")
        # else:
            # print(f"  Warning: Genome file not found: {genome_file}")

        # Calculate contig information
        contig_ids = sorted(sequences.keys())
        contig_lengths = [len(sequences[cid]) for cid in contig_ids]
        total_dna_size = sum(contig_lengths)

        # Get GC content from metadata (preferred) or calculate from sequences
        gc_content = 0.5  # Default
        if metadata.get('gc_content'):
            # Metadata has GC content as percentage (0-100), convert to fraction
            gc_content = float(metadata['gc_content']) / 100.0
        elif sequences:
            # Calculate GC content from sequences if not in metadata
            all_seq = ''.join(sequences.values()).upper()
            g_count = all_seq.count('G')
            c_count = all_seq.count('C')
            total = len(all_seq)
            if total > 0:
                gc_content = (g_count + c_count) / total

        # Calculate genome MD5 from sorted contig sequences
        genome_md5 = ''
        if sequences:
            sorted_seqs = [sequences[cid] for cid in contig_ids]
            genome_md5 = hashlib.md5(''.join(sorted_seqs).encode()).hexdigest()

        # Build taxonomy string from metadata (if not provided)
        if not taxonomy and metadata:
            # Use taxon_lineage_names if available (most complete)
            if 'taxon_lineage_names' in metadata and len(metadata['taxon_lineage_names']) > 1:
                # Skip first element ("cellular organisms") and join the rest
                taxonomy = '; '.join(metadata['taxon_lineage_names'][1:])
            else:
                # Construct from individual taxonomy fields
                tax_parts = []
                for level in ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                             'family', 'genus', 'species']:
                    if level in metadata and metadata[level]:
                        tax_parts.append(metadata[level])
                taxonomy = '; '.join(tax_parts) if tax_parts else ''

        # Get scientific name from metadata (if not provided)
        if not scientific_name and metadata:
            scientific_name = (metadata.get('genome_name') or
                             metadata.get('species') or
                             metadata.get('genus', genome_id))

        # Determine domain from metadata
        domain = 'Bacteria'  # Default
        if metadata:
            superkingdom = metadata.get('superkingdom', '').lower()
            if 'archaea' in superkingdom:
                domain = 'Archaea'
            elif 'eukaryot' in superkingdom or 'eukarya' in superkingdom:
                domain = 'Eukaryota'
            elif 'bacteria' in superkingdom:
                domain = 'Bacteria'

        # Use metadata for contig counts if sequences not available
        if not sequences and metadata:
            total_dna_size = metadata.get('genome_length', 0)
            num_contigs = metadata.get('contigs', 0)

        # Process features into KBase format
        kbase_features = []
        non_coding_features = []
        feature_counts = defaultdict(int)

        # print(f"Processing features...")
        for idx, feature in enumerate(features_data):
            kbase_feature = self._convert_bvbrc_feature_to_kbase(
                feature, idx, genome_id, sequences
            )

            if kbase_feature:
                feature_type = kbase_feature['type']
                feature_counts[feature_type] += 1

                # Categorize feature
                if feature_type in ['CDS', 'gene', 'protein_encoding_gene']:
                    kbase_features.append(kbase_feature)
                    if feature_type in ['CDS', 'protein_encoding_gene']:
                        feature_counts['protein_encoding_gene'] += 1
                else:
                    non_coding_features.append(kbase_feature)
                    feature_counts['non-protein_encoding_gene'] += 1

        # Create CDS features for all protein-coding genes
        cdss = []
        for feature in kbase_features:
            if feature.get('protein_translation'):  # Has protein sequence
                cds = feature.copy()
                cds['id'] = f"{feature['id']}_CDS"
                cds['type'] = 'CDS'
                cds['parent_gene'] = feature['id']
                feature['cdss'] = [cds['id']]
                cdss.append(cds)

        # Build genome object
        genome = {
            'id': genome_id,
            'scientific_name': scientific_name or genome_id,
            'domain': domain,
            'taxonomy': taxonomy or '',
            'genetic_code': 11,
            'dna_size': total_dna_size,
            'num_contigs': len(contig_ids) if sequences else metadata.get('contigs', 0),
            'contig_ids': contig_ids,
            'contig_lengths': contig_lengths,
            'gc_content': gc_content,
            'md5': genome_md5,
            'molecule_type': 'DNA',
            'source': 'PATRIC',
            'source_id': genome_id,
            'assembly_ref': '',
            'external_source_origination_date': metadata.get('completion_date', datetime.now().isoformat()),
            'notes': f'Imported from local BV-BRC files on {datetime.now().isoformat()}',
            'features': kbase_features,
            'non_coding_features': non_coding_features,
            'cdss': cdss,
            'mrnas': [],
            'feature_counts': dict(feature_counts),
            'publications': [],
            'genome_tiers': ['User'],
            'warnings': [],
            'taxon_ref': '',
        }

        # print(f"Genome object created:")
        # print(f"  - Features: {len(kbase_features)}")
        # print(f"  - Non-coding features: {len(non_coding_features)}")
        # print(f"  - CDS features: {len(cdss)}")
        # print(f"  - Contigs: {len(contig_ids)}")
        # print(f"  - DNA size: {total_dna_size:,} bp")
        # print(f"  - GC content: {gc_content:.2%}")

        return genome

    def _convert_bvbrc_feature_to_kbase(
        self,
        feature: Dict[str, Any],
        index: int,
        genome_id: str,
        sequences: Dict[str, str]
    ) -> Optional[Dict[str, Any]]:
        """
        Convert a BV-BRC feature from features/{genome_id}.json to KBase format.

        BV-BRC feature structure:
        {
            "patric_id": "fig|{genome_id}.{type}.{num}",
            "product": "gene product description",
            "feature_type": "CDS|rRNA|tRNA|misc",
            "pgfam_id": "PGF00000001",
            "plfam_id": "PLF00000123",
            "figfam_id": "FIG00000456",
            "annotation": "PATRIC"
        }
        """
        patric_id = feature.get('patric_id', '')
        product = feature.get('product', '')
        feature_type = feature.get('feature_type', 'gene')

        # Build functions list
        functions = []
        if product:
            functions.append(product)

        # Build aliases
        aliases = [['PATRIC_id', patric_id]]

        # Add family IDs as aliases
        for family_type, ont_key in [('figfam_id', 'FIGFAM'),
                                      ('pgfam_id', 'PGFAM'),
                                      ('plfam_id', 'PLFAM')]:
            family_id = feature.get(family_type, '')
            if family_id:
                aliases.append([ont_key, family_id])

        # Create feature ID
        feature_id = f"{genome_id}_{index}"

        # Build basic feature object
        kbase_feature = {
            'id': feature_id,
            'type': feature_type,
            'location': [],  # Will be filled if we have sequence location data
            'functions': functions,
            'aliases': aliases,
            'dna_sequence': '',
            'dna_sequence_length': 0,
            'md5': '',
        }

        # Note: The features/{genome_id}.json files from BV-BRC API don't include
        # sequence data or location data. For complete genomes, you would need to:
        # 1. Query the full BV-BRC API with more fields, or
        # 2. Use annotation tools to map features to sequences

        return kbase_feature

    def create_fasta_from_genome(self, genome: Dict[str, Any], output_file: str):
        """
        Create a FASTA file from genome features.
        Adapted from the notebook code.
        """
        print(f"Writing FASTA file to {output_file}...")

        with open(output_file, 'w') as f:
            for feature in genome.get('features', []):
                feature_id = feature.get('id', 'unknown')
                dna_sequence = feature.get('dna_sequence', '')

                if dna_sequence:
                    f.write(f">{feature_id}\n{dna_sequence}\n")

        print(f"FASTA file created with {len(genome.get('features', []))} sequences")

    def save_genome(self, genome: Dict[str, Any], output_file: str):
        """Save genome object to JSON file"""
        print(f"\nSaving genome to {output_file}...")

        with open(output_file, 'w') as f:
            json.dump(genome, f, indent=2)

        print(f"Genome saved successfully!")
        print(f"File size: {len(json.dumps(genome)) / 1024 / 1024:.2f} MB")


# ============================================================================
# Convenience Functions for Notebook/Library Usage
# ============================================================================

def fetch_genome_from_api(genome_id: str, output_file: Optional[str] = None) -> Dict[str, Any]:
    """
    Convenience function to fetch a genome from BV-BRC API.

    Args:
        genome_id: BV-BRC genome ID
        output_file: Optional path to save JSON output

    Returns:
        KBase Genome object dictionary

    Example:
        >>> genome = fetch_genome_from_api('511145.183')
        >>> genome = fetch_genome_from_api('511145.183', 'ecoli.json')
    """
    converter = BVBRCToKBaseConverter(genome_id)
    genome = converter.build_kbase_genome()

    if output_file:
        converter.save_genome(genome, output_file)

    return genome


def load_genome_from_json(json_file: str, output_file: Optional[str] = None) -> Dict[str, Any]:
    """
    Convenience function to load and validate a genome from JSON file.

    Args:
        json_file: Path to genome JSON file
        output_file: Optional path to save validated JSON output

    Returns:
        KBase Genome object dictionary

    Example:
        >>> genome = load_genome_from_json('my_genome.json')
        >>> genome = load_genome_from_json('my_genome.json', 'validated.json')
    """
    converter = LocalGenomeConverter()
    genome = converter.load_genome_json(json_file)
    genome = converter.validate_kbase_genome(genome)

    if output_file:
        converter.save_genome(genome, output_file)

    return genome


def load_genome_from_features(
    genome_id: str,
    features_dir: str = "features",
    genomes_dir: str = "genomes",
    metadata_dir: str = "genome_metadata",
    taxonomy: Optional[str] = None,
    scientific_name: Optional[str] = None,
    output_file: Optional[str] = None
) -> Dict[str, Any]:
    """
    Convenience function to load genome from local BV-BRC feature files.

    Args:
        genome_id: BV-BRC genome ID
        features_dir: Directory containing features JSON files
        genomes_dir: Directory containing genome FASTA files
        metadata_dir: Directory containing genome metadata JSON files
        taxonomy: Optional taxonomy string (overrides metadata)
        scientific_name: Optional scientific name (overrides metadata)
        output_file: Optional path to save JSON output

    Returns:
        KBase Genome object dictionary

    Example:
        >>> genome = load_genome_from_features('511145.183')
        >>> genome = load_genome_from_features('511145.183',
        ...                                     taxonomy='Bacteria; Proteobacteria',
        ...                                     output_file='ecoli.json')
    """
    converter = LocalGenomeConverter()
    genome = converter.load_genome_from_features_dir(
        genome_id=genome_id,
        features_dir=features_dir,
        genomes_dir=genomes_dir,
        metadata_dir=metadata_dir,
        taxonomy=taxonomy,
        scientific_name=scientific_name
    )

    if output_file:
        converter.save_genome(genome, output_file)

    return genome


def aggregate_taxonomies(
    genomes: List[Dict[str, Any]],
    asv_id: str,
    output_dir: str = "ASVset_taxonomies"
) -> tuple[str, Dict[str, List[str]]]:
    """
    Convenience function to aggregate taxonomies from multiple genomes.

    Args:
        genomes: List of genome dictionaries
        asv_id: Identifier for the ASV/synthetic genome
        output_dir: Directory to save taxonomy JSON (default: ASVset_taxonomies)

    Returns:
        Tuple of (consensus_taxonomy_string, taxonomy_dict)

    Example:
        >>> genomes = [load_genome_from_json('g1.json'), load_genome_from_json('g2.json')]
        >>> consensus, tax_dict = aggregate_taxonomies(genomes, 'ASV_001')
        >>> print(consensus)
        'Bacteria; Proteobacteria; Gammaproteobacteria'
        >>> print(tax_dict['Phylum'])
        ['Proteobacteria', 'Proteobacteria', 'Firmicutes']
    """
    converter = LocalGenomeConverter()
    return converter.aggregate_taxonomies(genomes, asv_id, output_dir)


def create_synthetic_genome(
    asv_id: str,
    genome_files: Optional[List[str]] = None,
    genome_ids: Optional[List[str]] = None,
    genomes: Optional[List[Dict[str, Any]]] = None,
    features_dir: str = "features",
    genomes_dir: str = "genomes",
    metadata_dir: str = "genome_metadata",
    taxonomy: Optional[str] = None,
    template_file: Optional[str] = None,
    save_taxonomy: bool = True,
    taxonomy_output_dir: str = "ASVset_taxonomies",
    output_file: Optional[str] = None
) -> Dict[str, Any]:
    """
    Convenience function to create synthetic genome from multiple sources.

    Args:
        asv_id: Identifier for synthetic genome
        genome_files: List of paths to genome JSON files (optional)
        genome_ids: List of genome IDs to load from features/{id}.json (optional)
        genomes: List of pre-loaded genome dictionaries (optional)
        features_dir: Directory containing features JSON files (default: features/)
        genomes_dir: Directory containing genome FASTA files (default: genomes/)
        metadata_dir: Directory containing genome metadata JSON files (default: genome_metadata/)
        taxonomy: Optional taxonomy string (if not provided, uses consensus from sources)
        template_file: Optional template genome JSON path
        save_taxonomy: If True, saves taxonomy aggregation to JSON (default: True)
        taxonomy_output_dir: Directory for taxonomy JSON (default: ASVset_taxonomies)
        output_file: Optional path to save JSON output

    Returns:
        KBase Genome object dictionary

    Note:
        Provide exactly one of: genome_files, genome_ids, or genomes

    Example:
        >>> # From genome IDs (features/{id}.json and genomes/{id}.fna)
        >>> genome = create_synthetic_genome('ASV_001',
        ...                                   genome_ids=['511145.183', '83332.133'])

        >>> # From JSON files
        >>> genome = create_synthetic_genome('ASV_001',
        ...                                   genome_files=['g1.json', 'g2.json'])

        >>> # From pre-loaded genomes
        >>> genomes = [load_genome_from_features(id) for id in ['511145.183', '83332.133']]
        >>> genome = create_synthetic_genome('ASV_001', genomes=genomes)

        >>> # With explicit taxonomy
        >>> genome = create_synthetic_genome('ASV_001',
        ...                                   genome_ids=['511145.183', '83332.133'],
        ...                                   taxonomy='Bacteria; Firmicutes',
        ...                                   output_file='asv_001.json')
    """
    converter = LocalGenomeConverter()
    genome = converter.create_synthetic_genome(
        asv_id=asv_id,
        genome_files=genome_files,
        genome_ids=genome_ids,
        genomes=genomes,
        features_dir=features_dir,
        genomes_dir=genomes_dir,
        metadata_dir=metadata_dir,
        taxonomy=taxonomy,
        template_file=template_file,
        save_taxonomy=save_taxonomy,
        taxonomy_output_dir=taxonomy_output_dir
    )

    if output_file:
        converter.save_genome(genome, output_file)

    return genome


def save_genome_to_json(genome: Dict[str, Any], output_file: str):
    """
    Convenience function to save a genome object to JSON.

    Args:
        genome: KBase Genome object dictionary
        output_file: Path to save JSON output

    Example:
        >>> save_genome_to_json(genome, 'my_genome.json')
    """
    converter = LocalGenomeConverter()
    converter.save_genome(genome, output_file)


def create_fasta_from_genome(genome: Dict[str, Any], output_file: str):
    """
    Convenience function to create FASTA file from genome features.

    Args:
        genome: KBase Genome object dictionary
        output_file: Path to save FASTA output

    Example:
        >>> create_fasta_from_genome(genome, 'features.fasta')
    """
    converter = LocalGenomeConverter()
    converter.create_fasta_from_genome(genome, output_file)


# ============================================================================
# Command-Line Interface
# ============================================================================

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Convert genome data to KBase format from BV-BRC API or local files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Mode selection (mutually exclusive)
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--api', metavar='GENOME_ID',
                           help='Fetch genome from BV-BRC API using genome ID')
    mode_group.add_argument('--local', metavar='JSON_FILE',
                           help='Load genome from local JSON file')
    mode_group.add_argument('--features', metavar='GENOME_ID',
                           help='Load genome from local features/{genome_id}.json and genomes/{genome_id}.fna')
    mode_group.add_argument('--synthetic', metavar='ASV_ID',
                           help='Create synthetic genome from multiple local genomes')

    # Additional arguments
    parser.add_argument('--genomes', metavar='FILES',
                       help='Comma-separated list of genome JSON files (for --synthetic mode)')
    parser.add_argument('--genome-ids', metavar='IDS',
                       help='Comma-separated list of genome IDs to load from features/ (for --synthetic mode)')
    parser.add_argument('--template', metavar='FILE',
                       help='Template genome JSON file (optional, for --synthetic mode)')
    parser.add_argument('--taxonomy', metavar='STRING',
                       help='Taxonomy string (optional, for --synthetic and --features modes)')
    parser.add_argument('--scientific-name', metavar='STRING',
                       help='Scientific name (optional, for --features mode)')
    parser.add_argument('--features-dir', metavar='DIR', default='features',
                       help='Directory containing feature JSON files (default: features/)')
    parser.add_argument('--genomes-dir', metavar='DIR', default='genomes',
                       help='Directory containing genome FASTA files (default: genomes/)')
    parser.add_argument('--metadata-dir', metavar='DIR', default='genome_metadata',
                       help='Directory containing genome metadata JSON files (default: genome_metadata/)')
    parser.add_argument('--taxonomy-dir', metavar='DIR', default='ASVset_taxonomies',
                       help='Directory to save taxonomy aggregation JSON (default: ASVset_taxonomies/)')
    parser.add_argument('--no-taxonomy', action='store_true',
                       help='Disable taxonomy aggregation for synthetic genomes')
    parser.add_argument('--output', '-o', metavar='FILE',
                       help='Output JSON file path')
    parser.add_argument('--fasta', metavar='FILE',
                       help='Also create FASTA file with feature sequences')

    # Legacy support: if no arguments with dashes, assume first arg is genome ID
    if len(sys.argv) >= 2 and not sys.argv[1].startswith('-'):
        # Legacy mode: python script.py genome_id [output_file]
        genome_id = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else f"{genome_id}_genome.json"

        print(f"BV-BRC to KBase Genome Converter (Legacy Mode)")
        print(f"=" * 50)
        print(f"Genome ID: {genome_id}")
        print(f"Output file: {output_file}")
        print()

        try:
            converter = BVBRCToKBaseConverter(genome_id)
            genome = converter.build_kbase_genome()
            converter.save_genome(genome, output_file)

            print("\n" + "=" * 50)
            print("Conversion completed successfully!")
            return

        except Exception as e:
            print(f"\nError: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    args = parser.parse_args()

    print(f"KBase Genome Converter")
    print(f"=" * 50)

    try:
        # Mode 1: Fetch from BV-BRC API
        if args.api:
            genome_id = args.api
            output_file = args.output or f"{genome_id}_genome.json"

            print(f"Mode: Fetch from BV-BRC API")
            print(f"Genome ID: {genome_id}")
            print(f"Output file: {output_file}")
            print()

            converter = BVBRCToKBaseConverter(genome_id)
            genome = converter.build_kbase_genome()
            converter.save_genome(genome, output_file)

            if args.fasta:
                local_conv = LocalGenomeConverter()
                local_conv.create_fasta_from_genome(genome, args.fasta)

        # Mode 2: Load from local file
        elif args.local:
            genome_file = args.local
            output_file = args.output or f"{os.path.splitext(os.path.basename(genome_file))[0]}_kbase.json"

            print(f"Mode: Load from local file")
            print(f"Input file: {genome_file}")
            print(f"Output file: {output_file}")
            print()

            local_conv = LocalGenomeConverter()
            genome = local_conv.load_genome_json(genome_file)
            genome = local_conv.validate_kbase_genome(genome)
            local_conv.save_genome(genome, output_file)

            if args.fasta:
                local_conv.create_fasta_from_genome(genome, args.fasta)

        # Mode 3: Load from features directory
        elif args.features:
            genome_id = args.features
            output_file = args.output or f"{genome_id}_genome.json"

            print(f"Mode: Load from local features directory")
            print(f"Genome ID: {genome_id}")
            print(f"Features dir: {args.features_dir}")
            print(f"Genomes dir: {args.genomes_dir}")
            print(f"Metadata dir: {args.metadata_dir}")
            print(f"Output file: {output_file}")
            print()

            local_conv = LocalGenomeConverter()
            genome = local_conv.load_genome_from_features_dir(
                genome_id=genome_id,
                features_dir=args.features_dir,
                genomes_dir=args.genomes_dir,
                metadata_dir=args.metadata_dir,
                taxonomy=args.taxonomy,
                scientific_name=args.scientific_name
            )
            local_conv.save_genome(genome, output_file)

            if args.fasta:
                local_conv.create_fasta_from_genome(genome, args.fasta)

        # Mode 4: Create synthetic genome
        elif args.synthetic:
            if not args.genomes and not args.genome_ids:
                parser.error("--synthetic mode requires either --genomes or --genome-ids argument")

            if args.genomes and args.genome_ids:
                parser.error("--synthetic mode: provide either --genomes or --genome-ids, not both")

            asv_id = args.synthetic
            output_file = args.output or f"{asv_id}_synthetic.json"

            # Parse input
            genome_files = None
            genome_ids = None
            if args.genomes:
                genome_files = [g.strip() for g in args.genomes.split(',')]
            elif args.genome_ids:
                genome_ids = [g.strip() for g in args.genome_ids.split(',')]

            print(f"Mode: Create synthetic genome")
            print(f"ASV ID: {asv_id}")
            if genome_files:
                print(f"Source genome files: {len(genome_files)}")
            elif genome_ids:
                print(f"Source genome IDs: {len(genome_ids)}")
                print(f"Features dir: {args.features_dir}")
                print(f"Genomes dir: {args.genomes_dir}")
                print(f"Metadata dir: {args.metadata_dir}")
            print(f"Output file: {output_file}")
            if not args.no_taxonomy:
                print(f"Taxonomy directory: {args.taxonomy_dir}")
            print()

            local_conv = LocalGenomeConverter()
            genome = local_conv.create_synthetic_genome(
                asv_id=asv_id,
                genome_files=genome_files,
                genome_ids=genome_ids,
                features_dir=args.features_dir,
                genomes_dir=args.genomes_dir,
                metadata_dir=args.metadata_dir,
                taxonomy=args.taxonomy,
                template_file=args.template,
                save_taxonomy=not args.no_taxonomy,
                taxonomy_output_dir=args.taxonomy_dir
            )
            local_conv.save_genome(genome, output_file)

            if args.fasta:
                local_conv.create_fasta_from_genome(genome, args.fasta)

        print("\n" + "=" * 50)
        print("Conversion completed successfully!")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
