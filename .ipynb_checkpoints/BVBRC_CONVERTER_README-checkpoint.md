# KBase Genome Converter

A Python script that creates KBase Genome objects from multiple sources:
1. Fetches genome data from the BV-BRC (formerly PATRIC) API
2. Loads and validates local genome JSON files
3. Creates synthetic genomes by merging multiple source genomes

## Features

### BV-BRC API Mode
- Fetches genome metadata, contigs, and features from BV-BRC API
- Converts to KBase Genome object format
- Handles pagination for large feature sets
- Retrieves DNA and protein sequences for all features
- Captures ontology terms (GO, FIGFAM, PGFAM, PLFAM, RefSeq)
- Creates CDS features automatically for protein-coding genes

### Local Genome Mode
- Loads genome from local JSON files
- Validates and ensures KBase format compliance
- Fills in missing required fields with sensible defaults

### Synthetic Genome Mode
- Merges multiple source genomes into a single synthetic genome
- Collects unique functions across all source genomes
- Calculates feature probabilities (frequency across source genomes)
- Creates one contig per unique function
- Computes average GC content
- Generates genome MD5 from sorted protein sequences
- Adapted from notebook code for ASV genome creation

### General Features
- Saves genome as JSON file locally
- Optional FASTA file generation with feature sequences
- Does not require KBase authentication or upload
- Backward compatible with original command-line interface

## Requirements

```bash
pip install requests numpy
```

Note: `numpy` is only required for synthetic genome mode (calculating average GC content).

## Installation

### For Command-Line Use
Simply download the script and run it directly.

### For Notebook/Library Use
```python
# Option 1: Add to Python path
import sys
sys.path.append('/path/to/GenomeImporter')
from bvbrc_to_kbase_genome import *

# Option 2: Run in notebook
%run bvbrc_to_kbase_genome.py

# Then use classes and functions directly
genome = load_genome_from_features('511145.183')
```

See [NOTEBOOK_USAGE.md](NOTEBOOK_USAGE.md) for detailed notebook examples.

## Usage

The script supports four modes of operation:

### Mode 1: Fetch from BV-BRC API

```bash
# Basic usage
python bvbrc_to_kbase_genome.py --api <genome_id>

# With custom output file
python bvbrc_to_kbase_genome.py --api <genome_id> --output my_genome.json

# Also generate FASTA file
python bvbrc_to_kbase_genome.py --api <genome_id> --fasta features.fasta
```

**Examples:**
```bash
# E. coli K-12 MG1655
python bvbrc_to_kbase_genome.py --api 511145.183

# Mycobacterium tuberculosis H37Rv with FASTA output
python bvbrc_to_kbase_genome.py --api 83332.133 --output mtb.json --fasta mtb.fasta
```

### Mode 2: Load from Local Genome File

```bash
# Load and validate existing genome JSON
python bvbrc_to_kbase_genome.py --local genome.json

# With custom output
python bvbrc_to_kbase_genome.py --local genome.json --output validated_genome.json
```

**Examples:**
```bash
# Validate a genome file
python bvbrc_to_kbase_genome.py --local my_genome.json

# Validate and generate FASTA
python bvbrc_to_kbase_genome.py --local my_genome.json --fasta features.fasta
```

### Mode 3: Load from Features Directory (BV-BRC Local Files)

Load genomes from local BV-BRC feature files without re-downloading from the API.

**File Structure Expected:**
```
features/
  └── {genome_id}.json    # Feature metadata from BV-BRC API
genomes/
  └── {genome_id}.fna     # Genome sequences in FASTA format
```

**Usage:**
```bash
# Basic usage (assumes features/ and genomes/ in current directory)
python bvbrc_to_kbase_genome.py --features <genome_id>

# With custom directories
python bvbrc_to_kbase_genome.py --features <genome_id> \
  --features-dir path/to/features \
  --genomes-dir path/to/genomes

# With taxonomy and scientific name
python bvbrc_to_kbase_genome.py --features 511145.183 \
  --taxonomy "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia" \
  --scientific-name "Escherichia coli K-12 MG1655" \
  --output ecoli_genome.json
```

**Examples:**
```bash
# Load genome from local features directory
python bvbrc_to_kbase_genome.py --features 511145.183

# With all options
python bvbrc_to_kbase_genome.py --features 1110693.3 \
  --features-dir /path/to/features \
  --genomes-dir /path/to/genomes \
  --taxonomy "Bacteria; Firmicutes" \
  --scientific-name "Bacillus subtilis" \
  --output bacillus.json \
  --fasta bacillus.fasta
```

**What's Loaded:**
- **Features**: JSON array from `features/{genome_id}.json` containing:
  - `patric_id` - PATRIC feature identifier
  - `product` - Gene product description
  - `feature_type` - CDS, rRNA, tRNA, etc.
  - `pgfam_id`, `plfam_id`, `figfam_id` - Family identifiers
  - `annotation` - Annotation source
- **Sequences**: FASTA from `genomes/{genome_id}.fna` containing genome contigs
- **Calculated**: GC content, genome MD5, DNA size, contig counts

**Note**: The features JSON files from BV-BRC API (downloaded with limited `select` fields) don't include sequence data or precise location coordinates. The script creates features with metadata only. For complete feature sequences, use Mode 1 (API) which fetches full data.

### Mode 4: Create Synthetic Genome (from Multiple Genomes)

```bash
# Create synthetic genome from multiple source genomes
python bvbrc_to_kbase_genome.py --synthetic <asv_id> --genomes genome1.json,genome2.json,genome3.json

# With taxonomy and template
python bvbrc_to_kbase_genome.py --synthetic ASV_001 \
  --genomes g1.json,g2.json,g3.json \
  --taxonomy "Bacteria; Firmicutes; Bacilli" \
  --template template_genome.json \
  --output asv_001.json
```

**Examples:**
```bash
# Create synthetic genome from 3 source genomes
python bvbrc_to_kbase_genome.py --synthetic ASV_12345 \
  --genomes genome1.json,genome2.json,genome3.json

# With all options
python bvbrc_to_kbase_genome.py --synthetic genus.species \
  --genomes /path/to/g*.json \
  --taxonomy "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales" \
  --output synthetic_genome.json \
  --fasta synthetic_features.fasta
```

### Legacy Mode (Backward Compatible)

The original command-line interface is still supported:

```bash
# Old style (automatically uses BV-BRC API)
python bvbrc_to_kbase_genome.py <genome_id>
python bvbrc_to_kbase_genome.py <genome_id> <output_file>
```

## Output Format

The script generates a KBase Genome object with the following structure:

```json
{
  "id": "genome_id",
  "scientific_name": "Organism name",
  "domain": "Bacteria|Archaea|Eukaryota",
  "taxonomy": "taxonomic lineage",
  "genetic_code": 11,
  "dna_size": 12345678,
  "num_contigs": 10,
  "contig_ids": ["contig1", "contig2", ...],
  "contig_lengths": [1000, 2000, ...],
  "gc_content": 0.52,
  "md5": "genome_md5_hash",
  "molecule_type": "DNA",
  "source": "PATRIC",
  "source_id": "original_genome_id",
  "assembly_ref": "",
  "features": [...],
  "non_coding_features": [...],
  "cdss": [...],
  "mrnas": [],
  "feature_counts": {...},
  "publications": [],
  "genome_tiers": ["ExternalDB", "User"],
  "warnings": [],
  "taxon_ref": ""
}
```

### Feature Structure

Each feature contains:

```json
{
  "id": "genome_id_0",
  "type": "CDS|gene|tRNA|rRNA|...",
  "location": [["contig_id", start, strand, length]],
  "functions": ["functional annotation"],
  "aliases": [["source", "value"]],
  "dna_sequence": "ATCG...",
  "dna_sequence_length": 1000,
  "md5": "sequence_md5",
  "protein_translation": "MKTL...",
  "protein_translation_length": 333,
  "protein_md5": "protein_md5"
}
```

## What's Included

### BV-BRC API Mode
1. **Genome Metadata**: Name, taxonomy, GC content, genetic code
2. **Contigs**: All genome sequences with IDs and lengths
3. **Features**: All genes, CDS, RNA features with:
   - DNA and protein sequences
   - Functional annotations
   - PATRIC and RefSeq identifiers
   - Location information (contig, start, strand, length)
   - Ontology terms (GO, FIGFAM, PGFAM, PLFAM)
4. **CDS Features**: Automatically generated for all protein-coding genes
5. **Feature Statistics**: Counts by feature type

### Synthetic Genome Mode
1. **Merged Features**: Unique functions collected from all source genomes
2. **Feature Probabilities**: Each function's frequency across source genomes
3. **Simplified Assembly**: One contig per unique function
4. **Aggregated Metadata**: Average GC content, combined DNA size
5. **Provenance**: Source genome IDs stored in `source_id` field
6. **Genome MD5**: Calculated from sorted protein MD5 hashes

## What's NOT Included

- `assembly_ref`: Left blank (set to empty string) as per requirements
- `taxon_ref`: Left blank (would require KBase taxonomy lookup)
- No upload to KBase workspace
- No Assembly object creation

## Synthetic Genome Workflow

The synthetic genome mode (adapted from the BV_BRC-Copy1.ipynb notebook) works as follows:

1. **Load Source Genomes**: Reads multiple genome JSON files
2. **Collect Unique Functions**:
   - Iterates through all features in all source genomes
   - For each unique function annotation, creates one feature in the synthetic genome
   - Tracks how many source genomes contain each function
3. **Calculate Probabilities**:
   - Each function's probability = (# genomes with function) / (total genomes)
   - Stored for downstream analysis (not in the genome object itself)
4. **Create Simplified Assembly**:
   - One contig per unique feature
   - Contig ID: `{asv_id}_{index}.contig`
   - Feature location: entire contig (start=1, length=DNA sequence length)
5. **Aggregate Metadata**:
   - Average GC content across all source genomes
   - Sum of all feature DNA lengths for total DNA size
   - Taxonomy from user input or first source genome
6. **Generate MD5**:
   - Collect protein MD5 from each feature
   - Sort MD5 list
   - Calculate MD5 hash of sorted, semicolon-joined MD5 strings

This approach creates a consensus genome representing the union of functions found across multiple related genomes, useful for:
- ASV (Amplicon Sequence Variant) genome reconstruction
- Pangenome or core genome representation
- Genus-level representative genomes
- Synthetic communities

## Finding BV-BRC Genome IDs

1. Visit https://www.bv-brc.org/
2. Search for your organism
3. Click on the genome
4. The genome ID is shown in the overview (format: `XXXXX.YYY`)

## API Endpoints Used

The script queries these BV-BRC API endpoints:

- `/genome/` - Genome metadata
- `/genome_sequence/` - Contig sequences
- `/genome_feature/` - Gene/feature annotations (paginated)
- `/feature_sequence/` - DNA and protein sequences (batched)

## Notes

- SSL verification is disabled (following the Perl implementation pattern)
- Features are retrieved in batches of 10,000
- Sequences are fetched in batches of 100 by MD5 hash
- Large genomes may take several minutes to process
- The output JSON file can be quite large (tens of MB for typical bacterial genomes)

## Troubleshooting

### Genome ID not found
- Verify the genome ID exists at https://www.bv-brc.org/
- Ensure the ID format is correct (e.g., `511145.183`)

### SSL/Certificate errors
- The script disables SSL verification by default
- If you encounter issues, check your network connection

### Memory issues with large genomes
- The script loads all data into memory
- For very large genomes, you may need to increase available RAM

## Integration with KBase

This script generates JSON that is compatible with KBase Genome objects but does NOT:
- Upload to KBase workspace
- Create Assembly objects
- Validate against KBase type specs
- Require KBase authentication

To upload the resulting JSON to KBase, you would need to:
1. Create an Assembly object separately
2. Update the `assembly_ref` field
3. Use KBase workspace client to save the object
4. Validate against the KBase.Genome type specification

## License

This script is part of the GenomeImporter repository and follows the same license.
