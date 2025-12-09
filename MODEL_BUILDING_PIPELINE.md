# Model Building Pipeline Documentation

## Overview

This document describes the complete pipeline for building genome-scale metabolic models using `GenomeLoader().build_model()`. The pipeline integrates KBase genome objects with ModelSEEDpy's model reconstruction capabilities.

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         GenomeLoader().build_model()                        │
│                     (load_genomes_to_modelseedpy.py)                        │
└───────────────────────────────────┬─────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           1. Load Genome JSON                               │
│                         (KBase Genome Object)                               │
└───────────────────────────────────┬─────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      2. Create MSGenome Object                              │
│                   (ModelSEEDpy genome representation)                       │
└───────────────────────────────────┬─────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         3. Load Template                                    │
│              (template_core.json from ModelSEEDpy)                          │
└───────────────────────────────────┬─────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    4. MSBuilder Model Construction                          │
│                     (ModelSEEDpy core builder)                              │
└───────────────────────────────────┬─────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       5. COBRApy Model Output                               │
│                    (Ready for FBA simulations)                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Layer 1: GenomeLoader Operations

**File:** `load_genomes_to_modelseedpy.py`

### GenomeLoader Class

The `GenomeLoader` class provides the high-level interface for loading genomes and building models.

#### Key Methods

| Method | Description |
|--------|-------------|
| `__init__(genome_dir)` | Initialize with directory containing genome JSON files (default: `genome_objects/`) |
| `list_available_genomes()` | Returns list of genome IDs available in the directory |
| `load_genome(genome_id)` | Loads KBase genome object from JSON file |
| `get_msgenome(genome_id)` | Converts KBase genome to ModelSEEDpy MSGenome object |
| `build_model(genome_id, template, model_id)` | Builds complete metabolic model |

### build_model() Function Flow

```python
def build_model(self, genome_id, template=None, model_id=None, printing=True):
    # Step 1: Get MSGenome object (creates if not cached)
    msgenome = self.get_msgenome(genome_id, printing=printing)

    # Step 2: Load template (hardcoded to template_core.json)
    template_path = '~/repos/ModelSEEDpy/modelseedpy/data/templates/template_core.json'
    template_obj = MSTemplateBuilder.from_dict(template_data).build()

    # Step 3: Create MSBuilder and build model
    builder = MSBuilder(msgenome, template_obj, name=model_id)
    model = builder.build(model_id, allow_all_non_grp_reactions=True, annotate_with_rast=False)

    return model  # COBRApy model
```

### get_msgenome() Function Flow

```python
def get_msgenome(self, genome_id, printing=True):
    # Step 1: Load raw genome JSON
    genome = self.load_genome(genome_id, printing=printing)

    # Step 2: Create empty MSGenome
    msgenome = MSGenome()

    # Step 3: Copy genome attributes (excluding features)
    for key, value in genome.items():
        if key != 'features':
            setattr(msgenome, key, value)

    # Step 4: Convert feature dicts to MSFeature objects
    for feat_dict in genome['features']:
        msfeature = MSFeature(
            feature_id=feat_dict.get('id', ''),
            sequence=feat_dict.get('dna_sequence', feat_dict.get('protein_translation', '')),
            description=feat_dict.get('function', feat_dict.get('type', '')),
            aliases=feat_dict.get('aliases', [])
        )
        # Add ontology terms if present
        if 'ontology_terms' in feat_dict:
            for term_type, terms in feat_dict['ontology_terms'].items():
                msfeature.add_ontology_term(term_type, term)
        msfeatures.append(msfeature)

    msgenome.features = msfeatures
    return msgenome
```

---

## Layer 2: ModelSEEDpy Core Operations

### MSGenome (modelseedpy/core/msgenome.py)

The `MSGenome` class represents a genome in ModelSEEDpy format.

#### MSGenome Structure

```python
class MSGenome:
    def __init__(self):
        self.features = DictList()  # List of MSFeature objects
        self.id = None
        self.annoont = None
        self.scientific_name = None
```

#### MSFeature Structure

```python
class MSFeature:
    def __init__(self, feature_id, sequence, description=None, aliases=None):
        self.id = feature_id
        self.seq = sequence
        self.description = description
        self.ontology_terms = {}  # Dict mapping ontology type -> list of terms
        self.aliases = aliases or []

    def add_ontology_term(self, ontology_term, value):
        # Add functional annotation (e.g., RAST, EC numbers)
```

### MSTemplate (modelseedpy/core/mstemplate.py)

The template defines the "universe" of possible reactions, metabolites, and gene-protein-reaction (GPR) associations.

#### Template Components

| Component | Description |
|-----------|-------------|
| `compartments` | Cellular compartments (cytosol, extracellular, etc.) |
| `compounds` | Metabolite definitions |
| `compcompounds` | Compartmentalized metabolites (species) |
| `roles` | Functional roles (enzyme functions) |
| `complexes` | Protein complexes that catalyze reactions |
| `reactions` | Reaction definitions with stoichiometry |
| `biomasses` | Biomass composition definitions |
| `drains` | Sink/demand reactions |

#### Template Reaction Types

```python
class TemplateReactionType(Enum):
    CONDITIONAL = "conditional"   # Requires gene evidence
    UNIVERSAL = "universal"       # Always included
    SPONTANEOUS = "spontaneous"   # Non-enzymatic reactions
    GAPFILLING = "gapfilling"     # For gap-filling only
```

### MSBuilder (modelseedpy/core/msbuilder.py)

The `MSBuilder` class is the core model construction engine.

#### MSBuilder Initialization

```python
class MSBuilder:
    def __init__(self, genome, template=None, name=None, ontology_term="RAST", index="0"):
        self.name = name
        self.genome = genome
        self.template = template
        self.genome_class = None

        # Build search index: normalized role name -> gene IDs
        self.search_name_to_genes, self.search_name_to_original = _aSearch(
            genome, ontology_term
        )
```

#### Role Normalization

The `normalize_role()` function standardizes functional role names for matching:

```python
def normalize_role(s):
    s = s.strip().lower()
    s = re.sub(r"[\W_]+", "", s)  # Remove non-alphanumeric characters
    return s
```

#### MSBuilder.build() Method

The main model construction method:

```python
def build(self, model_or_id, index="0", allow_all_non_grp_reactions=False,
          annotate_with_rast=True, biomass_classic=False, biomass_gc=0.5,
          add_reaction_from_rast_annotation=True):

    # 1. Optionally annotate genome with RAST
    if annotate_with_rast:
        rast = RastClient()
        res = rast.annotate_genome(self.genome)
        self.search_name_to_genes, self.search_name_to_original = _aSearch(
            self.genome, "RAST"
        )

    # 2. Auto-select template if not provided
    if not self.template:
        self.auto_select_template()

    # 3. Create base COBRA model
    cobra_model = Model(model_or_id)
    self.base_model = cobra_model

    # 4. Generate reaction-to-complex mappings (GPR rules)
    self.generate_reaction_complex_sets()
    complex_groups = self.build_complex_groups(...)

    # 5. Add metabolic reactions (gene-associated)
    if add_reaction_from_rast_annotation:
        metabolic_reactions = self.build_metabolic_reactions()
        cobra_model.add_reactions(metabolic_reactions)

    # 6. Add non-metabolite reactions (universal/spontaneous)
    non_metabolic_reactions = self.build_non_metabolite_reactions(
        cobra_model, allow_all_non_grp_reactions
    )
    cobra_model.add_reactions(non_metabolic_reactions)

    # 7. Add exchange reactions
    self.add_exchanges_to_model(cobra_model)

    # 8. Build and add biomass reactions
    for rxn_biomass in self.template.biomasses:
        reaction = rxn_biomass.build_biomass(cobra_model, index, ...)
        biomass_reactions.append(reaction)
    cobra_model.objective = biomass_reactions[0].id

    # 9. Add sink/demand reactions
    reactions_sinks = self.build_drains()
    cobra_model.add_reactions(reactions_sinks)

    return cobra_model
```

### Gene-Protein-Reaction (GPR) Mapping

The GPR mapping process:

1. **Extract functional roles from genome features** via ontology terms
2. **Normalize role names** using `normalize_role()`
3. **Match roles to template complexes** - each complex has associated roles
4. **Build GPR rules**:
   - Roles within a complex are ANDed (all required)
   - Multiple complexes for a reaction are ORed (any sufficient)

```python
def get_gpr_from_template_reaction(self, template_reaction, allow_incomplete_complexes=True):
    template_reaction_complexes = self._get_template_reaction_complexes(template_reaction)
    gpr_set = self._build_reaction_complex_gpr_sets(
        template_reaction_complexes, allow_incomplete_complexes
    )
    return gpr_set
```

#### GPR String Construction

```python
def build_gpr2(cpx_sets):
    list_of_ors = []
    for cpx in cpx_sets:
        list_of_ands = []
        for role_id in cpx_sets[cpx]:
            gene_ors = cpx_sets[cpx][role_id]
            if len(gene_ors) > 1:
                list_of_ands.append("(" + " or ".join(gene_ors) + ")")
            else:
                list_of_ands.append(list(gene_ors)[0])
        list_of_ors.append("(" + " and ".join(list_of_ands) + ")")
    if len(list_of_ors) > 1:
        return " or ".join(list_of_ors)
    return list_of_ors[0]
```

### Biomass Reaction Construction

The `MSTemplateBiomass.build_biomass()` method constructs the biomass reaction:

#### Biomass Component Types

| Type | Description |
|------|-------------|
| `dna` | DNA components |
| `rna` | RNA components |
| `protein` | Amino acids |
| `lipid` | Lipid components |
| `cellwall` | Cell wall components |
| `cofactor` | Cofactors |
| `energy` | ATP requirements |
| `other` | Other components |

#### Coefficient Types

```python
class TemplateBiomassCoefficientType(Enum):
    MOLFRACTION = "MOLFRACTION"  # Mole fraction based on abundance
    MOLSPLIT = "MOLSPLIT"        # Equal split among class
    MULTIPLIER = "MULTIPLIER"    # Multiply by abundance
    EXACT = "EXACT"              # Use exact coefficient
```

#### GC Content Adjustment

For DNA components, coefficients are adjusted based on GC content:
- `AT` type: `2 * coefficient * (1 - GC) * (abundance / total_mw)`
- `GC` type: `2 * coefficient * GC * (abundance / total_mw)`

---

## Layer 3: Template Selection

### Current Implementation Assumptions

**IMPORTANT:** The current `GenomeLoader.build_model()` implementation has a **hardcoded template path**:

```python
template_path = path.expanduser('~/repos/ModelSEEDpy/modelseedpy/data/templates/template_core.json')
```

This means:
- **All genomes use the same core template** regardless of taxonomy
- No automatic template selection based on Gram stain or domain
- The template parameter in `build_model()` is currently ignored

### MSBuilder.auto_select_template() (Not Currently Used)

If template auto-selection were enabled, it would:

1. **Load genome classifier**: KNN classifier trained on RAST annotations
2. **Classify genome** into one of four classes:
   - `A` = Archaea
   - `C` = Cyanobacteria
   - `N` = Gram Negative
   - `P` = Gram Positive

3. **Map class to template**:

```python
template_genome_scale_map = {
    "A": "template_gram_neg",
    "C": "template_gram_neg",
    "N": "template_gram_neg",
    "P": "template_gram_pos",
}
```

### Genome Classification

```python
class MSGenomeClass(Enum):
    P = "Gram Positive"
    N = "Gram Negative"
    C = "Cyano"
    A = "Archaea"
```

The classifier (`knn_ACNP_RAST_filter_01_17_2023`) uses:
- **Input**: RAST functional annotations from the genome
- **Method**: K-Nearest Neighbors classification
- **Output**: Single character class code (A, C, N, or P)

---

## Key Assumptions and Defaults

### 1. Template Selection

| Assumption | Value | Impact |
|------------|-------|--------|
| Template file | `template_core.json` | Same template for all organisms |
| Auto-selection | Disabled | Gram-specific templates not used |
| Archaea handling | Uses gram-negative template | May be suboptimal for archaeal metabolisms |

### 2. Model Building Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `allow_all_non_grp_reactions` | `True` | Includes universal/spontaneous reactions even without metabolite context |
| `annotate_with_rast` | `False` | Skips RAST re-annotation (assumes annotations already present) |
| `index` | `"0"` | Compartment index suffix |
| `ontology_term` | `"RAST"` | Ontology type to search for gene functions |

### 3. Reaction Inclusion Rules

| Reaction Type | Inclusion Rule |
|---------------|----------------|
| Conditional | Only if gene evidence exists via GPR matching |
| Universal | Always included (or if `allow_all_non_grp_reactions=True`) |
| Spontaneous | Always included (non-enzymatic) |

### 4. Biomass Construction

| Parameter | Default | Description |
|-----------|---------|-------------|
| `classic` | `False` | Uses modern biomass with DNA/RNA/Protein synthesis reactions |
| `GC` | `0.5` | GC content (50%) - not extracted from genome |
| `add_to_model` | `True` | Automatically adds biomass to model |

### 5. Exchange and Drain Reactions

- **Exchanges**: Created for all metabolites in extracellular compartment (`e0`)
- **Sinks**: Created for specific metabolites (e.g., biomass, glycogen)
- **Demands**: Created based on template drain specifications

**Default sinks** (if template doesn't specify drains):
```python
DEFAULT_SINKS = {
    "cpd02701_c": 1000,  # S-Adenosyl-4-methylthio-2-oxobutanoate
    "cpd11416_c": 1000,  # Biomass
    "cpd15302_c": 1000,  # glycogen(n-1)
    "cpd03091_c": 1000,  # 5'-Deoxyadenosine
    "cpd01042_c": 1000,  # p-Cresol
}
```

---

## Data Flow Summary

```
KBase Genome JSON
        │
        ▼
┌───────────────────┐
│  genome['features']│ ─────► MSFeature objects
│  genome['id']      │        with ontology_terms
│  genome['taxonomy']│
└───────────────────┘
        │
        ▼
┌───────────────────┐
│    MSGenome       │
│  .features        │ ─────► DictList of MSFeature
│  .id              │
│  .scientific_name │
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   MSBuilder       │
│  .genome          │
│  .template        │
│  .search_name_to_genes ──► {normalized_role: {gene_ids}}
└───────────────────┘
        │
        │  Match genome roles to template complexes
        │  Build GPR rules
        │  Add reactions with stoichiometry
        ▼
┌───────────────────┐
│  COBRA Model      │
│  .reactions       │ ─────► Reactions with GPRs, bounds
│  .metabolites     │ ─────► Metabolites with compartments
│  .genes           │ ─────► Genes from GPR rules
│  .objective       │ ─────► Biomass reaction
└───────────────────┘
```

---

## Potential Improvements

1. **Enable template auto-selection** based on genome classification
2. **Extract GC content** from genome object for accurate biomass composition
3. **Support multiple ontology sources** (EC, KO, etc.) for GPR matching
4. **Gap-filling integration** to improve model completeness
5. **Organism-specific biomass** based on experimental data

---

## Usage Example

```python
from load_genomes_to_modelseedpy import GenomeLoader

# Initialize loader
loader = GenomeLoader(genome_dir="genome_objects")

# List available genomes
genomes = loader.list_available_genomes()

# Build a model
model = loader.build_model("Acetobacterium.1", printing=True)

# Model statistics
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes: {len(model.genes)}")

# Run FBA
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")
```

---

## References

- **ModelSEEDpy**: https://github.com/ModelSEED/ModelSEEDpy
- **COBRApy**: https://github.com/opencobra/cobrapy
- **KBase**: https://www.kbase.us/
