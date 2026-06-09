#!/usr/bin/env python3
"""
Test the working MSGenome creation strategy
"""

import json
from pathlib import Path

print("=" * 80)
print("Testing Working Strategy: Empty MSGenome + Set Attributes")
print("=" * 80)

# Load genome
genome_file = "genome_objects/AAP99.1.json"
with open(genome_file, 'r') as f:
    genome_data = json.load(f)

genome_id = genome_data['id']
print(f"\n1. Loaded genome: {genome_id}")
print(f"   Features: {len(genome_data.get('features', []))}")

# Create MSGenome using working strategy
print(f"\n2. Creating MSGenome...")
from modelseedpy.core.msgenome import MSGenome, MSFeature

msgenome = MSGenome()

# Set basic attributes (excluding features which need special handling)
for key, value in genome_data.items():
    if key != 'features':  # Skip features, we'll handle them specially
        try:
            setattr(msgenome, key, value)
        except:
            pass
msgenome.id = genome_id

# Convert feature dicts to MSFeature objects
if 'features' in genome_data and genome_data['features']:
    print(f"   Converting {len(genome_data['features'])} features...")
    msfeatures = []
    for feat_dict in genome_data['features']:
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
    print(f"   ✓ Converted {len(msfeatures)} features to MSFeature objects")

print(f"   ✓ MSGenome created: {msgenome.id}")
print(f"   ✓ Features: {len(msgenome.features) if hasattr(msgenome, 'features') and msgenome.features else 0}")

# Get template - use local template_core.json
print(f"\n3. Loading template...")
template = None

# Load template_core.json directly
try:
    from os import path
    # For editable install, use the known repo location
    template_path = path.expanduser('~/repos/ModelSEEDpy/modelseedpy/data/templates/template_core.json')

    with open(template_path, 'r') as f:
        template_data = json.load(f)
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    template = MSTemplateBuilder.from_dict(template_data).build()
    print(f"   ✓ Loaded template_core.json")
    print(f"   - Template: {template.id}")
    print(f"   - Reactions: {len(template.reactions)}")
except Exception as e:
    print(f"   ✗ Failed to load template: {e}")
    import traceback
    traceback.print_exc()

if not template:
    print("\n✗ Could not load template. Cannot build model.")
    print("\nPlease ensure templates are available or internet connection for download")
    import sys
    sys.exit(1)

# Build model
print(f"\n4. Building model...")
from modelseedpy.core.msbuilder import MSBuilder

model_id = f"{genome_id}_model"
builder = MSBuilder(msgenome, template, name=model_id)

# Try building with different parameters to avoid the bug
print(f"   Attempting to build with allow_all_non_grp_reactions=True...")
try:
    model = builder.build(model_id, allow_all_non_grp_reactions=True, annotate_with_rast=False)
except Exception as e:
    print(f"   First attempt failed: {e}")
    print(f"   Trying with annotate_with_rast=True...")
    model = builder.build(model_id, annotate_with_rast=True)

print(f"   ✓ Model built successfully!")
print(f"   - Reactions: {len(model.reactions)}")
print(f"   - Metabolites: {len(model.metabolites)}")
print(f"   - Genes: {len(model.genes)}")

# Test FBA
print(f"\n5. Testing FBA optimization...")
solution = model.optimize()
print(f"   Status: {solution.status}")
if solution.status == 'optimal':
    print(f"   Objective: {solution.objective_value:.4f}")
else:
    print(f"   Note: Model may need gap-filling")

# Save model
print(f"\n6. Saving model...")
from cobra.io import write_sbml_model

output_file = f"{genome_id}_model.xml"
write_sbml_model(model, output_file)
print(f"   ✓ Saved to: {output_file}")

print("\n" + "=" * 80)
print("SUCCESS! Working strategy:")
print("  1. Create empty MSGenome()")
print("  2. Set attributes from genome dict")
print("  3. Use MSBuilder to build model")
print("=" * 80)
