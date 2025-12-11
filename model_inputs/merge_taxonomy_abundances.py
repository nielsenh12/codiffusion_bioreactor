import csv

# Read taxonomy file to get taxonomy info per seq
with open('taxonomy.csv', 'r') as f:
    reader = csv.DictReader(f)
    taxonomy_cols = reader.fieldnames
    taxonomy_data = {row['seq']: row for row in reader}

# Read abundances file
with open('abundances.csv', 'r') as f:
    reader = csv.DictReader(f)
    abundance_cols = reader.fieldnames
    abundances_data = {row['seq']: row for row in reader}

# Get sample columns from abundances
sample_columns = [col for col in abundance_cols if col != 'seq']

# Output columns same as taxonomy: seq, sample, rel_ab, Kingdom, ...
output_cols = taxonomy_cols

# Create expanded data: one row per seq/sample combination
expanded_data = []
for seq, tax_row in taxonomy_data.items():
    ab_row = abundances_data.get(seq, {})

    for sample in sample_columns:
        rel_ab = ab_row.get(sample, '0')

        output_row = {
            'seq': seq,
            'sample': sample,
            'rel_ab': rel_ab,
            'Kingdom': tax_row['Kingdom'],
            'Phylum': tax_row['Phylum'],
            'Class': tax_row['Class'],
            'Order': tax_row['Order'],
            'Family': tax_row['Family'],
            'Genus': tax_row['Genus'],
            'Species': tax_row['Species'],
            'date': tax_row['date'],
            'media': tax_row['media'],
            'timepoint': tax_row['timepoint'],
            'notes': tax_row['notes'],
        }
        expanded_data.append(output_row)

print(f"Original taxonomy rows: {len(taxonomy_data)}")
print(f"Samples: {len(sample_columns)}")
print(f"Expanded rows: {len(expanded_data)}")

# Write to total.csv
with open('total.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=output_cols)
    writer.writeheader()
    writer.writerows(expanded_data)

print(f"Saved to total.csv")
