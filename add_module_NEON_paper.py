import csv
import gzip
import os

# Load gene catalog data into a dictionary
def load_gene_catalog(filepath):
    gene_catalog = {}
    with gzip.open(filepath, 'rt') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if not row[0].startswith('#'):
                gene_catalog[row[0]] = row[10]  # gene_id and Module
    return gene_catalog

# Ensure output directory exists
def ensure_output_dir(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# Process a single TSV file
def process_tsv(input_tsv_path, output_tsv_path, gene_catalog):
    with open(input_tsv_path, 'r') as input_file, open(output_tsv_path, 'w', newline='') as output_file:
        reader = csv.reader(input_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')

        header = next(reader)
        header.append('Module')  # Add 'Module' column
        writer.writerow(header)

        for row in reader:
            gene_id = row[0]  # Assuming gene_id is the first column
            row.append(gene_catalog.get(gene_id, 'NA'))  # Append module or 'NA'
            writer.writerow(row)

# Process all TSV files in input directory
def process_directory(input_dir, output_dir, gene_catalog):
    ensure_output_dir(output_dir)
    
    for filename in os.listdir(input_dir):
        if filename.endswith(".tsv"):
            input_tsv_path = os.path.join(input_dir, filename)
            output_tsv_path = os.path.join(output_dir, filename)
            process_tsv(input_tsv_path, output_tsv_path, gene_catalog)
            print(f"Enriched TSV file created: {output_tsv_path}")

# Main execution
if __name__ == "__main__":
    catalog_path = '/home/projects/zeevid/Data/Databases/GMGC_v1/GMGC10.emapper2.annotations.tsv.gz'
    input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/kmer_counts'
    output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/kmer_module'

    gene_catalog = load_gene_catalog(catalog_path)
    process_directory(input_dir, output_dir, gene_catalog)
