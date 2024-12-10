import os
import pandas as pd
from collections import defaultdict

INPUT_DIR = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv'
OUTPUT_DIR = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv'
MAX_COUNT = 5000

def initialize_counts(levels):
    """Initialize the counts and read counts dictionaries for each taxonomic level."""
    return {level: defaultdict(lambda: [0] * 256) for level in levels}, {level: defaultdict(int) for level in levels}

def sum_tetramers(df, tax_level, tetramer_headers, counts, read_counts):
    """Sum tetramers up to MAX_COUNT for each taxonomic level."""
    for taxon, group in df.groupby(tax_level):
        remaining_reads = max(0, MAX_COUNT - read_counts[tax_level][taxon])
        if remaining_reads > 0:
            tetramer_sums = group[tetramer_headers].iloc[:remaining_reads].sum().values
            counts[tax_level][taxon] = [sum(x) for x in zip(counts[tax_level][taxon], tetramer_sums)]
            read_counts[tax_level][taxon] += min(remaining_reads, len(group))

def save_results(counts, read_counts, output_paths, tetramer_headers):
    """Save results to CSV for each taxonomic level."""
    for tax_level, output_path in output_paths.items():
        result_df = pd.DataFrame(counts[tax_level], index=tetramer_headers).T
        result_df['taxon_count'] = pd.Series(read_counts[tax_level])
        result_df.to_csv(output_path, sep="\t")
 
def process_file(file_path, output_paths):
    """Process a single file to aggregate tetramer counts by taxonomic level."""
    df = pd.read_csv(file_path, sep="\t", dtype={'Phylum': str, 'Class': str, 'Module': str})
    tetramer_headers = df.columns[1:257]
    df[tetramer_headers] = df[tetramer_headers].apply(pd.to_numeric, errors='coerce').astype(int)
    counts, read_counts = initialize_counts(['Phylum', 'Class', 'Module'])

    # Aggregate counts for each taxonomic level
    for tax_level in ['Phylum', 'Class', 'Module']:
        sum_tetramers(df, tax_level, tetramer_headers, counts, read_counts)
    save_results(counts, read_counts, output_paths, tetramer_headers)

def main():
    for file_name in filter(lambda f: f.endswith(".tsv"), os.listdir(INPUT_DIR)):
        output_paths = {
            'Phylum': os.path.join(OUTPUT_DIR, 'Phylum', file_name),
            'Class': os.path.join(OUTPUT_DIR, 'Class', file_name),
            'Module': os.path.join(OUTPUT_DIR, 'Module', file_name)
        }
        # Process each file
        process_file(os.path.join(INPUT_DIR, file_name), output_paths)

if __name__ == "__main__":
    main()



##optimized pandas - super fast and good

# import os
# import pandas as pd
# from collections import defaultdict

# INPUT_DIR = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv'
# OUTPUT_DIR = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv'
# MAX_COUNT = 5000

# def process_file(file_path, output_paths):
#     # Initialize counts and read counts for both phylum and class
#     counts = {key: defaultdict(lambda: [0] * 256) for key in output_paths}
#     read_counts = {key: defaultdict(int) for key in output_paths}
#     df = pd.read_csv(file_path, sep="\t", dtype={'phylum': str, 'class': str})
#     tetramer_headers = df.columns[1:257]
#     # Cast tetramer columns to integer to avoid invalid data types
#     df[tetramer_headers] = df[tetramer_headers].apply(pd.to_numeric, errors='coerce').astype(int)

#     # Iterate over taxonomy levels 'phylum' and 'class'
#     for tax_level in ['phylum', 'class']:
#         grouped_df = df.groupby(tax_level)
#         for taxon, group in grouped_df:
#             # Get the sum of tetramers up to MAX_COUNT reads
#             if read_counts[tax_level][taxon] < MAX_COUNT:
#                 num_reads = min(MAX_COUNT - read_counts[tax_level][taxon], len(group))
#                 tetramer_sums = group[tetramer_headers].iloc[:num_reads].sum().values
                
#                 counts[tax_level][taxon] = [sum(x) for x in zip(counts[tax_level][taxon], tetramer_sums)]
#                 read_counts[tax_level][taxon] += num_reads
#         # Create output DataFrame and save results
#         result_df = pd.DataFrame(counts[tax_level], index=tetramer_headers).T
#         result_df['taxon_count'] = pd.Series(read_counts[tax_level])
#         result_df.to_csv(output_paths[tax_level], sep="\t")

# def main():
#     for file_name in os.listdir(INPUT_DIR):
#         if file_name.endswith(".tsv"):
#             output_paths = {
#                 'phylum': os.path.join(OUTPUT_DIR, 'Phylum', f'opt_pd_{file_name}.tsv'),
#                 'class': os.path.join(OUTPUT_DIR, 'Class', f'opt_pd_{file_name}.tsv')
#             }

#             # Skip existing output files
#             if all(os.path.exists(path) for path in output_paths.values()):
#                 print(f"Output files already exist for {file_name}, skipping...")
#                 continue

#             process_file(os.path.join(INPUT_DIR, file_name), output_paths)

# if __name__ == "__main__":
#     main()

