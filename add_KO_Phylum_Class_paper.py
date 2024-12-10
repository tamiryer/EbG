##Marine - using polars instead of pandas - the best version
import polars as pl

# Load the input TSV file (4-mer counts)
def load_input_tsv(file_path):
    return pl.read_csv(file_path, separator='\t')

# Load the catalog TSV.GZ file
def load_catalog(file_path):
    # Load the catalog, keeping only necessary columns
    return pl.read_csv(file_path, separator='\t', has_header=True, 
                       columns=['OM-RGC_ID', 'NewPhylum', 'Class', 'KO'])

# Merge the 4-mer counts with the catalog information
def merge_data(input_df, catalog_df):
    merged_df = input_df.join(catalog_df, left_on='Read', right_on='OM-RGC_ID', how='left')
    # Drop the redundant 'gene' column after merging
    merged_df = merged_df.drop('OM-RGC_ID')
    merged_df = merged_df.with_columns([
        pl.col('NewPhylum').fill_null('NA'),
        pl.col('Class').fill_null('NA'),
        pl.col('KO').fill_null('NA')
    ])
    merged_df.rename({"NewPhylum":"Phylum"})
    return merged_df

def save_output(merged_df, output_file):
    merged_df.write_csv(output_file, separator='\t')

# Main function to tie everything together
def process_files(input_tsv, catalog_tsv_gz, output_file):
    input_df = load_input_tsv(input_tsv)            # Load input 4-mer TSV
    catalog_df = load_catalog(catalog_tsv_gz)       # Load catalog TSV.GZ
    merged_df = merge_data(input_df, catalog_df)    # Merge the two DataFrames
    save_output(merged_df, output_file)             # Save the result to a new TSV

# Example usage
input_tsv = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/GEOTRACES/full/kmer_count_csv/4mers/SRR5788468_1.tsv'            # Path to your input TSV file with 4-mers
catalog_tsv_gz = '/home/projects/zeevid/Data/Databases/TARA_OMRGC.v2/OM-RGC_v2_noseq_new_phyla.tsv.gz'  # Path to your gzipped catalog TSV file
output_file = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/GEOTRACES/full/kmer_count_csv/polars_chat_SRR5788468_1.tsv'         # Output file path
process_files(input_tsv, catalog_tsv_gz, output_file)


###Soil - using polars instead of soil - same as the script above

# import polars as pl

# # Load the input TSV file (4-mer counts)
# def load_input_tsv(file_path):
#     # Load the TSV file into a DataFrame
#     return pl.read_csv(file_path, separator='\t')

# # Load the catalog TSV.GZ file
# def load_catalog(file_path):
#     # Load the catalog, keeping only necessary columns
#     return pl.read_csv(file_path, separator='\t', has_header=True, 
#                        columns=['gene', 'phylum', 'class', 'KO'])

# # Merge the 4-mer counts with the catalog information
# def merge_data(input_df, catalog_df):
#     # Merge the input TSV with catalog using 'Read' and 'gene' columns
#     merged_df = input_df.join(catalog_df, left_on='Read', right_on='gene', how='left')
#     # Drop the redundant 'gene' column after merging
#     merged_df = merged_df.drop('gene')
#     # Fill missing 'phylum' and 'class' values with 'NA'
#     merged_df = merged_df.with_columns([
#         pl.col('phylum').fill_null('NA'),
#         pl.col('class').fill_null('NA'),
#         pl.col('KO').fill_null('NA')
#     ])
#     merged_df.rename({"phylum":"Phylum", "class":"Class"}))
#     return merged_df

# # Save the result to a new TSV file
# def save_output(merged_df, output_file):
#     merged_df.write_csv(output_file, separator='\t')

# # Main function to tie everything together
# def process_files(input_tsv, catalog_tsv_gz, output_file):
#     input_df = load_input_tsv(input_tsv)            # Load input 4-mer TSV
#     catalog_df = load_catalog(catalog_tsv_gz)       # Load catalog TSV.GZ
#     merged_df = merge_data(input_df, catalog_df)    # Merge the two DataFrames
#     save_output(merged_df, output_file)             # Save the result to a new TSV

# # Example usage
# input_tsv = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv/BMI_HWVWMBGX7_mms_R1_D-BMI_Plate88WellB6.tsv'            # Path to your input TSV file with 4-mers
# catalog_tsv_gz = '/home/projects/zeevid/Data/Databases/GMGC_v1/metadata_GMGC10.taxonomy_class_phylum.tsv.gz'  # Path to your gzipped catalog TSV file
# output_file = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv/polars_chatgpt _BMI_HWVWMBGX7_mms_R1_D-BMI_Plate88WellB6.tsv'         # Output file path

# process_files(input_tsv, catalog_tsv_gz, output_file)


#===============================================================================

# #works great! dont touch  - Ask Dudi if he prefers the pandas over the polars
# import pandas as pd

# # Load the input TSV file (4-mer counts)
# def load_input_tsv(file_path):
#     # Load the TSV file into a DataFrame
#     return pd.read_csv(file_path, sep='\t')

# # Load the catalog TSV.GZ file
# def load_catalog(file_path):
#     # Load the catalog, keeping only necessary columns
#     return pd.read_csv(file_path, sep='\t', compression='gzip', usecols=['gene', 'phylum', 'class'])

# # Merge the 4-mer counts with the catalog information
# def merge_data(input_df, catalog_df):
#     # Merge the input TSV with catalog using 'Read' and 'gene' columns
#     merged_df = input_df.merge(catalog_df, left_on='Read', right_on='gene', how='left')
#     # Drop the redundant 'gene' column after merging
#     merged_df = merged_df.drop(columns=['gene'])
#     merged_df = merged_df.rename(columns={'phylum': 'Phylum', 'class': 'Class'})
#     return merged_df

# # Save the result to a new TSV file
# def save_output(merged_df, output_file):
#     merged_df.to_csv(output_file, sep='\t', index=False)

# # Main function to tie everything together
# def process_files(input_tsv, catalog_tsv_gz, output_file):
#     input_df = load_input_tsv(input_tsv)            # Load input 4-mer TSV
#     catalog_df = load_catalog(catalog_tsv_gz)       # Load catalog TSV.GZ
#     merged_df = merge_data(input_df, catalog_df)    # Merge the two DataFrames
#     save_output(merged_df, output_file)             # Save the result to a new TSV


# # Example usage
# # input_tsv = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv/BMI_HWVWMBGX7_mms_R1_D-BMI_Plate88WellB6.tsv'            # Path to your input TSV file with 4-mers
# # catalog_tsv_gz = '/home/projects/zeevid/Data/Databases/GMGC_v1/metadata_GMGC10.taxonomy_class_phylum.tsv.gz'  # Path to your gzipped catalog TSV file
# # output_file = '/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/taxonomy_analysis_files/NEON/Soil/full/kmer_count_csv/paper_BMI_HWVWMBGX7_mms_R1_D-BMI_Plate88WellB6.tsv'         # Output file path
# # process_files(input_tsv, catalog_tsv_gz, output_file)
