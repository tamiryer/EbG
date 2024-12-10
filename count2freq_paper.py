import os
import pandas as pd
import sys

def convert_counts_to_frequencies(input_dir, output_dir):
    """Convert count data in TSV files to frequencies and save results as TSV."""
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    for file in filter(lambda f: f.endswith('.tsv'), os.listdir(input_dir)):
        df = pd.read_csv(os.path.join(input_dir, file), index_col=0)
        df.div(df.sum(axis=0), axis=1).to_csv(os.path.join(output_dir, file))  # Convert counts to frequencies and save
        print(f"Frequencies saved to {os.path.join(output_dir, file)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python convert_counts_to_frequencies.py <input_directory> <output_directory>")

    convert_counts_to_frequencies(sys.argv[1], sys.argv[2])


