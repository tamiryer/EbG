import os
import pandas as pd

def get_input_files(input_dir):
    """
    Get a list of all input files in the directory, sorted by prefix.
    """
    all_files = [f for f in os.listdir(input_dir) if f.endswith("_1.tsv") or f.endswith("_2.tsv")]
    all_files.sort(key=lambda x: x.split("_")[0])
    return all_files

def get_file_paths(input_dir, filename):
    """
    Get the input and output file paths for a given filename.
    """
    prefix = os.path.splitext(filename)[0].split("_")[0]
    input_file1 = os.path.join(input_dir, filename)
    input_file2 = os.path.join(input_dir, f'{prefix}_{"2" if "1" in filename else "1"}.tsv')
    output_file = os.path.join(output_dir, f'{prefix}.tsv')
    return input_file1, input_file2, output_file

def concatenate_files(input_file1, input_file2, output_file):
    """
    Concatenate two input TSV files and write the combined DataFrame to the output file.
    """
    df1 = pd.read_csv(input_file1, sep='\t')
    df2 = pd.read_csv(input_file2, sep='\t')
    combined_df = pd.concat([df1, df2], ignore_index=True)
    combined_df.to_csv(output_file, sep='\t', index=False)

def main():
    global input_dir, output_dir
    input_dir = "/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/temp"
    output_dir = "/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/kmer_counts_5k_paired"

    os.makedirs(output_dir, exist_ok=True)

    all_files = get_input_files(input_dir)

    for i, filename in enumerate(all_files):
        input_file1, input_file2, output_file = get_file_paths(input_dir, filename)

        # Check if output file already exists
        if os.path.exists(output_file):
            print(f"Skipping {os.path.splitext(filename)[0]}, output file already exists")
            continue

        # Check if both input files exist
        if not (os.path.exists(input_file1) and os.path.exists(input_file2)):
            print(f"Skipping {os.path.splitext(filename)[0]}, one or both input files do not exist")
            continue

        concatenate_files(input_file1, input_file2, output_file)
        print(f"Concatenated files for {os.path.splitext(filename)[0]} into {output_file}")

        # Alternate the order of the input files for the next iteration
        if i % 2 == 0:
            all_files[i], all_files[i+1] = all_files[i+1], all_files[i]

if __name__ == "__main__":
    main()