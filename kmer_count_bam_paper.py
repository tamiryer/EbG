import os
import sys
from collections import Counter
from itertools import product
import pysam

ALL_TETRAMERS = [''.join(tetramer) for tetramer in product('ACGT', repeat=4)] #Creates all 256 4-mers possibilities as headers

def calculate_tetramer_count(sequence):
    return Counter(sequence[i:i+4] for i in range(len(sequence) - 3) if 'N' not in sequence[i:i+4])

def process_bam(input_file, output_file):
    with pysam.AlignmentFile(input_file, "rb", require_index=False) as bamfile, open(output_file, 'w') as fout:
        # Write header
        fout.write("Read\t" + "\t".join(ALL_TETRAMERS) + "\n")
        
        # Process reads and write results in batches
        batch_size = 1000
        batch = []
        for read in bamfile:
            reference_name = bamfile.get_reference_name(read.reference_id)
            count = calculate_tetramer_count(read.query_sequence)
            batch.append((reference_name, count))
            
            if len(batch) >= batch_size:
                write_batch(fout, batch)
                batch = []
        
        # Write any remaining reads in the last batch
        if batch:
            write_batch(fout, batch)

def write_batch(fout, batch):
    for reference_name, count in batch:
        fout.write(f"{reference_name}\t" + "\t".join(str(count.get(t, 0)) for t in ALL_TETRAMERS) + "\n")

def create_done_flag(output_file):
    done_file = output_file + ".done"
    with open(done_file, 'w') as f:
        f.write("Processing completed successfully.")
    print(f"Created done flag: {done_file}")

def main(input_file, output_file):
    if os.path.exists(output_file):
        print(f"Skipping {input_file} as output already exists.")
        return
    
    try:
        process_bam(input_file, output_file)
        create_done_flag(output_file)
        print(f"Processing completed for {input_file}")
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        # Optionally, you can create an error flag file here
        with open(output_file + ".error", 'w') as f:
            f.write(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_bam_file> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
