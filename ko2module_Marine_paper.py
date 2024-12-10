import os
import csv

def read_ko_module_mapping(mapping_file_path):
    """
    Reads the KO-to-Module mapping from the provided file path and returns a dictionary.
    """
    ko_module_mapping = {}
    with open(mapping_file_path, "r", encoding='utf-8') as mapping_file:
        next(mapping_file)  # Skip the header
        for line in mapping_file:
            ko, module = line.rstrip().split("\t")
            ko_number = ko.split(":")[1]
            module_number = module.split(":")[1]
            ko_module_mapping[ko_number] = module_number
    return ko_module_mapping

def process_file(input_path, output_path, ko_module_mapping, ko_column_index=259):
    """
    Processes a single input TSV file, adds the corresponding Module column, 
    and removes the KO column before saving the result.
    """
    with open(input_path, "r", encoding='utf-8') as input_file, \
            open(output_path, "w", newline="", encoding='utf-8') as output_file:
        
        tsv_reader = csv.reader(input_file, delimiter='\t')
        tsv_writer = csv.writer(output_file, delimiter='\t')

        # Process the header
        header = next(tsv_reader)
        header.append("Module")  # Add "Module" column
        del header[ko_column_index]  # Remove the KO column
        tsv_writer.writerow(header)

        # Process each data row
        for row in tsv_reader:
            ko_number = row[ko_column_index]

            # Find the module corresponding to the KO number or set to "NA"
            module = "NA" if ko_number == "NA" else ko_module_mapping.get(ko_number, "NA")
            row.append(module)  # Add module column
            del row[ko_column_index]  # Remove KO column

            tsv_writer.writerow(row)

def process_directory(input_dir, output_dir, ko_module_mapping):
    """
    Processes all TSV files in the input directory, adds Module column, 
    removes KO column, and writes the output to the output directory.
    """
    for filename in os.listdir(input_dir):
        if filename.endswith(".tsv"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            process_file(input_path, output_path, ko_module_mapping)

def main():
    # File paths
    input_dir = '/home/projects/zeevid/Analyses/2023-Tamir/Geotraces/kmer_phylum_class_ko_trimmed/'
    output_dir = '/home/projects/zeevid/Analyses/2023-Tamir/Geotraces/kmer_phylum_class_ko_module_trimmed/'
    mapping_file_path = '/home/projects/zeevid/Data/Databases/KEGG/Processed_Files/ko_module.tsv'

    # Read KO-to-Module mapping
    ko_module_mapping = read_ko_module_mapping(mapping_file_path)
    process_directory(input_dir, output_dir, ko_module_mapping)

if __name__ == "__main__":
    main()
