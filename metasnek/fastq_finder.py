import os
import warnings
import glob
import csv
import re


def parse_directory(file_list):
    """Pairs samples from a list of files.

    Args:
        file_list (list): A list of file paths.

    Returns:
        tuple: A tuple containing two sets:
            - paired_files: A set of tuples with the sample name, unpaired file path, and paired file path.
            - unpaired_files: A set of tuples with the sample name and unpaired file path.
    """

    paired_files = set()
    unpaired_files = set()

    for file in file_list:
        file_path, file_ext = os.path.splitext(file)
        file_name = os.path.basename(file)
        file_ext = file_ext.lower()

        pattern = r"\.(fasta|fastq)(\.gz)?$"
        if re.search(pattern, file_name, re.IGNORECASE):
            R1_file = None
            R2_file = None

            if "_R1" in file_name:
                sample_name = file_name.rsplit("_R1", 1)[0]
                R2_file = file_path.replace("_R1", "_R2") + file_ext
                R1_file = file
            elif "_R2" in file_name:
                sample_name = file_name.rsplit("_R2", 1)[0]
                R1_file = file_path.replace("_R2", "_R1") + file_ext
                R2_file = file
            if R1_file and R2_file and R1_file in file_list and R2_file in file_list:
                paired_files.add((sample_name, R1_file, R2_file))
            else:
                sample_name = re.split(r"\.(fasta|fastq)(\.gz)?$", file_name)[0]
                if "_R1" in sample_name or "_R2" in sample_name:
                    warnings.warn(f"Orphaned paired read detected: {file_name}", Warning)
                unpaired_files.add((sample_name, file))

    return paired_files, unpaired_files


def parse_tsv_file(file_path):
    """Parses a 3-column TSV file of sample names and sequencing reads (column 3 is optional)

    Args:
        file_path (str): Path to the TSV file.

    Returns:
        tuple: A tuple containing two lists:
            - paired_reads: A list of tuples with the sample name, R1 file, and R2 file (if available).
            - unpaired_reads: A list of tuples with the sample name and R1 file (for unpaired reads).
    """

    paired_reads = set()
    unpaired_reads = set()

    with open(file_path, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')

        for row in reader:

            sample_name = row[0].strip()
            r1_file = row[1].strip()
            r2_file = row[2].strip() if len(row) >= 3 else None

            if not os.path.isfile(r1_file):
                raise FileNotFoundError(f"R1 file '{r1_file}' does not exist.")

            if r2_file and not os.path.isfile(r2_file) and not r2_file.lower() in ["none", "null"]:
                raise FileNotFoundError(f"R2 file '{r2_file}' does not exist.")

            if r2_file:
                paired_reads.add((sample_name, r1_file, r2_file))
            else:
                unpaired_reads.add((sample_name, r1_file))

    return paired_reads, unpaired_reads


def parse_samples(input_file_or_directory):
    """Work out if filepath is a file or directory and run appropriate parser

    Args:
        input_file_or_directory (str): filepath for TSV file or directory of reads

    Returns:
        tuple: A tuple containing two lists:
            - paired_reads: A list of tuples with the sample name, R1 file, and R2 file (if available).
            - unpaired_reads: A list of tuples with the sample name and R1 file (for unpaired reads).
    """
    if os.path.isdir(input_file_or_directory):
        file_list = glob.glob(os.path.join(input_file_or_directory, "*"))
        paired_files, unpaired_files = parse_directory(file_list)
    elif os.path.isfile(input_file_or_directory):
        try:
            paired_files, unpaired_files = parse_tsv_file(input_file_or_directory)
        except FileNotFoundError as e:
            raise ValueError("Parse_samples failed with error from parse_tsv_file: " + str(e))
    else:
        raise ValueError(f"{input_file_or_directory} is neither a file nor directory")

    if len(paired_files) == 0 and len(unpaired_files) == 0:
        raise ValueError(f"Failed to detect any reads files and samples for {input_file_or_directory}")

    return paired_files, unpaired_files


def convert_to_dictionary(paired_reads, unpaired_reads):
    """Converts paired and unpaired reads to a single dictionary.

    Args:
        paired_reads (set): A list of tuples with sample name, R1 file, and R2 file (if available).
        unpaired_reads (set): A list of tuples with sample name and R1 file (for unpaired reads).

    Returns:
        dict:
            - sample name (dict):
                - R1 (str): filepath of R1 reads file
                - R2 (str): filepath of R2 reads file or None for unpaired
    """

    reads_dictionary = {}

    for sample_name, r1_file, r2_file in paired_reads:
        if sample_name not in reads_dictionary:
            reads_dictionary[sample_name] = {'R1': r1_file, 'R2': r2_file}

    for sample_name, r1_file in unpaired_reads:
        if sample_name not in reads_dictionary:
            reads_dictionary[sample_name] = {'R1': r1_file, 'R2': None}

    return reads_dictionary


def parse_samples_to_dictionary(input_file_or_directory):
    """Convenience function to parse the samples directory or TSV and return the samples dictionary

    Args:
        input_file_or_directory (str): filepath of samples TSV or directory

    Returns:
        dict:
            - sample name (dict):
                - R1 (str): filepath of R1 reads file
                - R2 (str): filepath of R2 reads file or None for unpaired
    """
    paired, unpaired = parse_samples(input_file_or_directory)
    sample_dictionary = convert_to_dictionary(paired, unpaired)
    return sample_dictionary


def write_samples_tsv(dictionary, output_file):
    """Write the samples dictionary to a TSV file

    Args:
        dictionary:
            - sample name (dict):
                - R1 (str): filepath of R1 reads file
                - R2 (str): filepath of R2 reads file or None for unpaired
        output_file (str): filepath of output file for writing
    """

    with open(output_file, "w") as out:
        for sample in dictionary.keys():
            if dictionary[sample]['R2'] is not None and dictionary[sample]['R2'].lower() not in ["none", "null"]:
                out.write(f"{sample}\t{dictionary[sample]['R1']}\t{dictionary[sample]['R2']}\n")
            else:
                out.write(f"{sample}\t{dictionary[sample]['R1']}\n")
