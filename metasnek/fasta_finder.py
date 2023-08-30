import os
import warnings
import glob
import csv
import re


def fastas_from_directory(fasta_directory):
    """Find all the fasta files in a directory and return them as a dictionary

    Args:
        fasta_directory (str): filepath to directory

    Returns:
        fasta_files (dict):
            key (str): fasta file name
            value (str): complete file path
    """

    fasta_files = {}
    file_list = glob.glob(os.path.join(fasta_directory, "*"))
    for file_path in file_list:
        file_name = os.path.basename(file_path)
        if file_name.lower().endswith(
            (".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn")
        ):
            fasta_files[file_name] = os.path.join(fasta_directory, file_path)
    return fasta_files


def parse_tsv_file(file_path):
    """Parse a 2-column TSV of col 1: reference name and col 2: fasta-format filepath

    Args:
        file_path (str): filepath of TSV file

    Returns:
        fasta_files (dict):
            key (str): file name
            value (str): filepath
    """

    fasta_files = {}

    with open(file_path, "r") as tsv_file:
        for line in tsv_file:
            l = line.strip().split("\t")
            if len(l) == 2:
                fasta_files[l[0]] = l[1]

    return fasta_files


def parse_fastas(file_or_directory):
    """Work out if file_or_directory is a fasta-file, a tsv-file, or directory;
    parse with either parse_tsv_fasta() or fastas_from_directory()

    Args:
        file_or_directory (str): filepath for fasta, TSV file, or directory of FASTA files

    Returns:
        fasta_files (dict):
            key (str): file name minus extension
            value (str): filepath
    """
    fasta_files = {}

    if os.path.isfile(file_or_directory):
        if file_or_directory.lower().endswith(
            (".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn")
        ):
            fasta_files[
                os.path.splitext(os.path.basename(file_or_directory))[0]
            ] = file_or_directory
        elif file_or_directory.lower().endswith(".tsv"):
            fasta_files = parse_tsv_fasta(file_or_directory)
        else:
            print(f"Unsupported file format: {file_or_directory}")
    elif os.path.isdir(file_or_directory):
        fasta_files = fastas_from_directory(file_or_directory)
    else:
        print(f"Input not recognized: {file_or_directory}")

    return fasta_files


def write_fastas_tsv(fasta_dict, fastas_tsv):
    """Write a fasta_files dictionary to a TSV file

    Args:
        fasta_dict (dict):
            key (str): file name/prefix
            value (str): filepath
        fastas_tsv (str): filepath of TSV file for writing

    Returns:
        None
    """
    with open(fastas_tsv, "w") as tsv_file:
        for key, value in fasta_dict.items():
            tsv_file.write(f"{key}\t{value}\n")


def combine_fastas(fasta_dict, fasta_file):
    """Concatenate fasta files in fasta dictionary, adding the fasta ref name (key)
    as a prefix for the contig IDs

    Args:
        fasta_dict (dict):
            key (str): file name/prefix
            value (str): filepath
        fasta_file (str): Filepath for new concatenated fasta file

    Returns:
        None
    """
    with open(fasta_file, "w") as out_fasta:
        for ref_name, filepath in fasta_dict.items():
            with open(filepath, "r") as in_fasta:
                for line in in_fasta:
                    if line.startswith(">"):
                        line = line.replace(">", ">" + ref_name + ":")
                    out_fasta.write(line)
