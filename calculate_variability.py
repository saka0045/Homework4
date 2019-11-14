#!/usr/bin/env python3
"""
Calculates the variability in the DNA sequence for the fasta file

Discussed homework with Jeannette Rustin
"""

__author__ = "Yuta Sakai"

import argparse
import os
import scipy.signal as ssg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputFile", dest="input_file", required=True,
        help="full path to the Homework4-seqs.fna fasta file"
    )
    parser.add_argument(
        "-o", "--outputDir", dest="output_dir", required=True,
        help="Directory to save output file"
    )

    args = parser.parse_args()

    input_file = os.path.abspath(args.input_file)
    out_path = os.path.abspath(args.output_dir)

    # Add / at the end if it is not included in the output path
    if out_path.endswith("/"):
        out_path = out_path
    else:
        out_path = out_path + "/"

    fasta_file = open(input_file, "r")
    variability_file = open(out_path + "variability.csv", "w")
    variable_regions_file = open(out_path + "variable_regions.txt", "w")

    # Make headers for the variability file
    variability_file.write("Base,Variability,Smoothed Variability\n")

    # Collect the identifiers and sequences from the fasta file
    identifiers, sequences = parse_fasta(fasta_file)

    # Calculate the variability of the sequence
    sequence_identity_list = calculate_sequence_variability(identifiers, sequences)

    # Smooth the variability using Savitzky-Golay filter
    smoothed_sequence_identity_list = ssg.savgol_filter(sequence_identity_list, 61, 3)

    # Write the results to csv file
    make_csv(smoothed_sequence_identity_list, variability_file, sequence_identity_list)

    # Interrogate the smoothed_sequence_identity_list and identify the variable region
    end_variable_region_list, start_variable_region_list = identify_variable_regions(smoothed_sequence_identity_list)

    # Make the variability region tab delimited file
    for index, start_region in enumerate(start_variable_region_list):
        variable_regions_file.write(str(start_region) + "\t" + str(end_variable_region_list[index]) + "\n")

    fasta_file.close()
    variability_file.close()
    variable_regions_file.close()
    print("Script is done running")


def identify_variable_regions(smoothed_sequence_identity_list):
    """
    Function to identify the variable regions. Variable region is defined as a region where the sequence identity drops
    below 75%. With 4 nucleotides there is 25% chance of any nucleotide to be at that sequence by chance.
    Chance of any 3 nucleotides to be at that position is up to 75%, therefore if the sequence identity at the given
    position is greater than 75%, the sequence is conserved.
    :param smoothed_sequence_identity_list:
    :return:
    """
    start_variable_region = "not yet initialized"
    end_variable_region = "not yet initialized"
    start_variable_region_list = []
    end_variable_region_list = []
    for base, smoothed_sequence_identity in enumerate(smoothed_sequence_identity_list, start=1):
        # Ignore the beginning portion until the smoothed sequence identity becomes higher than 75%
        if smoothed_sequence_identity >= 75 and start_variable_region == "not yet initialized" and end_variable_region \
                == "not yet initialized":
            # Initialize the variable start region once sequence identity is higher than 75%
            start_variable_region = ""
        # Keep going down the list until the sequence identity drops below 75%
        if start_variable_region == "":
            if smoothed_sequence_identity < 75:
                # Once you enter the variable region, make this the start variable position and
                # initialize the end position
                start_variable_region = base
                end_variable_region = ""
            else:
                continue
        # Keep going down the list of sequence identity until the identity becomes greater than 75%
        if start_variable_region != "" and end_variable_region == "":
            # When you hit greater than 75% identity, you will be exiting the variable region
            if smoothed_sequence_identity >= 75:
                # Record the end position
                end_variable_region = base
                # Once you have a start and stop for the variable region, add them to the lists
                start_variable_region_list.append(start_variable_region)
                end_variable_region_list.append(end_variable_region)
                # Initialize the value and start over
                start_variable_region = ""
                end_variable_region = ""
            else:
                continue
    return end_variable_region_list, start_variable_region_list


def make_csv(smoothed_sequence_identity_list, variability_file, sequence_identity_list):
    index = 0
    while index in range(0, len(sequence_identity_list)):
        base = index + 1
        variability_file.write(str(base) + "," + str(sequence_identity_list[index]) + "," +
                               str(smoothed_sequence_identity_list[index]) + "\n")
        index += 1


def parse_fasta(fasta_file):
    identifiers = []
    sequences = []
    for line in fasta_file:
        line = line.rstrip()
        if line.startswith(">"):
            # Get rid of the ">"
            line = line[1:]
            identifiers.append(line)
        else:
            sequences.append(line)
    return identifiers, sequences


def calculate_sequence_variability(identifiers, sequences):
    index = 0
    sequence_identity_list = []
    while index in range(0, len(sequences[0])):
        sequence_identity_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for sequence in sequences:
            base = sequence[index]
            # Only keep track of A, C, G, T and ignore gaps, S, Y and Ns
            if base not in sequence_identity_dict.keys():
                continue
            else:
                sequence_identity_dict[base] += 1
        # Find the most common base and the number of occurance
        max_value = max(sequence_identity_dict.values())
        # Not sure if we need to know the base
        for base, value in sequence_identity_dict.items():
            if value == max_value:
                most_common_base = base
        sequence_identity = (max_value / len(identifiers)) * 100
        sequence_identity_list.append(sequence_identity)
        index += 1
    return sequence_identity_list


if __name__ == "__main__":
    main()
