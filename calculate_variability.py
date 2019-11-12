#!/usr/bin/env python3
"""
Calculates the variability in the DNA sequence for the fasta file

Discussed homework with Jeannette Rustin
"""

__author__ = "Yuta Sakai"

import scipy.signal as ssg


def main():
    fasta_file = open("/Users/m006703/Class/CSCI5481/Homework4/Homework4-seqs.fna", "r")
    variability_file = open("/Users/m006703/Class/CSCI5481/Homework4/variability.csv", "w")
    variability_file.write("Base,Variability,Smoothed Variability\n")

    # Collect the identifiers and sequences from the fasta file
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

    # Calculate the variability of the sequence
    index = 0
    variability_list = []
    while index in range(0, len(sequences[0])):
        variability_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for sequence in sequences:
            base = sequence[index]
            # Only keep track of A, C, G, T and ignore gaps, S, Y and Ns
            if base not in variability_dict.keys():
                continue
            else:
                variability_dict[base] += 1
        # Find the most common base and the number of occurance
        max_value = max(variability_dict.values())
        # FIXME
        # Not sure if we need to know the base
        for base, value in variability_dict.items():
            if value == max_value:
                most_common_base = base
        variability = (max_value / len(identifiers)) * 100
        variability_list.append(variability)
        index += 1

    # Smooth the variability using Savitzky-Golay filter
    smooted_variability_list = ssg.savgol_filter(variability_list, 51, 3)

    # Write the results to csv file
    index = 0
    while index in range(0, len(variability_list)):
        base = index + 1
        variability_file.write(str(base) + "," + str(variability_list[index]) + "," +
                               str(smooted_variability_list[index]) + "\n")
        index += 1

    fasta_file.close()
    variability_file.close()
    print("Script is done running")


if __name__ == "__main__":
    main()
