#!/usr/bin/env python3
"""
Calculates the variability in the DNA sequence for the fasta file
"""


def main():
    fasta_file = open("/Users/m006703/Class/CSCI5481/Homework4/Homework4-seqs.fna", "r")

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
    while index in range(0, len(sequences[0])):
        variability_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for sequence in sequences:
            base = sequence[index]
            # Only keep track of A, C, G, T and ignore gaps, S, Y and Ns
            if base not in variability_dict.keys():
                continue
            else:
                variability_dict[base] += 1
        print(variability_dict)
        # Find the most common base and the number of occurence
        max_value = max(variability_dict.values())
        # FIXME
        # Not sure if we need to know the base
        for base, value in variability_dict.items():
            if value == max_value:
                most_common_base = base
        print(most_common_base + " " + str(max_value))
        variability = (max_value / len(identifiers)) * 100
        print(variability)
        index += 1

    fasta_file.close()


if __name__ == "__main__":
    main()
