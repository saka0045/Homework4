#!/usr/bin/env python3

"""
Python script for Fall 2019 CSCI5481 Homework 3

Discussed homework with Jeannette Rustin
"""

__author__ = "Yuta Sakai"

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputFile", dest="input_file", required=True,
        help="input fasta file"
    )

    args = parser.parse_args()

    # Get the path to the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = script_dir + "/"

    # Get the path to the fasta file
    input_file_path = os.path.abspath(args.input_file)

    # Open the input fasta file
    fasta_file = open(input_file_path, "r")

    # Extract the identifiers and sequences from the fasta file
    identifiers = []
    sequences = []

    for line in fasta_file:
        if line.startswith(">"):
            line = line.rstrip()
            # Remove the ">" from the string and add it to the identifiers list
            identifiers.append(line[1:])
        else:
            line = line.rstrip()
            sequences.append(line)

    # Calculate the genetic distance between every pair of sequences and store in a dictionary
    genetic_distance_dict, genetic_distance_file = create_genetic_distance_dict_and_file(identifiers, script_dir,
                                                                                         sequences)

    # Initialize by starting at the last root possible (nTips + [nTips - 1]) being the first node to be joined
    # Changed this to (nTips + [nTips - 1]) for hw4
    root = len(identifiers) + (len(identifiers) - 2)
    edges_dict = {}

    # Loop until all of the neighbors are joined
    while root > len(identifiers):
        # Calculate the Q matrix and store the Q distances in a dictionary
        q_matrix_dict = calculate_q_distance(genetic_distance_dict)

        # Find the minimum distance in the Q matrix
        min_q_distance, min_q_distance_key, min_q_distance_partner = find_neighbors(q_matrix_dict)

        # Join the neighbors and calculate the distance to the newly formed root
        edges_dict = join_neighbors(edges_dict, genetic_distance_dict, min_q_distance_key, min_q_distance_partner, root)

        # Update genetic_distance_dict
        genetic_distance_dict = update_genetic_distance(genetic_distance_dict, min_q_distance_key,
                                                        min_q_distance_partner, root)
        # Subtract 1 from the root number and repeat
        root -= 1

    # For odd number of tips, there will be 2 nodes left over
    # Add the distance between those nodes and add it to the edges_dict
    if len(genetic_distance_dict) == 2:
        edges_dict[str(root + 1)][str(root + 2)] = genetic_distance_dict[str(root + 1)][str(root + 2)]

    # Make the edges.txt file
    edges_file = open(script_dir + "edges.txt", "w")

    make_edges_file(edges_dict, edges_file, identifiers)

    edges_file.close()

    # Construct the tree in NEWICK format
    final_newick_string = make_newick_tree(identifiers, script_dir)

    # Write the NEWICK format to the tree.txt file
    newick_tree_file = open(script_dir + "tree.txt", "w")
    newick_tree_file.write(final_newick_string + "\n")

    # Close out files
    fasta_file.close()
    genetic_distance_file.close()
    newick_tree_file.close()

    print("Script is done running")
    print("Genetic distances file is saved at: " + script_dir + "genetic_distances.txt")
    print("Edge file is saved at: " + script_dir + "edges.txt")
    print("Tree in NEWICK format is saved at: " + script_dir + "tree.txt")


def make_newick_tree(identifiers, script_dir):
    """
    Function to construce the tree in NEWICK format. This function will use the edges.txt file to traverse the tree
    in postorder. Splits each tree coming off of the root to create a sub tree in NEWICK format and concatenates
    the string for the final tree
    :param identifiers:
    :param script_dir:
    :return:
    """
    # Use the edges.txt file to traverse the tree in postorder
    edges_file = open(script_dir + "edges.txt", "r")
    final_newick_string = "("
    for line in edges_file:
        line = line.rstrip()
        line_item = line.split("\t")
        ancestral_node = line_item[0]
        descendant_node = line_item[1]
        node_distance = line_item[2]
        # Check to see if the ancestral node is the root
        if ancestral_node == str(len(identifiers) + 1):
            try:
                # Next time it comes across the root, one branch from the root would have been traversed in preorder
                if len(ancestral_node_list) > 0:
                    # Below portion is done for each tree from the root, except for the final tree visited
                    # Reverse the order of that was traversed in preorder to traverse in postorder
                    newick_string = ""
                    ancestral_node_list.reverse()
                    descendant_node_list.reverse()
                    node_distance_list.reverse()
                    # Keep count of if the ancestral node was visited for the first time or not
                    ancestral_node_count = []
                    # Keep track of the highest node (smallest number) that the other child hasn't been visited
                    open_node = ""
                    newick_string = make_newick_substring(ancestral_node_count, ancestral_node_list, descendant_node_list,
                                                          identifiers, newick_string, node_distance_list, open_node)
                    # Append the portion of the tree visited from the root to the entire tree newick string
                    final_newick_string += newick_string
                    # Reset the node lists and start the next tree from the root
                    ancestral_node_list = []
                    descendant_node_list = []
                    node_distance_list = []
                    ancestral_node_list.append(ancestral_node)
                    descendant_node_list.append(descendant_node)
                    node_distance_list.append(node_distance)
            # Since the first item's ancestral node is the root, it will error out
            except NameError:
                # Create the lists and add the root information
                ancestral_node_list = []
                descendant_node_list = []
                node_distance_list = []
                ancestral_node_list.append(ancestral_node)
                descendant_node_list.append(descendant_node)
                node_distance_list.append(node_distance)
        # If the ancestral node is not a root, add the node information to the lists
        else:
            ancestral_node_list.append(ancestral_node)
            descendant_node_list.append(descendant_node)
            node_distance_list.append(node_distance)
    # Traverse the last tree for the root in postorder and create the NEWICK format
    newick_string = ""
    ancestral_node_list.reverse()
    descendant_node_list.reverse()
    node_distance_list.reverse()
    ancestral_node_count = []
    open_node = ""
    newick_string = make_newick_substring(ancestral_node_count, ancestral_node_list, descendant_node_list,
                                          identifiers, newick_string, node_distance_list, open_node)
    # Append the NEWICK format of the final tree for the root to the entire tree string
    final_newick_string += newick_string
    # Close out the final_newick_string for the complete tree
    # If the final newick string ends with ",", replace it with ");"
    if final_newick_string.endswith(","):
        final_newick_string = final_newick_string[:-1] + ");"
    # Otherwise, close out the tree with ");"
    else:
        final_newick_string += ");"
    edges_file.close()
    return final_newick_string


def make_newick_substring(ancestral_node_count, ancestral_node_list, descendant_node_list,
                          identifiers, newick_string, node_distance_list, open_node):
    """
    Function to traverse one of the trees from the root in postorder and construct the results in NEWICK format
    :param ancestral_node_count:
    :param ancestral_node_list:
    :param descendant_node_list:
    :param identifiers:
    :param newick_string:
    :param node_distance_list:
    :param open_node:
    :return: Returns the NEWICK format string for the tree that traversed
    """
    for index in range(0, len(ancestral_node_list)):
        # If the descendant node is a tip
        if int(descendant_node_list[index]) <= len(identifiers):
            # Get the identifier name in identifiers list
            descendant_node_key = (int(descendant_node_list[index]) - 1)
            descendant_node_identifier = identifiers[descendant_node_key]
            # If this tip is the first child, open the parenthesis and add the tip information
            if ancestral_node_list[index] not in ancestral_node_count:
                newick_string += "(" + descendant_node_identifier + ":" + node_distance_list[
                    index] + ","
                # Add the ancestral node to the ancestral node list to
                # inform one of the child for this node has been visited
                ancestral_node_count.append(ancestral_node_list[index])
            else:
                # If this tip is the second child, add the tip information and close the parenthesis
                newick_string += descendant_node_identifier + ":" + node_distance_list[index] + ")"
        # If the descendant node is a node
        elif int(descendant_node_list[index]) > len(identifiers):
            # If this child node is the first branch of the ancestral node visited:
            if ancestral_node_list[index] not in ancestral_node_count:
                newick_string += ":" + node_distance_list[index] + ","
                # First time encountering a node with only one child visited (aka, the first node in the list),
                # add "(" at the beginning of the newick string
                if open_node == "":
                    newick_string = "(" + newick_string
                    # Keep track of the open node
                    open_node = ancestral_node_list[index]
                # If you come across a node with only one child visited (open node) that is higher in the
                # tree (less number) than the previous open node,
                # add "(" at the very beginning of the string and update the open node
                elif int(ancestral_node_list[index]) < int(open_node):
                    # Don't add any parenthesis for the root
                    if (ancestral_node_list[index]) == str((len(identifiers) + 1)):
                        continue
                    else:
                        newick_string = "(" + newick_string
                        open_node = ancestral_node_list[index]
                # If the open node is at a lower position in the tree (higher number),
                # add "(" at the last occurrence of an open node "("
                elif int(ancestral_node_list[index]) > int(open_node):
                    index_for_parenthesis = newick_string.rfind("(")
                    newick_string = (newick_string[:(index_for_parenthesis + 1)] + "(" +
                                     newick_string[(index_for_parenthesis + 1):])
                # Append the ancestral node to the list of nodes with only one child visited
                # Skip this step if the ancestral node is a root
                if ancestral_node_list[index] != str((len(identifiers) + 1)):
                    ancestral_node_count.append(ancestral_node_list[index])
                else:
                    continue
            else:
                # If this node already had the other child visited:
                newick_string += ":" + node_distance_list[index] + ")"
    return newick_string


def make_edges_file(edges_dict, edges_file, identifiers):
    """
    Make the edges.txt file
    See code for comments
    :param edges_dict:
    :param edges_file:
    :param identifiers:
    :return:
    """
    # Keep track of which descendant nodes and tips that was already written to the edges.txt
    descendant_nodes_written = []
    tips_written = []
    # Also keep track of how many times the ancesteral node was written
    # Each ancestral node should have two descendants so it needs to be written twice
    ancestral_nodes_only_written_once = []
    # Start at the internal node root, ntips + 1
    ancestral_node_to_write = str((len(identifiers) + 1))
    # Calculate the number of nodes to write, which is all of the nested keys in the edges_dict
    number_of_nodes_to_write = 0
    for key in edges_dict.keys():
        number_of_nodes_to_write += len(edges_dict[key])
    # Loop until all of the nodes are written
    while number_of_nodes_to_write > 0:
        edges_file.write(ancestral_node_to_write + "\t")
        # Append the ancestral node to the written ancestral node list
        ancestral_nodes_only_written_once.append(ancestral_node_to_write)
        # Obtain the descendant node information
        descendant_nodes_dict = edges_dict[ancestral_node_to_write]
        descendant_nodes_list = list(descendant_nodes_dict.keys())
        # See if the descendant tip/node was already written to the file as a descendant tip/node
        for index in range(0, len(descendant_nodes_list)):
            if descendant_nodes_list[index] in descendant_nodes_written:
                continue
            elif descendant_nodes_list[index] in tips_written:
                continue
            # If not, write to the file
            else:
                descendant_node_to_write = descendant_nodes_list[index]
                edges_file.write(descendant_node_to_write + "\t")
                # Write the distance from the ancestral node to the descendant tip/node
                edge_length = descendant_nodes_dict[descendant_node_to_write]
                edges_file.write(str(edge_length) + "\n")
                # Decrease the number of nodes to write by 1
                number_of_nodes_to_write -= 1
                # If the descendant node number is greater than nTips, it is a node
                if int(descendant_node_to_write) > len(identifiers):
                    descendant_nodes_written.append(descendant_node_to_write)
                    # If the ancestral node was written as many times as it has descendant tip/node,
                    # remove from the ancestral node list
                    if ancestral_nodes_only_written_once.count(ancestral_node_to_write) == \
                            len(edges_dict[ancestral_node_to_write]):
                        # Delete as many times as it needs to
                        # Most ancestral nodes have 2 descendant tips/nodes but for odd number of tips,
                        # An ancestral node can have 3 descendant tips/nodes
                        for loop_count in enumerate(edges_dict[ancestral_node_to_write].keys()):
                            ancestral_nodes_only_written_once.remove(ancestral_node_to_write)
                        # Make the next root the next ancestral node to write
                        ancestral_node_to_write = descendant_node_to_write
                    # If not written twice, make the descendant the next ancestral node to write
                    else:
                        ancestral_node_to_write = descendant_node_to_write
                    break
                # If the descendant node number is less than or equal to nTips, it is a tip
                else:
                    tips_written.append(descendant_node_to_write)
                    # If the ancestral node was written as many times as it has descendant tip/node,
                    # remove from the ancestral node list
                    if ancestral_nodes_only_written_once.count(ancestral_node_to_write) == \
                            len(edges_dict[ancestral_node_to_write]):
                        for loop_count in enumerate(edges_dict[ancestral_node_to_write].keys()):
                            ancestral_nodes_only_written_once.remove(ancestral_node_to_write)
                        # Make the next ancestral node to write the last item in the list and continue the loop
                        try:
                            ancestral_node_to_write = ancestral_nodes_only_written_once[-1]
                        except IndexError:
                            number_of_nodes_to_write = 0
                            break
                    break


def update_genetic_distance(genetic_distance_dict, min_q_distance_key, min_q_distance_partner, root):
    """
    Update the genetic_distance_dict with the newly formed root
    :param genetic_distance_dict:
    :param min_q_distance_key:
    :param min_q_distance_partner:
    :param root:
    :return: Updated genetic_distance_dict
    """
    # Add the new root to the dictionary
    genetic_distance_dict[str(root)] = {}
    # Calculate the distance to the new node
    for key in genetic_distance_dict.keys():
        # No need to update the tip/root that is being joined
        if key == min_q_distance_key or key == min_q_distance_partner:
            continue
        # No need to do this for the newly formed root also because we will be copying the distances as we go
        if key == str(root):
            continue
        else:
            genetic_distance_dict[key][str(root)] = (0.5 * (genetic_distance_dict[min_q_distance_key][key] +
                                                            genetic_distance_dict[min_q_distance_partner][key] -
                                                            genetic_distance_dict[min_q_distance_key][
                                                                min_q_distance_partner]))
            # Copy the same distance to the nested dictionary of the newly formed root
            genetic_distance_dict[str(root)][key] = genetic_distance_dict[key][str(root)]
            # Delete the joined tip/root from the nested dictionary
            del genetic_distance_dict[key][min_q_distance_key]
            del genetic_distance_dict[key][min_q_distance_partner]
    # Delete the joined tip/root from the genetic distance dictionary
    del genetic_distance_dict[min_q_distance_key]
    del genetic_distance_dict[min_q_distance_partner]
    return genetic_distance_dict


def join_neighbors(edges_dict, genetic_distance_dict, min_q_distance_key, min_q_distance_partner, root):
    """
    Join the neighbors, calculate the distance from the neighbors to the newly formed root
    and store the information in the edges_dict dictionary
    :param edges_dict: Updates the already existing edges_dict
    :param genetic_distance_dict:
    :param min_q_distance_key:
    :param min_q_distance_partner:
    :param root:
    :return: Updated edges_dict
    """
    # Calculate the distance to the new node from the first branch
    first_branch_distance = (0.5 * genetic_distance_dict[min_q_distance_key][min_q_distance_partner] +
                             (1 / (2 * (len(list(genetic_distance_dict.keys())) - 2))) *
                             (sum(genetic_distance_dict[min_q_distance_key].values()) -
                              sum(genetic_distance_dict[min_q_distance_partner].values())))
    # Calculate the second branch distance
    second_branch_distance = genetic_distance_dict[min_q_distance_key][min_q_distance_partner] - first_branch_distance
    # Store the joined node information, start at number 120 and decrease every cycle
    edges_dict[str(root)] = {min_q_distance_key: first_branch_distance,
                             min_q_distance_partner: second_branch_distance}
    return edges_dict


def find_neighbors(Q_matrix_dict):
    """
    Find the tip/root with the smallest Q distance (neighbors)
    :param Q_matrix_dict:
    :return:
    """
    # Initialize the overall minimum Q distance as 0
    min_q_distance = 0
    for key in Q_matrix_dict.keys():
        # Check the minimum Q distance for a given sequence
        min_q_distance_for_sequence = min(Q_matrix_dict[key].values())
        # If the minimum Q distance for the sequence is smaller than the overall min Q distance,
        # replace the overall minimum Q distance
        if min_q_distance_for_sequence < min_q_distance:
            min_q_distance = min_q_distance_for_sequence
            # Store the sequence (key) where this happens
            min_q_distance_key = key
            # See against which tip or root the min Q value was calculated relative to the key
            for nested_key, nested_value in Q_matrix_dict[key].items():
                if nested_value == min_q_distance:
                    min_q_distance_partner = nested_key
    return min_q_distance, min_q_distance_key, min_q_distance_partner


def calculate_q_distance(genetic_distance_dict):
    """
    Calculates the Q distance between all tip/root and store it in a nested dictionary.
    Uses the same tip/root number as the genetic_distance_dict
    :param genetic_distance_dict:
    :return:
    """
    q_matrix_dict = {}
    for first_key in genetic_distance_dict.keys():
        node_i = first_key
        # Created a nested dictionary inside the Q_matrix_dict, the key will be the tip or root number
        q_matrix_dict[node_i] = {}
        # Calculate the sum of distances to the ith sequence
        sum_distance_i_k = sum(genetic_distance_dict[node_i].values())
        for second_key in genetic_distance_dict.keys():
            node_j = second_key
            # No need to calculate Q for the same sequences
            if first_key == second_key:
                continue
            else:
                # Calculate the sum of distances to the jth sequence
                sum_distance_j_k = sum(genetic_distance_dict[node_j].values())
                # Calculate the Q distance
                Q_distance = ((len(list(genetic_distance_dict.keys())) - 2) * genetic_distance_dict[node_i][node_j] -
                              sum_distance_i_k - sum_distance_j_k)
                # Append the root or tip number as the key and Q distance as the value
                q_matrix_dict[node_i][node_j] = Q_distance
    return q_matrix_dict


def create_genetic_distance_dict_and_file(identifiers, script_dir, sequences):
    """
    The key for the genetic_distance_dict will represent the tip or root number and will start with 1,
    where 1 is a tip and represents the first sequence in the fasta file.
    Until the dictionary gets updated, it will start out with the 61 sequences in the fasta file.
    The internal root node will start with 62 since there are 61 sequences in the fasta file.
    There are 59 (61 - 2) nodes possible in this phylogenetic tree, so the first pair to be joined will be
    root 120 (61 + 59).
    :param identifiers:
    :param script_dir:
    :param sequences:
    :return:
    """
    genetic_distance_dict = {}
    for i in range(0, len(identifiers)):
        # Create a nested dictionary that will store the genetic distances with the tip number as key
        nested_genetic_distance_dict = {}
        for j in range(0, len(sequences)):
            genetic_distance = calculate_genetic_distance(sequences[i], sequences[j])
            nested_genetic_distance_dict[str(j + 1)] = genetic_distance
        genetic_distance_dict[str(i + 1)] = nested_genetic_distance_dict
    # Create a tab-delimited file of genetic distance in the directory where this script lives
    genetic_distance_file = open(script_dir + "genetic_distances.txt", "w")
    # Make the header for the file, start with a blank "cell"
    genetic_distance_file.write("\t")
    genetic_distance_file.write("\t".join(identifiers))
    genetic_distance_file.write("\n")
    # Iterate through all the items in genetic_distance_dict and write it to the result file
    for key, val in genetic_distance_dict.items():
        sequence = identifiers[int(key) - 1]
        genetic_distance_file.write(sequence + "\t")
        # Iterate through the nested dictionary to write the genetic distances to a file
        for nested_value in genetic_distance_dict[key].values():
            genetic_distance_file.write(str(nested_value) + "\t")
        genetic_distance_file.write("\n")
    return genetic_distance_dict, genetic_distance_file


def calculate_genetic_distance(first_sequence, second_sequence):
    """
    Calculates the genetic distance between two sequences. If the two sequences is not the same length,
    the program will exit with an error
    :param first_sequence:
    :param second_sequence:
    :return: genetic_distance
    """
    similarity_score = 0
    if len(first_sequence) != len(second_sequence):
        sys.exit("Length of two sequences needs to be the same!")
    elif len(first_sequence) == len(second_sequence):
        for position in range(0, len(first_sequence)):
            if first_sequence[position] == second_sequence[position]:
                similarity_score += 1
        genetic_distance = 1 - (similarity_score / len(first_sequence))
        return genetic_distance


if __name__ == "__main__":
    main()
