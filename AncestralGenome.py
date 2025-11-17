#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 12:53:14 2024

@author: user
"""

### this script will reconstruct ancestral genomes
import os
import ete3
import re
import numpy as np
import time
import argparse
from multiprocessing import Pool, cpu_count
import csv
from common_functions import FilterHogs, FindDuplications

### functions ################################################################

## function to plot a tree
# needs a tree, a path to orthofind folder, and a filename
def PlotTree(t, path, filename):
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    ts = ete3.TreeStyle()
    lstyle = ete3.NodeStyle()
    lstyle["fgcolor"] = "blue"
    lstyle["size"] = 1.5
    nstyle = ete3.NodeStyle()
    nstyle["fgcolor"] = "red"
    nstyle["size"] = 3
    for n in t.traverse():
        if n.is_leaf():
            n.set_style(lstyle)
        else:
            n.add_face(ete3.TextFace(n.name), column=0)
            n.set_style(nstyle)
    t.render(os.path.join(path, filename + ".png"), w=1920, units='px', tree_style=ts)

### functions for processing species tree, and mapping to gene tree  ##########

#ortho_folder_path = "/Users/user/Dropbox/ORTHOFINDER/PROJECT_GAIN_LOSS/Feng2024/Results_Jan08"
#gene_tree_file = "N0.HOG0000065_tree.txt"


def load_gene_tree_from_big_file(og, tree_file_path):
    with open(tree_file_path, "r") as f:
        for line in f:
            if line.startswith(og + ": "):  # OG ID followed by a space
                _, tree_str = line.strip().split(maxsplit=1)
                return ete3.Tree(tree_str, quoted_node_names=True, format=1)
    raise ValueError(f"OG {og} not found in tree file.")


## function to get gene names for ancestral genomes
## does it for one OG
def GetAncestralGenes(OG, node, species_tree, ortho_folder_path, species_names):
    # to store all selected sequences
    all_selected_sequences = []
    og = OG

    ## load the gene tree
    tree_file_path = os.path.join(ortho_folder_path, "Resolved_Gene_Trees/Resolved_Gene_Trees.txt")
    gene_tree = load_gene_tree_from_big_file(og, tree_file_path)
    #gene_tree_file = og + "_tree.txt" # fungi8
    #gene_tree = ete3.Tree(os.path.join(ortho_folder_path, "HOG_Gene_Trees/", gene_tree_file), quoted_node_names=True, format=1)
    
    # remove species that we dont need
    focal_node = node.name
    target_node = species_tree.search_nodes(name=focal_node)[0]
    target_species = target_node.get_leaves()
    #### make this
    target_species = np.array(list(leaf.name for leaf in target_species))
    # get species below target
   
    # get list of all leaves in gene tree
    all_leavesN = gene_tree.get_leaves() # in node form
    all_leaves = np.array(list(leaf.name for leaf in all_leavesN)) # in array form
    # get corresponding list of species
    #leaves1 = set("_".join(leaf.name.split("_")[:2]) for leaf in node.children[0].get_leaves())
    #leaves2 = set(re.sub(r"_[^_]*$", "", leaf.name) for leaf in node.children[1].get_leaves())
    #all_species = np.array(list(re.sub(r"_[^_]*$", "", leaf.name) for leaf in gene_tree.get_leaves()))
    all_species = np.array(["_".join(leaf.name.split("_")[:2]) for leaf in gene_tree.get_leaves()])
    #print(all_species)
    # find the leaves which match target species
    indicheese = np.where(np.isin(all_species, target_species))
    leaves_to_stay = [all_leavesN[i] for i in indicheese[0]]
    #print("all leaves")
    #print(all_leaves)
    #print("all species")
    #print(all_species)
    #print("target species")
    #print(target_species)
    #print("leaves to stay")
    #print(leaves_to_stay)
    # if there are no leaves in the species below target node, skip
    if not leaves_to_stay:
        return None
    
    # prune the tree
    # but calculate duplications first, as pruning changes the node
    duplications = FindDuplications(gene_tree, species_tree, species_names)
    # only account for duplications with >=0.5 support
    # if not duplications.empty:
    #     duplications = duplications[duplications['support'] >= 0.5]
    duplications = [d for d in duplications if d['support'] >= 0.5]
    gene_tree.prune(leaves_to_stay)
    
    ## deal with duplications that happened after target node
    # I will randomly pick one child branch to remove
    # define nodes that are after target node
    desc_nodes = target_node.get_descendants()
    desc_nodes = np.array(list(leaf.name for leaf in desc_nodes))
    
    # find which duplications are after target node
    dupes_to_go = []
    # Only do if there are duplications
    if duplications:
        index1 = [i for i, d in enumerate(duplications) if d['speciestree_node'] in desc_nodes]
        # Randomly pick to drop leaves1 or leaves2
        choice_column = np.random.choice(['leaves1', 'leaves2'], size=len(index1))
        for idx, col in zip(index1, choice_column):
            dupes_to_go.append(duplications[idx][col])        
    
    # get list of leaves to keep
    all_leavesN = gene_tree.get_leaves() # in node form
    all_leaves = np.array(list(leaf.name for leaf in all_leavesN)) # in array form
    keep_leaves = [all_leavesN[i] for i in np.where(~np.isin(all_leaves, np.array(dupes_to_go, dtype=object)))[0]]
    gene_tree.prune(keep_leaves)
 
    # deal with duplications that happened before focal node
    # we need to include the right number of genes in our ancestral genome
    
    # get duplications that are in ancestor of target node
    ancestors = [focal_node]
    ## look at the parent node
    up_node = target_node.up
    ## while the parent node exists, grab the node names of ancestors
    while up_node:
        ancestors.append(up_node.name)
        up_node = up_node.up
    
    expected_copies = 1
    # Only do if there are duplications
    if duplications:
        it = [i for i, d in enumerate(duplications) if d['speciestree_node'] in ancestors and d['support'] >= 0.5]
        duplications_anc = [duplications[i] for i in it]
        # Count how many duplications there are on the way to the target node
        expected_copies += len(duplications_anc)
    
    # we will calculate tip to root distance for each leaf, and pick the median
    # we then repeat until we have the right number of sequences
    
    # Calculate root-to-tip distance for each leaf
    # make this a loop to get distance for all leaves
    # get list of all leaves in gene tree
    all_leavesN = gene_tree.get_leaves() # in node form
    all_leaves = np.array(list(leaf.name for leaf in all_leavesN)) # in array form
    distance = []
    for leaf in gene_tree.get_leaves():
        distance.append(leaf.get_distance("n0"))

    ## get the sequence with median branch length,
    ## and repeat until you have the expected number of copies
    selected_sequences = []  # Store selected medians here
    # Make mutable copies of the arrays
    current_distances = distance.copy()
    current_leaves = all_leaves.copy()
    # if we have enough expected copies, use the median approach
    # if we just fall short, so be it
    for _ in range(expected_copies):
        if len(current_distances) > 0:
            # Get index of the current median distance
            median_index = np.argpartition(current_distances, len(current_distances) // 2)[len(current_distances) // 2]
            # Append the median leaf to the results
            selected_sequences.append(current_leaves[median_index])
            # Remove the selected median from both arrays for the next iteration
            current_distances = np.delete(current_distances, median_index)
            current_leaves = np.delete(current_leaves, median_index)
        else:
            temp=1
    all_selected_sequences.append(selected_sequences)
    #print(len(selected_sequences))
    return all_selected_sequences

## function to open species fastas, and make a dictionary mapping geneID to seq
### currently requires input fastas to be in a folder input_fastas - but can switch
## to using the ones in working directory? (numeric coded)
def MakeAllSpeciesFastaDict(all_species, ortho_folder_path):
    # all_fasta_dicts contains each species' gene names and sequences dictionaries
    all_fasta_dicts = {}
    # Loop through all species
    for species_name in all_species:
        #print(species_name)
        # Define the file path for the current species
        for filename in os.listdir(os.path.join(ortho_folder_path, "../../")):
            normalized_file = filename.replace(".", "_")
            last_under = normalized_file.rfind('_')
            normalized_file = normalized_file[:last_under]
            if normalized_file == f"{species_name}":
                #print("match found")
                species_partial_path = f"../../{filename}"
                
        species_path = os.path.join(ortho_folder_path, species_partial_path)
        # Open the FASTA file and read it into a string
        with open(species_path, 'r') as file:
            fasta_data = file.read()
        # Split the data into entries (each entry starts with '>')
        entries = fasta_data.split('>')
        # Initialize a dictionary to store gene names and sequences
        fasta_dict = {}
        # Process each entry
        for entry in entries:
            if entry.strip():  # Check if the entry is not just whitespace
                lines = entry.splitlines()  # Split the entry into lines
                gene_name = lines[0].split()[0]  # The first line is the gene name
                sequence = ''.join(lines[1:])  # Join the remaining lines to form the sequence
                fasta_dict[gene_name] = sequence  # Add to dictionary
        # Store the dictionary using the species name as key
        all_fasta_dicts[species_name] = fasta_dict
    return all_fasta_dicts       

## function to open species fastas, and make a dictionary mapping geneID to seq
def MakeAncestralSpeciesFastaDict(all_species, ortho_folder_path):
    # all_fasta_dicts contains each species' gene names and sequences dictionaries
    all_fasta_dicts = {}
    gene_counts = {}
    total_lengths = {}
    orthogroup_counts = {}
    # Loop through all species
    for species_name in all_species:
        # Define the file path for the current species
        species_partial_path = f"{species_name}.fasta"
        species_path = os.path.join(ortho_folder_path, species_partial_path)
        # Open the FASTA file and read it into a string
        with open(species_path, 'r') as file:
            fasta_data = file.read()
        # Split the data into entries (each entry starts with '>')
        entries = fasta_data.split('>')
        # Initialize a dictionary to store gene names and sequences
        fasta_dict = {}
        total_length = 0
        gene_count = 0
        # Process each entry
        for entry in entries:
            if entry.strip():  # Check if the entry is not just whitespace
                lines = entry.splitlines()  # Split the entry into lines
                gene_name = lines[0].split()[0]  # The first line is the gene name
                sequence = ''.join(lines[1:])  # Join the remaining lines to form the sequence
                fasta_dict[gene_name] = sequence  # Add to dictionary
                total_length += len(sequence)
                gene_count += 1
                # count how many times an orthogroup appers
                orthogroup = gene_name.split('_')[1]
                #print(f"Gene Name: {gene_name}, Orthogroup: {orthogroup}")  
                # add this orthogroup, if not already there
                if orthogroup not in orthogroup_counts:
                    orthogroup_counts[orthogroup] = {}
                # either set count to 0 (if its a new OG), or add 1 to count
                if species_name not in orthogroup_counts[orthogroup]:
                    orthogroup_counts[orthogroup][species_name] = 0
                orthogroup_counts[orthogroup][species_name] += 1
        # Store the dictionary using the species name as key
        all_fasta_dicts[species_name] = fasta_dict
        gene_counts[species_name] = gene_count
        total_lengths[species_name] = total_length
    return all_fasta_dicts, gene_counts, total_lengths, orthogroup_counts 

## function to build ancestral genome, once we know which genes are in it
def BuildAncestralGenome(species_list, gene_list, orthogroup_list, fasta_dicts, focal_node):
    fasta_content = ""
    # Create a dictionary to count the occurrences of each orthogroup
    orthogroup_counts = {}
    # First, count the number of genes in each orthogroup
    for orthogroup in orthogroup_list:
        if orthogroup in orthogroup_counts:
            orthogroup_counts[orthogroup] += 1
        else:
            orthogroup_counts[orthogroup] = 1
    # Create a dictionary to keep track of the current gene number for each orthogroup
    current_gene_num = {}
    for species_id, gene_id, orthogroup in zip(species_list, gene_list, orthogroup_list):
        if orthogroup in current_gene_num:
            current_gene_num[orthogroup] += 1
        else:
            current_gene_num[orthogroup] = 1
        
        ####
        try:
            gene_sequence = fasta_dicts[species_id][gene_id]
        except KeyError as e:
            if species_id not in fasta_dicts:
                raise KeyError(f"Species ID '{species_id}' not found in fasta_dicts.") from e
            elif gene_id not in fasta_dicts[species_id]:
                raise KeyError(f"Gene ID '{gene_id}' not found under species ID '{species_id}'.") from e
            else:
                # This case should never be reached, but it's here for completeness
                raise KeyError(f"Unexpected KeyError with species ID '{species_id}' and gene ID '{gene_id}'.") from e
                
        ####
        
        gene_sequence = fasta_dicts[species_id][gene_id]
        # Get the total number of genes in this orthogroup
        gene_num2 = orthogroup_counts[orthogroup]
        # Get the current gene number
        gene_num = current_gene_num[orthogroup]
        fasta_content += f">{focal_node}_{orthogroup}_{gene_num}_{gene_num2}_{gene_id}\n"
        fasta_content += f"{gene_sequence}\n"
    return fasta_content

### function to save a fasta
def WriteAncestralFasta(focal_node, ancestral_genome, ortho_folder_path):
    output_fasta_path = os.path.join(ortho_folder_path,"AncestralGenomes/", focal_node + ".fasta")
    with open(output_fasta_path, 'w') as fasta_file:
        # Write the header line
        fasta_file.write(ancestral_genome)    

# Process orthogroups
def ProcessOrthogroupCurrent(index, gains_current, node, species_tree, ortho_folder_path, species_names):
    OG = gains_current[index]['Orthogroup']
    return OG, GetAncestralGenes(OG, node, species_tree, ortho_folder_path, species_names)

# the main function, which pulls everything together
# allowing ancestral genomes to be build and saved
# for all nodes in the species tree
def AncestralGenome(node, species_tree, OG_file, ortho_folder_path, gains, species_names):
    ## define a focal node in species tree
    focal_node = node.name
    # find nodes that are ancestor to our focal node
    ## set the node we are working with
    target_node = species_tree.search_nodes(name=focal_node)[0]
    ## variable to store ancestor nodes
    ancestors = [focal_node]
    ## look at the parent node
    up_node = target_node.up
    ## while the parent node exists, grab the node names of ancestors
    while up_node:
        ancestors.append(up_node.name)
        up_node = up_node.up
    
    ## filter gains to only include HOGs which have a gain_node in ancestors
    #gains_current = gains[gains['Gain Node'].isin(ancestors)].reset_index(drop=True)
    gains_current = [gain for gain in gains if gain['Gain Node'] in ancestors]
    
    # get species below target
    target_species = [leaf.name for leaf in target_node.get_leaves()]    
    with Pool(processes=cpu_count()-2) as pool:
        results = pool.starmap(ProcessOrthogroupCurrent, [(index, gains_current, node, species_tree, ortho_folder_path, species_names) for index in range(len(gains_current))])
        
    ancestral_sequence = []
    orthogroup_names = []

    for OG, result in results:
        if result is not None:
            ancestral_sequence.append(result)
            orthogroup_names.extend([OG] * len(result[0]))
        
    ### manipulate format to make useful, try flatten
    if ancestral_sequence:
        ancestral_gene_ids = np.concatenate(ancestral_sequence, axis=None).ravel().tolist()
    else:
        print("Error: ancestral_sequence is empty")
        #print(orthogroup_names)
        #print(ancestral_sequence)
 
    ### get seperate lists of species and genes
    #species = [item.rsplit("_", 1)[0] for item in ancestral_gene_ids]
    species = ["_".join(item.split("_")[:2]) for item in ancestral_gene_ids]
    #genes = [item.rsplit("_", 1)[1] for item in ancestral_gene_ids]
    genes = ["_".join(item.split("_")[2:]) for item in ancestral_gene_ids]
    all_species = set(species)
    all_species = list(all_species)
    
    ## open species fasta, and make a dictionary of gene to sequence
    all_fasta_dicts = MakeAllSpeciesFastaDict(all_species, ortho_folder_path)
    
    # build the ancestral genome
    ancestral_genome = BuildAncestralGenome(species, genes, orthogroup_names, all_fasta_dicts, focal_node)
    
    # save it as a fasta
    WriteAncestralFasta(focal_node, ancestral_genome, ortho_folder_path)
        
    #print("Built ancestral genome for node " + focal_node)

#######



########
def main(ortho_folder_path):
    
    #print("---Analysing", len(species_names)-1, "ancestral nodes", "from", len(species_names), "species" )
    
    ### Construct ancestral genome ###
    #print("---Traversing the species tree...")
    
    # make ancestral genomes folder, if not exist
    os.makedirs(os.path.join(ortho_folder_path,"AncestralGenomes/"), exist_ok=True)
    ## load orthogroup file
    HOG_file = []
    #HOG_filepath = os.path.join(ortho_folder_path, "Phylogenetic_Hierarchical_Orthogroups", f"{N}.tsv")
    HOG_filepath = os.path.join(ortho_folder_path, "Orthogroups", "Orthogroups.tsv")
    with open(HOG_filepath, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        reader.fieldnames = [fieldname.replace('.', '_') for fieldname in reader.fieldnames]
        for row in reader:
            HOG_file.append(row)

    ## pre-processing
    # remove orthogroups that have <4 genes
    OG_file = FilterHogs(HOG_file)

    processed_OG_file = []
    for record in OG_file:
        # Drop 2nd and 3rd columns
        if 'OG' in record: del record['OG']
        if 'Gene Tree Parent Clade' in record: del record['Gene Tree Parent Clade']
        # Rename the first column to 'Orthogroup'
        record['Orthogroup'] = record.pop('Orthogroup')
        # Add the modified record to the new list
        processed_OG_file.append(record)
    OG_file = processed_OG_file

    ## load gains file
    gains_file_path = os.path.join(ortho_folder_path, "GainsLossDuplication", "Gains.tsv")
    #gains = pd.read_csv(gains_file_path, sep='\t', names=['Gain Node', 'Parent Node', 'Orthogroup'])
    gains = []
    with open(gains_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['Gain Node', 'Parent Node', 'Orthogroup'], delimiter='\t')
        next(reader)
        for row in reader:
            gains.append(row)

    ## load and plot species tree
    species_tree_path = os.path.join(ortho_folder_path, f"Species_Tree/SpeciesTree_rooted_node_labels.txt")
    species_tree = ete3.Tree(species_tree_path, quoted_node_names=True, format=1)
    # record old species names
    species_namesO = species_tree.get_leaf_names()
    # Rename leaves by replacing '.' with '_'
    for leaf in species_tree.iter_leaves():
        leaf.name = leaf.name.replace('.', '_')
    # Save the modified tree
    species_tree.name = "N0"
    modified_tree_path = os.path.join(ortho_folder_path, f"Species_Tree/SpeciesTree_N0_rooted_altnamed_node_labels.txt")
    species_tree.write(outfile=modified_tree_path, format=1)
    # use the modified tree
    # define species names
    species_names = species_tree.get_leaf_names()
    
    # Traverse the tree in postorder
    node_list = []
    for node in species_tree.traverse("postorder"):
        # Check if the node is neither the root nor a leaf
        #if not node.is_leaf() and node.up is not None:
        if not node.is_leaf():
            AncestralGenome(node, species_tree, OG_file, ortho_folder_path, gains, species_names)
            node_list.append(node.name)
    
    ## make a dict of the ancestral fastas
    anc_dic, gene_counts, total_lengths, orthogroup_counts = MakeAncestralSpeciesFastaDict(node_list, os.path.join(ortho_folder_path, "AncestralGenomes"))
    anc_data = {
        "Node": list(gene_counts.keys()),
        "Gene_Count": list(gene_counts.values()),
        "Number_of_amino_acids": list(total_lengths.values())
    }
   
    anc_data_list = [
        {"Node": node, "Gene_Count": gene_count, "Number_of_amino_acids": num_amino_acids}
        for node, gene_count, num_amino_acids in zip(anc_data["Node"], anc_data["Gene_Count"], anc_data["Number_of_amino_acids"])
    ]
    # Sort the list of dictionaries by 'Node'
    sorted_anc_data = sorted(anc_data_list, key=lambda x: x['Node'])
    # Save as a txt file
    output_file_path = os.path.join(ortho_folder_path, "AncestralGenomes", "AncestralGenomes.txt")
    with open(output_file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=sorted_anc_data[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted_anc_data)
    
    # Convert orthogroup_counts dictionary to a list of dictionaries
    og_counts_list = [{'Orthogroup': key, **value} for key, value in orthogroup_counts.items()]
    
    # Fill missing values with 0
    species_names = {species for counts in orthogroup_counts.values() for species in counts}
    for og in og_counts_list:
        for species in species_names:
            if species not in og:
                og[species] = 0

    # Sort by 'Orthogroup'
    og_counts_list = sorted(og_counts_list, key=lambda x: x['Orthogroup'])
    # Get all column names for the CSV file
    all_columns = sorted({key for d in og_counts_list for key in d.keys()})
    # Write to CSV
    output_file_path = os.path.join(ortho_folder_path, "AncestralGenomes", "Ancestral_HOG_counts.csv")
    with open(output_file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=all_columns)
        writer.writeheader()
        for row in og_counts_list:
            writer.writerow(row)
    
    ## save species tree
    PlotTree(species_tree, os.path.join(ortho_folder_path, "AncestralGenomes"), "species_tree")
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Generate Gene Trees for Hierarchical Orthogroups")
    parser.add_argument('folder', type=str, help="Path to the Orthofinder results folder")
    parser.add_argument('node', type=str, help="Node of the species tree")

    args = parser.parse_args()
    ortho_folder_path = args.folder
    N = args.node
    
    # print("---------------------------------------------------")
    # print("------- Generating Gene Trees ---------------------")
    # print("---------------------------------------------------")
        
    # print("---Building gene trees. This may take a few minutes")
    
    main(ortho_folder_path, N)