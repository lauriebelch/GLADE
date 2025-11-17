#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:07:56 2024

@author: user
"""

import csv
import ete3
import os
import re
import time



### functions for filtering and pre-processing#################################
## filter to get rid of hierarchical orthogroups with < 4 genes
'''
def FilterHogs(HOG_file):
    # Gene Tree Parent Clade is - for OG with <4 genes
    HOG_file = [row for row in HOG_file if 'n' in row['Gene Tree Parent Clade']]
    # some tiny HOGs remain, we need to count genes and remove them
    filtered_HOG_file = []
    for row in HOG_file:
        genes = ', '.join([value for key, value in row.items() if key not in {'Gene Tree Parent Clade', 'HOG'} and value]).split(', ')
        if len(genes) >= 4:
            filtered_HOG_file.append(row)
    return filtered_HOG_file
'''

def FilterHogs(HOG_file, min_genes=4):
    filtered_HOG_file = []
    for row in HOG_file:
        # Count all genes across species columns
        total_genes = 0
        for species, gene_list in row.items():
            if species == 'Orthogroup':
                continue
            if gene_list:
                total_genes += len(gene_list.split(', '))
        if total_genes >= min_genes:
            filtered_HOG_file.append(row)
    return filtered_HOG_file

## function to save tsv
def ListDictToTSV(data, filename, column_order=None):
    if not data:
        raise ValueError(f"The data list for '{os.path.basename(filename)}' is empty. Cannot write to TSV file.")
    with open(filename, mode='w', newline='', encoding='utf-8') as file:
       writer = csv.writer(file, delimiter='\t')
       # Determine the header order
       if column_order:
           header = column_order
       else:
           header = data[0].keys()
       writer.writerow(header)
       # Write the data rows
       for row in data:
           # Write the row values in the specified order
            if column_order:
                writer.writerow([row[col] for col in column_order])
            else:
                writer.writerow(row.values())
                
### functions for processing species tree, and mapping to gene tree  ##########
## function to make a dataframe showing leaves under each child of a node in species tree
## needed to identify gene losses at speciation
def SpeciesTreeToDataFrame(species_tree):
    data = []
    for node in species_tree.traverse():
        if node.is_leaf():
            leafnames = [sub.replace('.', '_') for sub in node.get_leaf_names()]
            data.append({
                "Node": node.name,
                "Child1": leafnames,
                "Child2": leafnames})
            continue
        leaves1 = node.children[0].get_leaf_names()
        leaves2 = node.children[1].get_leaf_names()
        # replace . with _
        # needed because species tree preserves fasta file names
        # whereas gene trees force _
        leaves1 = [sub.replace('.', '_') for sub in leaves1]
        leaves2 = [sub.replace('.', '_') for sub in leaves2]
        data.append({
            "Node": node.name,
            "Child1": leaves1,
            "Child2": leaves2})
    return data
#s_df = SpeciesTreeToDataFrame(species_tree)
       
# function to map gene tree and species tree nodes
# find the species under a node in the gene tree
# and go down the species dataframe to find the node where.         
## .. we have representatives in both children
## needed to find losses after duplication
def MapGeneSpeciesTree(gene_tree, species_tree, species_names):
    # make species dataframe
    s_df = SpeciesTreeToDataFrame(species_tree)
    # make dataframe to store map between gene and species tree
    data = []
    # make sure root of gene tree is named
    gene_tree.name = "n0"
    # traverse gene tree
    for node in gene_tree.traverse():
        if node.is_leaf(): continue # skip if node is leaf
        # if node is not a leaf, get leaves
        #leaves = set(re.sub(r"_[^_]*$", "", leaf.name) for leaf in node.get_leaves())  
        #leaves = set("_".join(leaf.name.split("_")[:2]) for leaf in node.get_leaves())
        # Find the species name that matches the beginning of leaf.name
        leaves = set()
        for leaf in node.get_leaves():
            # Find the species name that matches the beginning of leaf.name
            matching_species = next((species for species in species_names if leaf.name.startswith(species + "_")), None)
            if matching_species:
                leaves.add(matching_species)
            else:
                raise ValueError(f"No matching species found for leaf name '{leaf.name}'. Check species_names list.")

        for row in s_df:
            has_leaves_in_child1 = bool(set(row['Child1']).intersection(leaves))
            has_leaves_in_child2 = bool(set(row['Child2']).intersection(leaves))
            if has_leaves_in_child1 and has_leaves_in_child2:
                data.append({
                    "Gene Node": node.name,
                    "Species Node":row['Node']})
                break
    return data
#mdf = MapGeneSpeciesTree(example_gene_tree, species_tree)

# Function to extract gene names from species_gene format
def extract_species_name(leaf_name, species_names):
    # Find the matching species name at the start of leaf_name
    matching_species = next((species for species in species_names if leaf_name.startswith(species + "_")), None)
    # If a matching species is found, extract the gene name
    if matching_species:
        return matching_species
    else:
        return None  # or handle as appropriate if no match is found

## function to find duplications, label them, and record them
def FindDuplications(gene_tree, species_tree, species_names):
    # gene tree is a gene tree
    duplications = []
    # map gene tree to species tree
    mdf = MapGeneSpeciesTree(gene_tree, species_tree, species_names)
    # traverse the tree
    # make sure root of gene tree is named
    gene_tree.name = "n0"
    for node in gene_tree.traverse('postorder'):
        if node.is_leaf(): continue # skip if node is leaf
        # get leaf names for child clades 1 and 2
        #leaves1 = set(re.sub(r"_[^_]*$", "", leaf.name) for leaf in node.children[0].get_leaves())
        #leaves2 = set(re.sub(r"_[^_]*$", "", leaf.name) for leaf in node.children[1].get_leaves())
        #leaves1 = set("_".join(leaf.name.split("_")[:2]) for leaf in node.children[0].get_leaves())
        #leaves2 = set("_".join(leaf.name.split("_")[:2]) for leaf in node.children[1].get_leaves())
        leaves1 = set(extract_species_name(leaf.name, species_names) for leaf in node.children[0].get_leaves())
        leaves2 = set(extract_species_name(leaf.name, species_names) for leaf in node.children[1].get_leaves())
        #print(leaves1)
        # if any intersections, its a duplication
        if leaves1.intersection(leaves2):
            # define species node
            species_node = None
            for record in mdf:
                if record['Gene Node'] == node.name:
                    species_node = record['Species Node']
                    break
                    #print(mdf)
            # define species node in treenode form
            #print(species_node)
            #print(species_tree.search_nodes(name=species_node))
            species_nodeT = species_tree.search_nodes(name=species_node)[0]
            # define expected species
            expected_species = {leaf.name for leaf in species_nodeT.get_leaves()}
            # calculate support
            shared_species = leaves1.intersection(leaves2)
            support = len(shared_species) / len(expected_species)
            # add duplication label if support is high
            # otherwise, add low support duplication label
            if support >= 0.5:
                node.add_features(duplication_label="D")
            else:
                node.add_features(low_support_duplication_label="lD")
            duplications.append({
                "genetree_node": node.name,  # Make sure nodes have names or identifiers
                "leaves1": node.children[0].get_leaf_names(),  # Convert set to list for better compatibility with DataFrame
                "leaves2": node.children[1].get_leaf_names(),
                "speciestree_node": species_node,
                "support": support
            })
    df = duplications
    return df
#df1 = FindDuplications(example_gene_tree)              