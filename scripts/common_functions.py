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

def FilterHogs(HOG_file, min_genes=4):
    filtered_HOG_file = []
    for row in HOG_file:
        total = 0
        for species, gene_list in row.items():
            if species == 'Orthogroup':
                continue
            if gene_list:
                total += len(gene_list.split(", "))
        if total >= min_genes:
            filtered_HOG_file.append(row)
    return filtered_HOG_file

## function to save tsv
def ListDictToTSV(data, filename, column_order=None):
    if not data:
        raise ValueError(f"No data to write: {filename}")
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        # header
        if column_order:
            writer.writerow(column_order)
        else:
            writer.writerow(list(data[0].keys()))
        # rows
        for row in data:
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
            # leaf has itself as both children
            leaf = node.name
            data.append({
                "Node": node.name,
                "Child1": [leaf],
                "Child2": [leaf]
            })
            continue
        child1 = node.children[0].get_leaf_names()
        child2 = node.children[1].get_leaf_names()
        data.append({
            "Node": node.name,
            "Child1": child1,
            "Child2": child2
        })
    return data

#s_df = SpeciesTreeToDataFrame(species_tree)
       
# function to map gene tree and species tree nodes
# find the species under a node in the gene tree
# and go down the species dataframe to find the node where.         
## .. we have representatives in both children
## needed to find losses after duplication
def MapGeneSpeciesTree(gene_tree, species_tree, species_names):
    s_df = SpeciesTreeToDataFrame(species_tree)
    gene_tree.name = "n0"
    mapping = []
    for node in gene_tree.traverse():
        if node.is_leaf():
            continue
        # species in this gene-tree node
        leaves = { leaf.name.split("_")[0] for leaf in node.get_leaves() }
        # find the species_tree node whose two children both intersect leaves
        for row in s_df:
            child1 = set(row["Child1"])
            child2 = set(row["Child2"])
            if leaves & child1 and leaves & child2:
                mapping.append({
                    "Gene Node": node.name,
                    "Species Node": row["Node"]
                })
                break
    return mapping
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
    duplications = []
    mapping = MapGeneSpeciesTree(gene_tree, species_tree, species_names)
    gene_tree.name = "n0"
    # Convert mapping list to dict for quick lookup
    mapdict = { m["Gene Node"]: m["Species Node"] for m in mapping }
    for node in gene_tree.traverse("postorder"):
        if node.is_leaf():
            continue
        # species in children
        leaves1 = { leaf.name.split("_")[0] for leaf in node.children[0].get_leaves() }
        leaves2 = { leaf.name.split("_")[0] for leaf in node.children[1].get_leaves() }
        # duplication if overlap
        shared = leaves1 & leaves2
        if not shared:
            continue
        # species node
        species_node = mapdict.get(node.name)
        # compute support
        st_node = species_tree.search_nodes(name=species_node)[0]
        expected_species = { leaf.name for leaf in st_node.get_leaves() }
        support = len(shared) / len(expected_species)
        # label event
        if support >= 0.5:
            node.add_features(duplication_label="D")
        else:
            node.add_features(low_support_duplication_label="lD")
        duplications.append({
            "genetree_node": node.name,
            "leaves1": node.children[0].get_leaf_names(),
            "leaves2": node.children[1].get_leaf_names(),
            "speciestree_node": species_node,
            "support": support
        })
    return duplications
#df1 = FindDuplications(example_gene_tree)              