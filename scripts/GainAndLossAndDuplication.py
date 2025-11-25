#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 08:33:46 2024

@author: OrthoFinder Laurie
"""

## this script has functions that calculate gains, losses, and dupliations

# import libraries
import ete3
import os
import re
import csv
import argparse
from multiprocessing import Pool, cpu_count
from common_functions import ListDictToTSV, SpeciesTreeToDataFrame, MapGeneSpeciesTree, FindDuplications, FilterHogs

#### FUNCTIONS ####
##############################################################################    
### functions to find gains, losses, and duplications
## function to check if OG has a member in both children of a node
# needed for FindGainNode
# node is the node from a tree traversal
# OG is the Orthogroups.tsv file
# row is the row of that file (when we loop through orthogroups)
def CheckChildMembers(node, row):
    # get leaf names for both children
    leaves1 = node.children[0].get_leaf_names()
    leaves2 = node.children[1].get_leaf_names()
    # get OG_file lists to look in
    og1 = any(leaf in row and row[leaf] not in [None, ''] for leaf in leaves1)
    og2 = any(leaf in row and row[leaf] not in [None, ''] for leaf in leaves2)
    return og1 and og2

## function to find the node where a gain happened
# species_tree is the species tree
# OG is the Orthogroups.tsv file
# row is the row of that file (when we loop through orthogroups)
def FindGainNode(species_tree, row):
    # traverse the tree
    # if OG has members in both children of node0, gain was at N0
    # otherwise find point where OG has members in both child nodes
    # if we get to a leaf and find it, the gain was at that leaf
    # gain happens at branch between node and its parent, so record parent too
    for node in species_tree.traverse():
        # if node is a leaf and it has species, thats the gain node
        if node.is_leaf():
            if node.name in row and row[node.name] not in [None, '']:
                return node.name, node.up.name
            continue
        # otherwise, use my checkchild function to find gain node
        if CheckChildMembers(node, row):
            return (node.name, node.up.name if node.up else "root")
    return None

# function to identify a gene loss based on being in one clade but not other
## needed for FindLossesAfterDuplications
def IdentifyGeneLoss(leaves1, leaves2, combined_children):
    # loss is if a species is in one clade but not another 
   loss_to_child1 = []
   loss_to_child2 = []
   # Set operations to find differences
   set_leaves1 = set(leaves1)
   set_leaves2 = set(leaves2)
   for child in combined_children:
       if child in set_leaves1 and child not in set_leaves2:
           loss_to_child2.append(child)
       elif child in set_leaves2 and child not in set_leaves1:
           loss_to_child1.append(child)
   return {
       'loss_to_child1': loss_to_child1,
       'loss_to_child2': loss_to_child2
   }
    
# function to find loss after duplication
# for a gene tree node...
# find it in mdf, and get the corresponding species node
# make s_df
# if its a duplication node
# find the species node in s_df to define the species we expect
# define the two clades under our node
# if any of the species are in one clade but not the other...
# there is a lost in the branch from node to the child node it is lost in
# if a species is not there at all
def FindLossesAfterDuplications(gene_tree, species_tree, species_names):
    s_df = SpeciesTreeToDataFrame(species_tree)
    mdf = MapGeneSpeciesTree(gene_tree, species_tree, species_tree.get_leaf_names())

    loss_records = []

    # map GeneNode -> SpeciesNode for quick lookup
    mapdict = {m["Gene Node"]: m["Species Node"] for m in mdf}

    for node in gene_tree.traverse():
        if node.is_leaf():
            continue
        if not hasattr(node, "duplication_label") and not hasattr(node, "low_support_duplication_label"):
            continue

        # species under each child
        leaves1 = {leaf.name.split("_")[0] for leaf in node.children[0].get_leaves()}
        leaves2 = {leaf.name.split("_")[0] for leaf in node.children[1].get_leaves()}

        # species node this gene-tree node maps to
        species_node = mapdict[node.name]

        # find this species_node in species dataframe
        for record in s_df:
            if record["Node"] == species_node:
                expected_species = record["Child1"] + record["Child2"]
                break

        loss_info = IdentifyGeneLoss(leaves1, leaves2, expected_species)

        # LOST IN CHILD1
        for lost in loss_info["loss_to_child1"]:
            loss_records.append({
                "Focal Node": node.name,
                "Lost Species": lost,
                "Child Node": node.children[0].name,
                "Species Node": species_node
            })

        # LOST IN CHILD2
        for lost in loss_info["loss_to_child2"]:
            loss_records.append({
                "Focal Node": node.name,
                "Lost Species": lost,
                "Child Node": node.children[1].name,
                "Species Node": species_node
            })

    return loss_records
#loss_df = FindLossesAfterDuplications(example_gene_tree, species_tree)

########################################################3
## 2) find duplications and losses after duplications ##
# to save duplications
# function to do this in parallel
def FindDuplicationsParallel(index, row, gain_nodes, ortho_folder_path, species_tree, gene_trees, species_names):
    duplications = []
    # to save losses after duplications
    postduplication_loss = []
    #gene_tree_file = list(gain_nodes.keys())[index] + "_tree.txt"
    #gene_tree = ete3.Tree(os.path.join(ortho_folder_path, "HOG_Gene_Trees/", gene_tree_file), quoted_node_names=True, format=1)
    orthogroup_name = list(gain_nodes.keys())[index]
    gene_tree= ete3.Tree(gene_trees[orthogroup_name], quoted_node_names=True, format=1)
    dupe_list = FindDuplications(gene_tree, species_tree, species_names)
    for dupe in dupe_list:
        dupe['Orthogroup'] = row['Orthogroup']
        duplications.append(dupe)
    # only do FindLosses if there is a duplication
    if dupe_list:
        loss_list = FindLossesAfterDuplications(gene_tree, species_tree, species_names)
        for loss in loss_list:
            loss['Orthogroup'] = row['Orthogroup']
            postduplication_loss.append(loss)
    return duplications, postduplication_loss

## function to find loss on species tree
# start at gain node of gene
# look for cases where gene is in one child but not its sister
# this means the loss is between node and child which lacks the species
def FindLossNode(species_tree, row, gain_nodes):
    ## a gene is lost on a branch leading to a clade if a gene has ..
    ## no homologs within that clade, but homologs in sister
    # find the gain node, which is where I want to start
    orthogroup = row['Orthogroup']
    gnn = gain_nodes.get(orthogroup, {'Gain Node': 'NA'})['Gain Node']
    start_node = species_tree.search_nodes(name = gnn)[0]
    # if start_node is a leaf, we don't need to do anything
    loss_events = []
    counter = 0
    # traverse the species tree, starting at gain node
    for node in start_node.traverse():
        counter += 1
        if counter ==1 & node.is_leaf(): break
        if node.is_leaf(): continue # skip if leaf
        # get leaf names for both children
        leaves1 = node.children[0].get_leaf_names()
        leaves2 = node.children[1].get_leaf_names()
        # get OG_file lists to look in
        og1 = any(leaf in row and row[leaf] not in [None, ''] for leaf in leaves1)
        og2 = any(leaf in row and row[leaf] not in [None, ''] for leaf in leaves2)
        # a loss is if a gene has rep in og1 but not og2 (or vice versa)
        # meaning og1 and og2 are one true and one false
        # the false is the one with the loss
        if og1 != og2:
            lost_clade_leaves = leaves2 if og1 else leaves1
            child_node_name = node.children[1].name if og1 else node.children[0].name
            loss_events.append({
                'Node': node.name,
                'Lost species': lost_clade_leaves,
                'Child node': child_node_name})
    return loss_events     

## 3) find losses after speciation
# function to do this in parallel
def FindSpeciationLossParallel(index, row, ortho_folder_path, species_tree, OG_file, gain_nodes):
    loss_nodes = []
    loss_node = FindLossNode(species_tree, row, gain_nodes)
    for event in loss_node:
       event['Orthogroup'] = row['Orthogroup']
    loss_nodes.extend(loss_node)
    return loss_nodes

# parallel processing for dupes
def ParaDupes(ortho_folder_path, species_tree, OG_file, gain_nodes, gene_trees, species_names):
    all_duplications = []
    all_postduplication_loss = []
    with Pool(processes=cpu_count()-2) as pool:
        results = pool.starmap(FindDuplicationsParallel, [(index, row, gain_nodes, ortho_folder_path, species_tree, gene_trees, species_names) for index, row in enumerate(OG_file)])
    # Aggregate results
    for duplications, postduplication_loss in results:
        all_duplications.extend(duplications)
        all_postduplication_loss.extend(postduplication_loss)
    return all_duplications, all_postduplication_loss

def ParaLoss(ortho_folder_path, species_tree, OG_file, gain_nodes):
    with Pool(processes=cpu_count()-2) as pool:
        speciation_loss_results = pool.starmap(FindSpeciationLossParallel, [(index, row, ortho_folder_path, species_tree, OG_file, gain_nodes) for index, row in enumerate(OG_file)])
    # Aggregate results
    all_loss_nodes = [item for sublist in speciation_loss_results for item in sublist]
    # Rename columns
    new_column_names = {
        'Node': 'Node',
        'Lost species': 'Species',
        'Child node': 'Child Node',
        'Orthogroup': 'Orthogroup'
    }
    for event in all_loss_nodes:
        for old_name, new_name in new_column_names.items():
            if old_name in event:
                event[new_name] = event.pop(old_name)
    # Reorder columns
    new_order = ['Orthogroup', 'Node', 'Species', 'Child Node']
    speciation_loss = [{col: event[col] for col in new_order if col in event} for event in all_loss_nodes]
    return speciation_loss    

     
########################
def main(ortho_folder_path):

    #######################################################
    ## load required files ##
    
    # Load all gene trees once into a dictionary
    gene_trees = {}
    with open(os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/Resolved_Gene_Trees.txt"), 'r') as f:
        for line in f:
            if ':' in line:
                orthogroup, newick_str = line.split(":", 1)
                gene_trees[orthogroup.strip()] = newick_str.strip()

    
    ## load orthogroup file
    HOG_file = []
    HOG_filepath = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD", "Orthogroups.tsv")
    with open(HOG_filepath, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        reader.fieldnames = [fieldname.replace('.', '_') for fieldname in reader.fieldnames]
        for row in reader:
            HOG_file.append(row)
            
    ## load and plot species tree
    species_tree = ete3.Tree(os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/SpeciesTree_rooted_node_labels.txt"), quoted_node_names=True, format=1)
    species_tree.name = "N0"

    for node in species_tree.traverse():
        if node.is_leaf():
            node.name = node.name.replace('.', '_')
    species_names = species_tree.get_leaf_names()
    species_names.sort(key=len, reverse=True)

    ## save node-parent details
    node_dict = {}

    node_dict = {}
    for node in species_tree.traverse("preorder"):
        parent = node.up
        node_dict[node.name] = parent.name if parent else "root"
        
    parent_output_file = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "node_parent_info.txt")
    os.makedirs(os.path.dirname(parent_output_file), exist_ok=True)
    with open(parent_output_file, "w") as file:
        for node, parent in node_dict.items():
            file.write(f"{node},{parent}\n")
    
    #######################################################
    ## pre-processing
    # remove orthogroups that have <4 genes
    OG_file = FilterHogs(HOG_file)
    
    processed_OG_file = []
    for record in OG_file:
        # Drop 2nd and 3rd columns
        if 'OG' in record: del record['OG']
        if 'Gene Tree Parent Clade' in record: del record['Gene Tree Parent Clade']
        # Rename the first column to 'Orthogroup'
        #record['Orthogroup'] = record.pop('Orthogroup')
        # Add the modified record to the new list
        processed_OG_file.append(record)
    OG_file = processed_OG_file
    
    gain_nodes = {}
    for row in OG_file:
         gain_node, parent_node = FindGainNode(species_tree, row)
         gain_nodes[row['Orthogroup']] = {'Gain Node':gain_node, 'Parent Node': parent_node}
    gains = gain_nodes
     
    gains_list = []
    for key, value in gains.items():
         entry = {'Orthogroup': key}
         entry.update(value)
         gains_list.append(entry)
    
    duplications, postduplication_loss = ParaDupes(ortho_folder_path, species_tree, OG_file, gain_nodes, gene_trees, species_names)
    speciation_loss = ParaLoss(ortho_folder_path, species_tree, OG_file, gain_nodes)
    
    ## make ancestral genomes folder, if not exist
    os.makedirs(os.path.join(ortho_folder_path,"WorkingDirectory/GladeWD/GainsLossDuplication/"), exist_ok=True)
    
    ListDictToTSV(duplications, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Duplications.tsv"))
    ListDictToTSV(speciation_loss, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Loss_speciation.tsv"))
    
    column_order = ['Gain Node', 'Parent Node', 'Orthogroup']
    ListDictToTSV(gains_list, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Gains.tsv"), column_order)
    column_order = ['Orthogroup', 'Focal Node', 'Lost Species', 'Child Node', 'Species Node']
    #print(postduplication_loss)
    ListDictToTSV(postduplication_loss, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Loss_postduplication.tsv"), column_order)
    

########
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Generate Gene Trees for Hierarchical Orthogroups")
    parser.add_argument('folder', type=str, help="Path to the Orthofinder results folder")

    args = parser.parse_args()
    ortho_folder_path = args.folder
    
    # print("---------------------------------------------------")
    # print("------- Generating Gene Trees ---------------------")
    # print("---------------------------------------------------")
        
    # print("---Building gene trees. This may take a few minutes")
    
    main(ortho_folder_path)
       
    
