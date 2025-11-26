#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:04:20 2024

@author: user
"""

## this script has functions that calculate gains, losses, and dupliations on branches

# import libraries
import ete3
import os
import csv
import sys
import argparse
from common_functions import ListDictToTSV

csv.field_size_limit(sys.maxsize)


## functions ###################################################################
# function to traverse species trees and get info on nodes and branches
def SpeciesTreeTraverse(species_tree):
    # first traversal to get leaf names under each node
    node_leaves = {}
    for node in species_tree.traverse("postorder"):
        # record leaves
        node_leaves[node.name] = {'leaves': [leaf.name for leaf in node.get_leaves()]}
        # Record child nodes
        if node.is_leaf():
            # Leaf nodes have no children
            node_leaves[node.name]['children'] = []
        else:
            # Non-leaf nodes: capture child node names
            node_leaves[node.name]['children'] = [child.name for child in node.children]
    
    ## second traversal to define branch names, and calculate branch lengths
    ## also use this traversal to get parent nodes
    node_parent_list = []
    BRANCHES = {} # a dictionary to store info on branches
    for node in species_tree.traverse("postorder"):
            # stuff for parent nodes
            parent_name = node.up.name if node.up else None
            parent_name = node.up.name if node.up else "root"
            node_parent_list.append((node.name, parent_name))
            #if not node.is_root():
            if node:
                # Naming the branch as 'parentnode__childnode'
                node_name = f"root___{node.name}" # special case
                if not node.is_root():
                    node_name = f"{node.up.name}___{node.name}"
                node.name = node_name
                # record branch length, and define other variables
                BRANCHES[node_name] = {
                    "branch_length": node.dist,
                    "N_gains": 0,
                    "Orthogroups_gained": [],
                    "N_speciation_losses": [],
                    "Orthogroups_lost_speciation": [],
                    "N_duplications": [],
                    "Orthogroups_duplicated": [],
                    "N_postduplication_losses": [],
                    "Orthogroups_lost_postduplication": []
                    }
    return node_leaves, node_parent_list, BRANCHES

# Find branch where postduplication loss is
def FindLossBranch(loss_list, node_leaves):
    # I have the Species Node. The loss is in the child that doesn't have the species
    for row in loss_list:
        # lost species
        lost_spp = row["Lost Species"]
        # species node
        sn = row["Species Node"]
        # child node
        children = node_leaves[sn]['children']
        # which child node DOESN'T have lost Species
        #print(node_leaves[children[0]])
        #print(row)
        #print("###")
        ch0 = node_leaves[children[0]]['leaves']
        ch0 = [s.replace('.', '_') for s in ch0]
        ch1 = node_leaves[children[1]]['leaves']
        ch1 = [s.replace('.', '_') for s in ch1]
        if lost_spp in ch0:
            loss_node = children[1]
        elif lost_spp in ch1:
            loss_node = children[0]
        else:
            loss_node = None
        row['Species Child Node'] = loss_node
        row['Branch_name'] = sn + "___" + loss_node if loss_node else None
    return loss_list

### code ##########################################

## remove after testing
#ortho_folder_path = "/Users/user/Dropbox/ORTHOFINDER/PROJECT_GAIN_LOSS/Feng2024/Results_Jan08"
#N="N0"
###

def main(ortho_folder_path, n_threads):

    ## species tree
    #species_tree_path = os.path.join(ortho_folder_path, f"Species_Tree/SpeciesTree_rooted_node_labels.txt")
    species_tree_path = os.path.join(
        ortho_folder_path, "WorkingDirectory", "GladeWD", "SpeciesTree_rooted_node_labels.txt"
    )
    species_tree = ete3.Tree(species_tree_path, quoted_node_names=True, format=1)
    species_tree.name = "N0"
    
    # traverse species tree and extract info
    node_leaves, node_parent_list, BRANCHES = SpeciesTreeTraverse(species_tree)
    
    ## load gains file
    gains_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Gains.tsv")
    #gains = pd.read_csv(gains_file_path, sep='\t', names=['Gain Node', 'Parent Node', 'Orthogroup'])
    Gains = []
    with open(gains_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['Gain Node', 'Parent Node', 'Orthogroup'], delimiter='\t')
        next(reader)
        for row in reader:
            Gains.append(row)
    
    loss_s_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Loss_speciation.tsv")
    Loss_s = []
    with open(loss_s_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['Orthogroup', 'Node','Species', 'Child Node'], delimiter='\t')
        next(reader)
        for row in reader:
            Loss_s.append(row)
    
    Dupes_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Duplications.tsv")
    Dupes = []
    with open(Dupes_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['genetree_node', 'leaves1','leaves2','speciestree_node','support','Orthogroup'], delimiter='\t')
        next(reader)
        for row in reader:
            Dupes.append(row)
    
    loss_pd_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Loss_postduplication.tsv")
    Loss_pd = []
    with open(loss_pd_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['Orthogroup', 'Focal Node','Lost Species', 'Child Node', 'Species Node'], delimiter='\t')
        next(reader)
        for row in reader:
            Loss_pd.append(row)
    
    # Only count high support dupes
    high_support_dupes = [dupe for dupe in Dupes if float(dupe['support']) >= 0.5]
    
    # Add branch names
    for gain in Gains:
        gain['Branch_name'] = gain['Parent Node'] + "___" + gain['Gain Node']
    for loss in Loss_s:
        loss['Branch_name'] = loss['Node'] + "___" + loss['Child Node']
    # Get parent node of duplication node
    node_parent_dict = {node[0]: node[1] for node in node_parent_list}
    for dupe in high_support_dupes:
        parent_node = node_parent_dict.get(dupe['speciestree_node'], None)
        dupe['speciestree_parentnode'] = parent_node
        dupe['Branch_name'] = parent_node + "___" + dupe['speciestree_node'] if parent_node else None
    Dupes = high_support_dupes
    
    # Find branch where postduplication loss is
    Loss_pd = FindLossBranch(Loss_pd, node_leaves)
    
    for key in BRANCHES:
        # Find gains that match that branch
        gains_mini = [gain for gain in Gains if gain['Branch_name'] == key]
        BRANCHES[key]["N_gains"] = len(gains_mini)
        BRANCHES[key]["Orthogroups_gained"] = [gain['Orthogroup'] for gain in gains_mini]
        
        # Find speciation losses that match
        loss_s_mini = [loss for loss in Loss_s if loss['Branch_name'] == key]
        BRANCHES[key]["N_speciation_losses"] = len(loss_s_mini)
        BRANCHES[key]["Orthogroups_lost_speciation"] = [loss['Orthogroup'] for loss in loss_s_mini]
        
        # Find duplications that match
        dupes_mini = [dupe for dupe in Dupes if dupe['Branch_name'] == key]
        BRANCHES[key]["N_duplications"] = len(dupes_mini)
        BRANCHES[key]["Orthogroups_duplicated"] = [dupe['Orthogroup'] for dupe in dupes_mini]
        
        # Find postduplication losses that match
        loss_pd_mini = [loss for loss in Loss_pd if loss['Branch_name'] == key]
        BRANCHES[key]["N_postduplication_losses"] = len(loss_pd_mini)
        BRANCHES[key]["Orthogroups_lost_postduplication"] = [loss['Orthogroup'] for loss in loss_pd_mini]
    
    branch_list = []
    for key, value in BRANCHES.items():
        branch_entry = value.copy()
        branch_entry['branch'] = key
        branch_list.append(branch_entry)    
    
    # Create a folder if it doesn't exist
    os.makedirs(os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/"), exist_ok=True)
    # Specify columns to keep
    columns_to_keep = ['branch', 'branch_length', 'N_gains', 'N_speciation_losses', 'N_duplications', 'N_postduplication_losses']
    # Filter the branch_list to keep only the specified columns
    filtered_list = []
    for entry in branch_list:
        filtered_entry = {col: entry[col] for col in columns_to_keep if col in entry}
        filtered_list.append(filtered_entry)
    # Specify the file name
    file_name = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Branch_statistics.tsv")
    # Write to TSV file
    with open(file_name, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=columns_to_keep, delimiter='\t')
        writer.writeheader()
        for entry in filtered_list:
            writer.writerow(entry)
       
    ### save orthogroups gained/lost/duplicated
    
    # Create and save a CSV for each branch
    branch_keys = list(BRANCHES.keys())
    orthogroups_gained = [BRANCHES[branch]['Orthogroups_gained'] for branch in branch_keys]
    orthogroups_lost_speciation = [BRANCHES[branch]['Orthogroups_lost_speciation'] for branch in branch_keys]
    orthogroups_duplicated = [BRANCHES[branch]['Orthogroups_duplicated'] for branch in branch_keys]
    orthogroups_lost_postduplication = [BRANCHES[branch]['Orthogroups_lost_postduplication'] for branch in branch_keys]
    
    # Collect data into a list of dictionaries
    data = []
    for i, branch in enumerate(branch_keys):
        data.append({
            'Branch': branch,
            'Orthogroups_gained': orthogroups_gained[i],
            'Orthogroups_lost_speciation': orthogroups_lost_speciation[i],
            'Orthogroups_duplicated': orthogroups_duplicated[i],
            'Orthogroups_lost_duplication': orthogroups_lost_postduplication[i]
        })
    
    # Create and save exploded data for each category
    def explode_and_filter(data, key, subkey):
        exploded_data = []
        for entry in data:
            for item in entry[subkey]:
                if item:  # This will skip None or empty items
                    exploded_data.append({
                        'Branch': entry['Branch'],
                        subkey: item
                    })
        return exploded_data
    
    # Explode the lists and filter out None values
    exploded_gained = explode_and_filter(data, 'Branch', 'Orthogroups_gained')
    exploded_lost_speciation = explode_and_filter(data, 'Branch', 'Orthogroups_lost_speciation')
    exploded_duplicated = explode_and_filter(data, 'Branch', 'Orthogroups_duplicated')
    exploded_lost_postduplication = explode_and_filter(data, 'Branch', 'Orthogroups_lost_duplication')
    
    # Save each exploded data to separate CSV files
    columns_map = {
        'Orthogroups_gained': exploded_gained,
        'Orthogroups_lost_speciation': exploded_lost_speciation,
        'Orthogroups_duplicated': exploded_duplicated,
        'Orthogroups_lost_duplication': exploded_lost_postduplication
    }
    
    #Save each exploded data to separate TSV files
    ListDictToTSV(exploded_gained, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Gains_bybranch.tsv"), column_order=['Branch', 'Orthogroups_gained'])
    ListDictToTSV(exploded_lost_speciation, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Loss_speciation_bybranch.tsv"), column_order=['Branch', 'Orthogroups_lost_speciation'])
    ListDictToTSV(exploded_duplicated, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Duplications_bybranch.tsv"), column_order=['Branch', 'Orthogroups_duplicated'])
    ListDictToTSV(exploded_lost_postduplication, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/Loss_postduplication_bybranch.tsv"), column_order=['Branch', 'Orthogroups_lost_duplication'])

########
