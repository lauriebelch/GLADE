#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:34:16 2024

@author: user
"""

import os
import csv
import argparse
import numpy as np
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from common_functions import ListDictToTSV

#ortho_folder_path = "/Users/user/Dropbox/ORTHOFINDER/example_data/Proteomes16/OrthoFinder/Results_May24_2"
#N = "N0"

def NA_INT(value):
    return int(float(value)) if value != 'NA' else 'NA'

def process_data_item(d, ortho_dict, BS_dict):
    family_size = "NA"
    parent_size = "NA"
    Branch_length = "NA"
    
    match = ortho_dict.get(d['Orthogroup'])
    if match:
        family_size = match.get(d['Focal_Node'], "NA")
        if d['Parent_Node'] != "root":
            parent_size = match.get(d['Parent_Node'], "NA")
    match2 = BS_dict.get(d['Branch'])
    if match2:
        Branch_length = match2.get('branch_length', "NA")
    
    d['family_size'] = NA_INT(family_size)
    d['parent_size'] = NA_INT(parent_size)
    d['Branch_length'] = float(Branch_length) if Branch_length != "NA" else "NA"
    
    if d['family_size'] == "NA" or d['parent_size'] == "NA":
        d['change'] = 'NA'
        d['change_per_time'] = 'NA'
    else:
        d['change'] = d['family_size'] - d['parent_size']
        if d['Branch_length'] == 0 or d['Branch_length'] == "NA":
            d['change_per_time'] = 'NA'
        else:
            d['change_per_time'] = abs(d['change']) / d['Branch_length']
    
    return d

def main(ortho_folder_path, n_threads):

    ##### load some files
    ## load gains file
    gains_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Gains.tsv")
    #gains = pd.read_csv(gains_file_path, sep='\t', names=['Gain Node', 'Parent Node', 'Orthogroup'])
    Gains = []
    with open(gains_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['Gain Node', 'Parent Node', 'Orthogroup'], delimiter='\t')
        next(reader)
        for row in reader:
            Gains.append(row)
            
    loss_s_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Loss_speciation_bybranch.tsv")
    Loss_s = []
    with open(loss_s_file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        next(reader)
        for row in reader:
            Loss_s.append(row)

    Dupes_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Duplications_bybranch.tsv")
    Dupes = []
    with open(Dupes_file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        next(reader)
        for row in reader:
            Dupes.append(row)

    loss_pd_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Loss_postduplication_bybranch.tsv")
    Loss_pd = []
    with open(loss_pd_file_path, mode='r') as file:
        reader = csv.DictReader(file , delimiter='\t')
        next(reader)
        for row in reader:
            Loss_pd.append(row)
            
    Node_info_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "node_parent_info.txt")
    Node_info = []
    with open(Node_info_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['Node', 'Parent Node'], delimiter=',')
        #next(reader)
        for row in reader:
            Node_info.append(row)

    Ortho_ancest_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/AncestralGenomes", "Ancestral_HOG_counts.csv")
    Ortho_ancest = []
    with open(Ortho_ancest_file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter=',')
        next(reader)
        for row in reader:
            Ortho_ancest.append(row)
            
    BS_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Branch_statistics.tsv")
    BS = []
    with open(BS_file_path, mode='r') as file:
        reader = csv.DictReader(file, fieldnames=['branch', 'branch_length','N_gains', 'N_speciation_losses', 'N_duplications', 'N_postduplication_losses'], delimiter='\t')
        next(reader)
        for row in reader:
            BS.append(row)

    HOG_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/", "Orthogroups.tsv")
    HOG = []
    with open(HOG_file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        #next(reader)
        for row in reader:
            HOG.append(row)


    # loop through species, store family size of each OG in each species
    ## make variables
    Desc = []
    Family_ID = []
    for hog in HOG:
        Desc.append("OG_" + hog['Orthogroup'])
        Family_ID.append(hog['Orthogroup'])
    # list of species
    species = list(HOG[0].keys())[1:]
    # Create output list of dictionaries to store results
    output = [{"Desc": desc, "Family_ID": family_id} for desc, family_id in zip(Desc, Family_ID)]
    
    # do the loop
    mynum = 1  # start at 1, as this is the column where species start
    for SPECIES in species:
        # Initialize count_data list
        count_data = []
        # Smaller loop through hogs
        for hog in HOG:
            chosen_row = hog[list(hog.keys())[mynum]]
            count = len(chosen_row.split(",")) if chosen_row else 0
            count_data.append(count)
        # Add count_data to output
        for idx in range(len(output)):
            output[idx][SPECIES] = count_data[idx]
        mynum += 1  # Move to the next species column
    # save output
    ListDictToTSV(output, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/extant_OG_counts.tsv"))
    
    ## make lists of combinations for duplications/losses
    Duplications_list = []
    Loss_pd_list = []
    Loss_s_list = []
    All_branches = []
    All_orthogroup = []
    for dupe in Dupes:
        Duplications_list.append( dupe['Branch'] + " " + dupe['Orthogroups_duplicated'])
    for loss in Loss_pd:
        Loss_pd_list.append( loss['Branch'] + " " + loss['Orthogroups_lost_duplication'])
    for loss in Loss_s:
        Loss_s_list.append( loss['Branch'] + " " + loss['Orthogroups_lost_speciation'])
    for node in Node_info:
        All_branches.append( node['Parent Node'] + "___" + node['Node'] )
    for gain in Gains:
        All_orthogroup.append( gain['Orthogroup'] )
    All_orthogroup = list(set(All_orthogroup))
    
    # all orthogroup-branch combos
    combos1 = np.array([(x, y) for x in All_branches for y in All_orthogroup])
    combos1.shape
    combos = np.array([f"{y}_{x}" for x in All_branches for y in All_orthogroup])
    
    # Split the branches into Parent_Node and child
    Parent_Node = [x.split("___")[0] for x, y in combos1]
    child = [x.split("___")[1] for x, y in combos1]
    
    # Create the data dictionary
    data = []
    for i in range(len(combos)):
        orthogroup_branch = combos[i]
        branch, orthogroup = combos1[i]
        data.append({
            "Orthogroup_Branch": orthogroup_branch,
            "Orthogroup": orthogroup,
            "Branch": branch,
            "Parent_Node": Parent_Node[i],
            "Focal_Node": child[i]
        })
    
    # Convert output column to orthogroup
    for ortho in output:
        ortho['Orthogroup'] = ortho.pop('Family_ID')
    
    ## Initialize a defaultdict of dictionaries
    merged_dict = defaultdict(dict)
    # Merge ortho_ancest into merged_dict
    for item in Ortho_ancest:
        merged_dict[item['Orthogroup']].update(item)
    # Merge output into merged_dict, ensuring all keys from output are present
    for item in output:
        if 'Orthogroup' in item:
            merged_dict[item['Orthogroup']].update(item)
    # Gather all possible keys from both ortho_ancest and output
    all_keys = set()
    for item in Ortho_ancest:
        all_keys.update(item.keys())
    for item in output:
        all_keys.update(item.keys())
    # Convert merged_dict back to a list of dictionaries
    new_ortho_ancest = list(merged_dict.values())
    # Fill missing values with 0 for all possible keys
    for item in new_ortho_ancest:
        for key in all_keys:
            if key not in item:
                item[key] = 0
    # Ensure all keys from output are in new_ortho_ancest with full columns
    for item in output:
        if item['Orthogroup'] not in [d['Orthogroup'] for d in new_ortho_ancest]:
            new_item = {k: 0 for k in all_keys}  # Initialize new item with all keys and default 0
            new_item.update(item)  # Update with actual values from output
            new_ortho_ancest.append(new_item)
     
    # update Ortho_ancest
    Ortho_ancest = new_ortho_ancest
    
    # calculate size at branches
    ortho_dict = {oa['Orthogroup']: oa for oa in Ortho_ancest}
    BS_dict = {oa['branch']: oa for oa in BS}

    n_threads = max(1, min(n_threads, cpu_count()))
    with Pool(processes=n_threads) as pool:
        results = pool.starmap(process_data_item, [(d, ortho_dict, BS_dict) for d in data])
        
    # family_size = ["NA"] * len(data)
    # parent_size = ["NA"] * len(data)
    # Branch_length = ["NA"] * len(data)
    # counter = -1
    # for d in data:
    #     counter += 1
    #     # get matching og from ortho_ancest
    #     #print(d['Orthogroup'])
    #     match = ortho_dict.get(d['Orthogroup'])
    #     family_size[counter] = match[d['Focal_Node']]
    #     if d['Parent_Node'] != "root":
    #         parent_size[counter] = match[d['Parent_Node']]
    #     match2 = BS_dict.get(d['Branch'])
    #     Branch_length[counter] = match2['branch_length']
    ## end loop through data
    
    # checking if dict has all the keys we need
    #value = "N0.HOG0002441"
    #[entry for entry in Ortho_ancest if entry.get('Orthogroup') == value]
        
    # # Add variables to each dictionary in data
    # for i in range(len(data)):
    #     data[i]['family_size'] = NA_INT(family_size[i])
    #     data[i]['parent_size'] = NA_INT(parent_size[i])
    #     data[i]['Branch_length'] = float(Branch_length[i])
    
    # # Calculate change
    # for d in data:
    #     #d['change'] = d['family_size'] - d['parent_size']
    #     if d['family_size'] == "NA" or d['parent_size'] == "NA":
    #         d['change'] = 'NA'
    #     else:
    #         d['change'] = d['family_size'] - d['parent_size']
    #         #print( d['change'] )
    #     bl = d['Branch_length']
    #     if bl == 0 or bl == "NA" or d['change'] == "NA":
    #         d['change_per_time'] = 'NA'
    #     else:
    #         d['change_per_time'] = abs(d['change']) / bl
    #         #print( d['change_per_time'] )
            
        #print(d['change_per_time'])
        
    #data = [entry for entry in data if entry['change_per_time'] != 'NA']
    data_processed = []
    for d in results:
        if d['change_per_time'] != 'NA':
           data_processed.append(d)
    data = data_processed
    
    ListDictToTSV(data, os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication/OrthogroupBranchChange.tsv"))

