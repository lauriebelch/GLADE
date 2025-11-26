#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reconstruct ancestral genomes using numeric species and gene IDs.

Assumptions in NUMERIC MODE:
- Species tree leaves are numeric strings: "0", "1", "2", ...
- Gene tree leaves are "<speciesCode>_<geneCode>", e.g. "1_198"
- Orthogroups.tsv, gains, duplications, etc. are already numeric.
- Per-species FASTAs are in Orthofinder WorkingDirectory as:
    WorkingDirectory/Species0.fa, Species1.fa, ...
  with headers like: >0_3, >1_27, ...
"""

import os
import ete3
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import csv
from common_functions import FilterHogs, FindDuplications

import sys
csv.field_size_limit(sys.maxsize)


# SPECIES TREE (for summary)
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
    t.render(os.path.join(path, filename + ".png"), w=1920, units="px", tree_style=ts)


# Load a single gene tree for a given orthogroup from the BIG file
def load_gene_tree_from_big_file(og, tree_file_path):
    with open(tree_file_path, "r") as f:
        for line in f:
            if line.startswith(og + ": "):
                _, tree_str = line.strip().split(":", 1)
                return ete3.Tree(tree_str.strip(), quoted_node_names=True, format=1)
    raise ValueError(f"OG {og} not found in gene tree file {tree_file_path}.")


# For a given orthogroup + species-tree node, pick representative leaf genes
# for that ancestral genome (numeric version).
def GetAncestralGenes(OG, node, species_tree, ortho_folder_path, species_names):
    og = OG
    all_selected_sequences = []

    # Load numeric gene tree from GladeWD
    tree_file_path = os.path.join(
        ortho_folder_path, "WorkingDirectory", "GladeWD", "Resolved_Gene_Trees.txt"
    )
    gene_tree = load_gene_tree_from_big_file(og, tree_file_path)

    # Ensure the gene tree root is named "n0" for distance calculations
    gene_tree.name = "n0"

    # Species below the focal node in the species tree
    target_node = species_tree.search_nodes(name=node.name)[0]
    target_species = {leaf.name for leaf in target_node.get_leaves()}  # numeric codes: "0","1",...

    # All leaves in gene tree (node objects + names)
    all_leaves_nodes = gene_tree.get_leaves()
    all_leaves_names = np.array([leaf.name for leaf in all_leaves_nodes])

    # Extract species codes from gene tree leaves: "<species>_<geneCode>"
    all_species = np.array([name.split("_")[0] for name in all_leaves_names])

    # Keep only genes whose species is in the target clade
    idx = np.where(np.isin(all_species, list(target_species)))[0]
    if len(idx) == 0:
        # No genes from this OG in this clade → no contribution to this ancestor
        return None

    leaves_to_stay = [all_leaves_nodes[i] for i in idx]

    # Identify duplications on the full gene_tree before pruning
    duplications = FindDuplications(gene_tree, species_tree, species_names)
    duplications = [d for d in duplications if d["support"] >= 0.5]

    # Prune down to leaves in the target clade
    gene_tree.prune(leaves_to_stay, preserve_branch_length=True)

    # Handle duplications after the target node:
    # randomly drop one child clade (leaves1 or leaves2) for those dupes
    # whose speciestree_node is a descendant of target_node.
    desc_nodes = {n.name for n in target_node.get_descendants()}
    dupes_to_go = set()

    if duplications:
        # Indices of duplications that occur below the target node
        idx_after = [i for i, d in enumerate(duplications) if d["speciestree_node"] in desc_nodes]
        choice_cols = np.random.choice(["leaves1", "leaves2"], size=len(idx_after))
        for dup_idx, col in zip(idx_after, choice_cols):
            for leaf_name in duplications[dup_idx][col]:
                dupes_to_go.add(leaf_name)

    # Remove chosen duplicate leaves
    all_leaves_nodes = gene_tree.get_leaves()
    all_leaves_names = np.array([leaf.name for leaf in all_leaves_nodes])
    keep_mask = ~np.isin(all_leaves_names, list(dupes_to_go))
    keep_leaves = [all_leaves_nodes[i] for i in np.where(keep_mask)[0]]
    gene_tree.prune(keep_leaves, preserve_branch_length=True)

    # Handle duplications before the target node:
    # count how many duplications on the path from the root to target_node.
    # Each duplication adds 1 expected copy.
    ancestors = [target_node.name]
    up_node = target_node.up
    while up_node:
        ancestors.append(up_node.name)
        up_node = up_node.up

    expected_copies = 1
    if duplications:
        dupes_anc = [d for d in duplications if d["speciestree_node"] in ancestors]
        expected_copies += len(dupes_anc)

    # Choose representative genes:
    # compute root to leaf distances
    # iteratively choose medians until expected_copies or leaf list is exhausted
    leaf_nodes = gene_tree.get_leaves()
    leaf_names = np.array([leaf.name for leaf in leaf_nodes])

    distances = [leaf.get_distance("n0") for leaf in leaf_nodes]
    distances = np.array(distances)

    selected_sequences = []
    cur_dist = distances.copy()
    cur_leaves = leaf_names.copy()

    for _ in range(expected_copies):
        if len(cur_dist) == 0:
            break
        med_idx = np.argpartition(cur_dist, len(cur_dist) // 2)[len(cur_dist) // 2]
        selected_sequences.append(cur_leaves[med_idx])
        cur_dist = np.delete(cur_dist, med_idx)
        cur_leaves = np.delete(cur_leaves, med_idx)

    all_selected_sequences.append(selected_sequences)
    return all_selected_sequences

# Read numeric species FASTAs from Orthofinder WorkingDirectory
def MakeAllSpeciesFastaDict(all_species, ortho_folder_path):
    """
    all_species = iterable of numeric species codes, e.g. ["0","1","2"]
    FASTA files expected at:
        <ortho_folder_path>/WorkingDirectory/Species{code}.fa
    with headers like:
        >0_3
        >0_4
    """
    all_fasta_dicts = {}
    wd = os.path.join(ortho_folder_path, "WorkingDirectory")

    for species_code in all_species:
        fasta_path = os.path.join(wd, f"Species{species_code}.fa")
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"Expected file not found: {fasta_path}")

        with open(fasta_path, "r") as f:
            fasta_data = f.read()

        entries = fasta_data.split(">")
        fasta_dict = {}
        for entry in entries:
            if entry.strip():
                lines = entry.splitlines()
                gene_name = lines[0].split()[0]   # e.g. "0_3"
                seq = "".join(lines[1:])
                fasta_dict[gene_name] = seq

        all_fasta_dicts[species_code] = fasta_dict

    return all_fasta_dicts

# Read ancestral FASTAs and gather stats
def MakeAncestralSpeciesFastaDict(all_nodes, fasta_folder):
    """
    all_nodes = list of ancestral node names, e.g. ["N1","N2",...]
    fasta_folder = path to AncestralGenomes folder
    Reads <fasta_folder>/<node>.fasta for each node.
    """
    all_fasta_dicts = {}
    gene_counts = {}
    total_lengths = {}
    orthogroup_counts = {}

    for node_name in all_nodes:
        species_partial_path = f"{node_name}.fasta"
        species_path = os.path.join(fasta_folder, species_partial_path)
        if not os.path.exists(species_path):
            continue

        with open(species_path, "r") as f:
            fasta_data = f.read()

        entries = fasta_data.split(">")
        fasta_dict = {}
        total_len = 0
        gene_count = 0

        for entry in entries:
            if entry.strip():
                lines = entry.splitlines()
                gene_name = lines[0].split()[0]
                seq = "".join(lines[1:])
                fasta_dict[gene_name] = seq
                total_len += len(seq)
                gene_count += 1

                # gene_name format:
                # >{focal_node}_{orthogroup}_{gene_num}_{gene_num2}_{gene_id}
                parts = gene_name.split("_")
                if len(parts) < 2:
                    continue
                orthogroup = parts[1]

                if orthogroup not in orthogroup_counts:
                    orthogroup_counts[orthogroup] = {}
                if node_name not in orthogroup_counts[orthogroup]:
                    orthogroup_counts[orthogroup][node_name] = 0
                orthogroup_counts[orthogroup][node_name] += 1

        all_fasta_dicts[node_name] = fasta_dict
        gene_counts[node_name] = gene_count
        total_lengths[node_name] = total_len

    return all_fasta_dicts, gene_counts, total_lengths, orthogroup_counts


# Build a single ancestral genome FASTA text from selected genes
def BuildAncestralGenome(species_list, gene_list, orthogroup_list, fasta_dicts, focal_node):
    fasta_content = ""

    # count genes per orthogroup
    og_counts = {}
    for og in orthogroup_list:
        og_counts[og] = og_counts.get(og, 0) + 1

    current_gene_num = {}

    for species_id, gene_id, orthogroup in zip(species_list, gene_list, orthogroup_list):
        current_gene_num[orthogroup] = current_gene_num.get(orthogroup, 0) + 1
        gene_num = current_gene_num[orthogroup]
        gene_num2 = og_counts[orthogroup]

        if species_id not in fasta_dicts:
            raise KeyError(f"Species ID '{species_id}' not found in fasta_dicts.")
        if gene_id not in fasta_dicts[species_id]:
            raise KeyError(f"Gene ID '{gene_id}' not found under species '{species_id}'.")

        seq = fasta_dicts[species_id][gene_id]

        header = f">{focal_node}_{orthogroup}_{gene_num}_{gene_num2}_{gene_id}\n"
        fasta_content += header + seq + "\n"

    return fasta_content


def WriteAncestralFasta(focal_node, ancestral_genome, ortho_folder_path):
    outdir = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/AncestralGenomes")
    os.makedirs(outdir, exist_ok=True)
    output_fasta_path = os.path.join(outdir, f"{focal_node}.fasta")
    with open(output_fasta_path, "w") as fasta_file:
        fasta_file.write(ancestral_genome)

# Wrapper for multiprocessing OG by OG processing
def ProcessOrthogroupCurrent(index, gains_current, node, species_tree, ortho_folder_path, species_names):
    OG = gains_current[index]["Orthogroup"]
    return OG, GetAncestralGenes(OG, node, species_tree, ortho_folder_path, species_names)

# Build ancestral genome for a single node
def AncestralGenome(node, species_tree, ortho_folder_path, gains, species_names, n_threads):
    focal_node = node.name
    target_node = species_tree.search_nodes(name=focal_node)[0]

    # ancestors (including focal node)
    ancestors = [focal_node]
    up = target_node.up
    while up:
        ancestors.append(up.name)
        up = up.up

    # gains whose Gain Node is on the path from root to this node
    gains_current = [g for g in gains if g["Gain Node"] in ancestors]

    if not gains_current:
        return  # no OGs gained on path → no ancestral genes to reconstruct

    # multiprocessing over orthogroups
    n_threads = max(1, min(n_threads, cpu_count()))
    with Pool(processes=n_threads) as pool:
        results = pool.starmap(
            ProcessOrthogroupCurrent,
            [
                (idx, gains_current, node, species_tree, ortho_folder_path, species_names)
                for idx in range(len(gains_current))
            ],
        )

    ancestral_sequence = []
    orthogroup_names = []

    for OG, result in results:
        if result is not None:
            ancestral_sequence.append(result)
            # result is [ [gene1, gene2, ...] ], so result[0] is list of genes
            orthogroup_names.extend([OG] * len(result[0]))

    if not ancestral_sequence:
        return

    # Flatten to 1D list of gene IDs
    ancestral_gene_ids = np.concatenate(ancestral_sequence, axis=None).ravel().tolist()

    # species = part before '_' ; gene_id = full ID (e.g. "1_198")
    species = [g.split("_")[0] for g in ancestral_gene_ids]
    genes = [g for g in ancestral_gene_ids]

    all_species = sorted(set(species))

    # Build dictionary of all species FASTAs
    all_fasta_dicts = MakeAllSpeciesFastaDict(all_species, ortho_folder_path)

    # Build ancestral genome FASTA content and write to disk
    ancestral_genome = BuildAncestralGenome(species, genes, orthogroup_names, all_fasta_dicts, focal_node)
    WriteAncestralFasta(focal_node, ancestral_genome, ortho_folder_path)

def main(ortho_folder_path, n_threads):

    # Make ancestral genomes folder
    os.makedirs(os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/AncestralGenomes"), exist_ok=True)

    # Load numeric species tree
    species_tree_path = os.path.join(
        ortho_folder_path, "WorkingDirectory", "GladeWD", "SpeciesTree_rooted_node_labels.txt"
    )
    species_tree = ete3.Tree(species_tree_path, quoted_node_names=True, format=1)
    species_tree.name = "N0"
    species_names = species_tree.get_leaf_names()  # numeric species codes

    # Load gains (numeric)
    gains_file_path = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/GainsLossDuplication", "Gains.tsv")
    gains = []
    with open(gains_file_path, mode="r") as f:
        reader = csv.DictReader(
            f, fieldnames=["Gain Node", "Parent Node", "Orthogroup"], delimiter="\t"
        )
        next(reader)
        for row in reader:
            gains.append(row)

    # Reconstruct ancestral genome for every internal node
    node_list = []
    for node in species_tree.traverse("postorder"):
        if node.is_leaf():
            continue
        AncestralGenome(node, species_tree, ortho_folder_path, gains, species_names, n_threads)
        node_list.append(node.name)

    # Summaries from ancestral FASTAs
    fasta_folder = os.path.join(ortho_folder_path, "WorkingDirectory/GladeWD/AncestralGenomes")
    anc_dic, gene_counts, total_lengths, orthogroup_counts = MakeAncestralSpeciesFastaDict(
        node_list, fasta_folder
    )

    # Node-level summary
    anc_data_list = [
        {
            "Node": node,
            "Gene_Count": gene_counts.get(node, 0),
            "Number_of_amino_acids": total_lengths.get(node, 0),
        }
        for node in node_list
    ]
    sorted_anc_data = sorted(anc_data_list, key=lambda x: x["Node"])

    out_txt = os.path.join(fasta_folder, "AncestralGenomes.txt")
    with open(out_txt, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=sorted_anc_data[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(sorted_anc_data)

    # Orthogroup x node counts
    og_counts_list = [{"Orthogroup": og, **counts} for og, counts in orthogroup_counts.items()]

    # Fill missing entries with 0
    all_nodes = set(node_list)
    for og in og_counts_list:
        for n in all_nodes:
            og.setdefault(n, 0)

    og_counts_list = sorted(og_counts_list, key=lambda x: x["Orthogroup"])
    all_columns = sorted({k for d in og_counts_list for k in d.keys()})

    out_csv = os.path.join(fasta_folder, "Ancestral_HOG_counts.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_columns)
        writer.writeheader()
        for row in og_counts_list:
            writer.writerow(row)

    # Plot species tree for convenience
    PlotTree(species_tree, fasta_folder, "species_tree")

