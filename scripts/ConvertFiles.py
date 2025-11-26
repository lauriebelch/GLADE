# -*- coding: utf-8 -*-

import argparse
import os
import re
import math
import multiprocessing as mp
import tempfile
import ete3


def File_Dictionaries(Input):
    """
    read orthofinder files for numeric conversion
    speciesdict codes species
    sequenceidsdict codes genes
    """
    # Identify input file paths
    OG_Path = os.path.join(Input, "Orthogroups", "Orthogroups.tsv")
    SeqIDs = os.path.join(Input, "WorkingDirectory", "SequenceIDs.txt")
    SpeciesIDs = os.path.join(Input, "WorkingDirectory", "SpeciesIDs.txt")
    Output = os.path.join(Input, "WorkingDirectory", "GladeWD", "Orthogroups.tsv")

    # Remove an old numeric file if present
    if os.path.exists(Output):
        os.remove(Output)

    # Read SpeciesIDs.txt
    SpeciesDict = {}         # maps species base name → species code
    Alt_SpeciesDict = {}     # maps species code → species base name

    with open(SpeciesIDs) as Species:
        for line in Species:
            if ": " not in line:
                continue
            key, _, value = line.partition(": ")
            raw_species = value.strip()  # remove newline/spaces
            # Remove file extension (any extension)
            species_base = os.path.splitext(raw_species)[0]
            species_code = key.strip()
            SpeciesDict[species_base] = species_code
            Alt_SpeciesDict[species_code] = species_base

    # Read SequenceIDs.txt
    SequenceIDsDict = {code: {} for code in Alt_SpeciesDict}

    with open(SeqIDs) as SeqID:
        for line in SeqID:
            if ":" not in line:
                continue
            key, _, value = line.partition(":")   # "1_198", " gi|..| description"
            sp_code = key.split("_")[0]           # species code, e.g. "1"
            coded_gene = key.strip()              # full coded gene like "1_198"

            # Extract FIRST token of the FASTA header
            first_token = value.strip().split()[0]
            original_gene = first_token

            if sp_code in SequenceIDsDict:
                SequenceIDsDict[sp_code][original_gene] = coded_gene

    # Convert Orthogroups.tsv to numeric-coded version
    with open(OG_Path) as OG_file, open(Output, "w") as outfile:

        header = next(OG_file).rstrip("\n")
        colnames = header.split("\t")[1:]
        # Convert species names -> numeric codes
        numeric_cols = [SpeciesDict[s] for s in colnames]
        # write new header
        outfile.write("Orthogroup\t" + "\t".join(numeric_cols) + "\n")
        # species_order for row processing
        species_order = numeric_cols[:]   # already numeric
        
        # Process each orthogroup row ----
        for line in OG_file:
            if not line.startswith("OG"):
                # Non-OG lines (summary rows etc.) copied as-is
                outfile.write(line)
                continue
            # Split row into OG name + species gene lists
            parts = line.rstrip("\n").split("\t")
            og_name = parts[0]
            species_fields = parts[1:]
            new_fields = []
            for pos, field in enumerate(species_fields):
                # Empty species column
                if field == "":
                    new_fields.append("")
                    continue

                species_code = species_order[pos]
                # Genes in this column separated by ", "
                genes = [g for g in field.split(", ") if g != ""]
                replaced_genes = []
                for gene in genes:
                    try:
                        new_gene = SequenceIDsDict[species_code][gene]
                        replaced_genes.append(new_gene)
                    except KeyError:
                        raise KeyError(
                            f"Gene '{gene}' not found in SequenceIDs for species code '{species_code}'.\n"
                            f"Column species: {colnames[pos]}"
                        )

                new_fields.append(", ".join(replaced_genes))
            outfile.write(og_name + "\t" + "\t".join(new_fields) + "\n")
    # Return dictionaries for use in next conversion steps
    return SpeciesDict, SequenceIDsDict

def Convert_Orthogroups_TXT(Input, SequenceIDsDict):
    """
    Convert Orthogroups.txt (1 OG per line, gene list format)
    """

    og_in  = os.path.join(Input, "Orthogroups", "Orthogroups.txt")
    og_out = os.path.join(Input, "WorkingDirectory","GladeWD", "Orthogroups.txt")

    if os.path.exists(og_out):
        os.remove(og_out)

    # Flatten SequenceIDsDict for easier lookup
    gene_to_code = {}
    for species_code, mapping in SequenceIDsDict.items():
        for original, coded in mapping.items():
            gene_to_code[original] = coded

    # Extract FASTA gene ID (first token)
    def normalize_gene(g):
        return g.strip().split()[0]

    # Process file
    with open(og_in) as infile, open(og_out, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            og_name, genes_str = line.split(":")
            genes = genes_str.strip().split()

            new_genes = []
            for g in genes:
                norm = normalize_gene(g)
                if norm not in gene_to_code:
                    raise KeyError(
                        f"Gene '{norm}' not found in SequenceIDs. Original line: {line}"
                    )
                new_genes.append(gene_to_code[norm])

            outfile.write(og_name + ": " + " ".join(new_genes) + "\n")



def convert_leaf(full_leaf, SpeciesDict, SequenceIDsDict):
    """
    Convert a leaf from a gene tree.
    Species names and gene IDs may contain underscores
    so identify species by checking which SpeciesDict key
    is the longest matching prefix of the leaf.
    """
    matches = []
    for species in SpeciesDict.keys():
        prefix = species + "_"
        if full_leaf.startswith(prefix):
            matches.append(species)
    # No species match → internal node (e.g., "n1")
    if not matches:
        return full_leaf
    # Choose longest match to avoid partial species names
    species = max(matches, key=len)
    # Extract gene ID (everything after "<species>_")
    gene = full_leaf[len(species) + 1:]
    species_code = SpeciesDict[species]
    if gene not in SequenceIDsDict[species_code]:
        raise KeyError(
            f"Gene '{gene}' not found for species '{species}' (code={species_code}). "
            f"Full leaf: '{full_leaf}'"
        )
    return SequenceIDsDict[species_code][gene]

def _gene_tree_worker(args):
    chunk_file, out_file, SpeciesDict, SequenceIDsDict = args

    with open(chunk_file) as infile, open(out_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            og_name, tree = line.split(":", 1)
            tree = tree.strip()

            leaves = set(re.findall(r"[A-Za-z0-9_\|\.\-]+", tree))

            for leaf in sorted(leaves, key=len, reverse=True):
                try:
                    new = convert_leaf(leaf, SpeciesDict, SequenceIDsDict)
                except KeyError:
                    if leaf.startswith("n"):
                        continue
                    raise

                tree = tree.replace(leaf, new)

            outfile.write(f"{og_name}: {tree}\n")


def Convert_Gene_Trees(Input, SpeciesDict, SequenceIDsDict, n_threads):
    """
    Convert Resolved Gene Trees using simple temp-file multiprocessing.
    """

    tree_in  = os.path.join(Input, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    tree_out = os.path.join(Input, "WorkingDirectory", "GladeWD", "Resolved_Gene_Trees.txt")

    if os.path.exists(tree_out):
        os.remove(tree_out)

    # Read all lines once
    with open(tree_in) as f:
        lines = [l for l in f if l.strip()]

    if not lines:
        open(tree_out, "w").close()
        return

    n_threads = max(1, min(n_threads, mp.cpu_count()))
    chunk_size = math.ceil(len(lines) / n_threads)

    # Temp directory for chunks
    tmp_dir = tempfile.mkdtemp(prefix="glade_trees_")

    chunk_files = []
    out_files   = []

    # Write chunk input files
    for i in range(n_threads):
        start = i * chunk_size
        end   = start + chunk_size
        chunk = lines[start:end]

        if not chunk:
            break

        chunk_path = os.path.join(tmp_dir, f"chunk_{i}.txt")
        out_path   = os.path.join(tmp_dir, f"chunk_{i}.out")

        with open(chunk_path, "w") as f:
            f.writelines(chunk)

        chunk_files.append(chunk_path)
        out_files.append(out_path)

    # Launch workers
    args = [
        (chunk_files[i], out_files[i], SpeciesDict, SequenceIDsDict)
        for i in range(len(chunk_files))
    ]

    with mp.Pool(processes=len(chunk_files)) as pool:
        pool.map(_gene_tree_worker, args)

    # Merge outputs in correct order
    with open(tree_out, "w") as final_out:
        for out_file in out_files:
            with open(out_file) as f:
                final_out.writelines(f)

    # Cleanup temp files
    for f in chunk_files + out_files:
        os.remove(f)
    os.rmdir(tmp_dir)

def Convert_Species_Tree(Input, SpeciesDict):
    """
    Convert Orthofinder species tree to numeric-coded version,
    while preserving ALL internal node labels exactly (N0, N1, ...).

    Only leaf names (species names) are replaced using SpeciesDict.
    """

    st_in  = os.path.join(Input, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    st_out = os.path.join(Input, "WorkingDirectory", "GladeWD", "SpeciesTree_rooted_node_labels.txt")

    # Load original tree with ETE — safest method
    tree = ete3.Tree(st_in, quoted_node_names=True, format=1)

    # Replace leaf names using SpeciesDict
    for leaf in tree.iter_leaves():

        # Orthofinder sometimes outputs leaf names with dots; normalize like SpeciesDict
        leaf_clean = os.path.splitext(leaf.name.replace(".", "_"))[0]

        if leaf_clean not in SpeciesDict:
            raise KeyError(
                f"Leaf name '{leaf.name}' (normalized '{leaf_clean}') "
                f"not found in SpeciesDict keys: {list(SpeciesDict.keys())}"
            )

        leaf.name = SpeciesDict[leaf_clean]   # numeric code ("0", "1", "2", ...)

    # Ensure folder exists
    os.makedirs(os.path.dirname(st_out), exist_ok=True)

    # Write the numeric version — internal node labels remain untouched
    tree.write(outfile=st_out, format=1)




def main(ortho_folder_path, n_threads):
    parent_output_file = os.path.join(ortho_folder_path, "WorkingDirectory", "GladeWD","GLADEfiles.tsv")
    os.makedirs(os.path.dirname(parent_output_file), exist_ok=True)
    SpeciesDict, SequenceIDsDict = File_Dictionaries(ortho_folder_path)
    Convert_Orthogroups_TXT(ortho_folder_path, SequenceIDsDict)
    Convert_Gene_Trees(ortho_folder_path, SpeciesDict, SequenceIDsDict, n_threads)
    Convert_Species_Tree(ortho_folder_path, SpeciesDict)

