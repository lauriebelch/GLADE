# -*- coding: utf-8 -*-

import argparse
import os
import re
import math
import multiprocessing as mp
import tempfile
import ete4
import shutil

_NUM_TOKEN_RE = re.compile(r'^[0-9]+(\.[0-9]+)?([eE][+\-]?[0-9]+)?$')

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

def Build_OG_Leaf_Map(Input, SpeciesDict, SequenceIDsDict):
    """
    For -X runs: build per-orthogroup mapping from gene ID to coded IDs, using Orthogroups/Orthogroups.tsv
    Returns:
        OGLeafMap: dict[og_name] -> dict[raw_gene] -> coded_gene
    """
    OG_Path = os.path.join(Input, "Orthogroups", "Orthogroups.tsv")
    OGLeafMap = {}
    with open(OG_Path) as og_file:
        header = next(og_file).rstrip("\n")
        colnames = header.split("\t")[1:]  # species names
        # Convert species names -> species codes (same as in File_Dictionaries)
        try:
            species_codes = [SpeciesDict[s] for s in colnames]
        except KeyError as e:
            raise KeyError(
                f"Species name {e!s} from Orthogroups.tsv header not found in SpeciesDict. "
                f"Header species: {colnames}"
            )
        for line in og_file:
            if not line.startswith("OG"):
                continue
            parts = line.rstrip("\n").split("\t")
            og_name = parts[0]
            species_fields = parts[1:]
            og_map = {}
            for pos, field in enumerate(species_fields):
                if field == "":
                    continue
                sp_code = species_codes[pos]
                genes = [g for g in field.split(", ") if g != ""]
                for gene in genes:
                    try:
                        coded = SequenceIDsDict[sp_code][gene]
                    except KeyError:
                        raise KeyError(
                            f"Gene '{gene}' not found in SequenceIDs for species code '{sp_code}'. "
                            f"OG: {og_name}, species column: {colnames[pos]}"
                        )
                    # Guard against weird duplicates inside the same OG
                    if gene in og_map and og_map[gene] != coded:
                        raise ValueError(
                            f"Ambiguous gene '{gene}' within {og_name}: "
                            f"{og_map[gene]} vs {coded}"
                        )
                    og_map[gene] = coded
            OGLeafMap[og_name] = og_map
    return OGLeafMap

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
    chunk_file, out_file, SpeciesDict, SequenceIDsDict, used_X, OGLeafMap = args
    with open(chunk_file) as infile, open(out_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            og_name, tree = line.split(":", 1)
            tree = tree.strip()
            leaves = set(re.findall(r"[A-Za-z0-9_\|\.\-]+", tree))
            # Fetch per-OG map once (only used in -X mode)
            og_map = None
            if used_X:
                og_map = OGLeafMap.get(og_name)
                if og_map is None:
                    raise KeyError(
                        f"OG '{og_name}' not found in Orthogroups.tsv mapping. "
                        f"Cannot convert -X gene tree."
                    )
            for leaf in sorted(leaves, key=len, reverse=True):
                # Skip internal node labels
                if leaf.startswith("n"):
                    continue
                # Skip numeric tokens (branch lengths/support values)
                if _NUM_TOKEN_RE.match(leaf):
                    continue
                if used_X:
                    if leaf in og_map:
                        new = og_map[leaf]
                    else:
                        # If it looks like a real gene label (letters or |), treat as an error.
                        # Otherwise ignore (acts like "leave unchanged").
                        if re.search(r"[A-Za-z\|]", leaf):
                            raise KeyError(
                                f"Leaf '{leaf}' not found in OG map for {og_name} "
                                f"(OrthoFinder was run with -X)."
                            )
                        continue
                else:
                    try:
                        new = convert_leaf(leaf, SpeciesDict, SequenceIDsDict)
                    except KeyError:
                        # Keep your existing behavior
                        raise
                tree = tree.replace(leaf, new)
            outfile.write(f"{og_name}: {tree}\n")

def Convert_Gene_Trees(Input, SpeciesDict, SequenceIDsDict, n_threads, used_X=False, OGLeafMap=None):
    """
    Convert Resolved Gene Trees using simple temp-file multiprocessing.
    If used_X=True, converts leaves using OGLeafMap[og_name][leaf].
    """
    tree_in  = os.path.join(Input, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
    tree_out = os.path.join(Input, "WorkingDirectory", "GladeWD", "Resolved_Gene_Trees.txt")
    if os.path.exists(tree_out):
        os.remove(tree_out)
    with open(tree_in) as f:
        lines = [l for l in f if l.strip()]
    if not lines:
        open(tree_out, "w").close()
        return
    if used_X and OGLeafMap is None:
        raise ValueError("used_X=True but OGLeafMap was not provided")
    n_threads = max(1, min(n_threads, mp.cpu_count()))
    chunk_size = math.ceil(len(lines) / n_threads)
    tmp_dir = tempfile.mkdtemp(prefix="glade_trees_")
    chunk_files = []
    out_files   = []
    try:
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
        args = [
            (chunk_files[i], out_files[i], SpeciesDict, SequenceIDsDict, used_X, OGLeafMap)
            for i in range(len(chunk_files))
        ]
        with mp.Pool(processes=len(chunk_files)) as pool:
            pool.map(_gene_tree_worker, args)
        with open(tree_out, "w") as final_out:
            for out_file in out_files:
                with open(out_file) as f:
                    final_out.writelines(f)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

def Convert_Species_Tree(Input, SpeciesDict):
    """
    Convert Orthofinder species tree to numeric-coded version,
    while preserving ALL internal node labels exactly (N0, N1, ...).

    Only leaf names (species names) are replaced using SpeciesDict.
    """

    st_in  = os.path.join(Input, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    st_out = os.path.join(Input, "WorkingDirectory", "GladeWD", "SpeciesTree_rooted_node_labels.txt")

    # Load original tree with ETE — safest method
    with open(st_in) as fh:
        tree = ete4.Tree(fh, parser=1)

    # Replace leaf names using SpeciesDict
    for leaf in tree.leaves():

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
    tree.write(outfile=st_out, parser=1)

# check for -X flag
_X_FLAG_RE = re.compile(r'(^|\s)-X(\s|$)')

def orthofinder_used_X(ortho_folder_path):
    log_path = os.path.join(ortho_folder_path, "Log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(f"Log.txt not found in {ortho_folder_path}")
    with open(log_path) as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise ValueError("Log.txt does not contain a command line")
    line = lines[1].strip()
    if not line.startswith("Command Line:"):
        raise ValueError("Unexpected Log.txt format (no 'Command Line:' on line 2)")
    cmd = line[len("Command Line:"):]
    return bool(_X_FLAG_RE.search(cmd))


def main(ortho_folder_path, n_threads):
    parent_output_file = os.path.join(ortho_folder_path, "WorkingDirectory", "GladeWD","GLADEfiles.tsv")
    os.makedirs(os.path.dirname(parent_output_file), exist_ok=True)
    used_X = orthofinder_used_X(ortho_folder_path)
    SpeciesDict, SequenceIDsDict = File_Dictionaries(ortho_folder_path)
    Convert_Orthogroups_TXT(ortho_folder_path, SequenceIDsDict)
    OGLeafMap = None
    if used_X:
        OGLeafMap = Build_OG_Leaf_Map(ortho_folder_path, SpeciesDict,SequenceIDsDict)
    Convert_Gene_Trees(
        ortho_folder_path,
        SpeciesDict,
        SequenceIDsDict,
        n_threads,
        used_X=used_X,
        OGLeafMap=OGLeafMap
    )
    Convert_Species_Tree(ortho_folder_path, SpeciesDict)
