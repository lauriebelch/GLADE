# -*- coding: utf-8 -*-

import argparse
import os
import re
import math
import multiprocessing as mp
import tempfile
import ete4
import shutil
from common_functions import CleanGeneName, CleanSpeciesName, CheckBifurcating

## look up a species code from a species name (as written, or cleaned e.g. dots -> _)
def FindSpecies(SpeciesDict, name):
    if name in SpeciesDict:
        return SpeciesDict[name]
    if CleanSpeciesName(name) in SpeciesDict:
        return SpeciesDict[CleanSpeciesName(name)]
    raise KeyError(f"Species '{name}' not found in SpeciesIDs.txt. Known species: {list(SpeciesDict.keys())}")

## look up a gene code from a gene name (as written, or cleaned e.g. ( ) -> _)
## returns None if not found
def FindGene(gene_dict, gene):
    if gene in gene_dict:
        return gene_dict[gene]
    return gene_dict.get(CleanGeneName(gene))

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
            # also store the cleaned name (e.g. Rozella.allomycis -> Rozella_allomycis)
            SpeciesDict[CleanSpeciesName(species_base)] = species_code
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
                # also store the cleaned name, as OrthoFinder writes it in Orthogroups.tsv and the trees
                SequenceIDsDict[sp_code][CleanGeneName(original_gene)] = coded_gene

    # Convert Orthogroups.tsv to numeric-coded version
    with open(OG_Path) as OG_file, open(Output, "w") as outfile:

        header = next(OG_file).rstrip("\n")
        colnames = header.split("\t")[1:]
        # Convert species names -> numeric codes
        numeric_cols = [FindSpecies(SpeciesDict, s) for s in colnames]
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
                    new_gene = FindGene(SequenceIDsDict[species_code], gene)
                    if new_gene is None:
                        raise KeyError(
                            f"Gene '{gene}' not found in SequenceIDs for species code '{species_code}'.\n"
                            f"Column species: {colnames[pos]}"
                        )
                    replaced_genes.append(new_gene)

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
        species_codes = [FindSpecies(SpeciesDict, s) for s in colnames]
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
                    coded = FindGene(SequenceIDsDict[sp_code], gene)
                    if coded is None:
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
                    og_map[CleanGeneName(gene)] = coded
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
                code = FindGene(gene_to_code, norm)
                if code is None:
                    raise KeyError(
                        f"Gene '{norm}' not found in SequenceIDs. Original line: {line}"
                    )
                new_genes.append(code)

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
    new = FindGene(SequenceIDsDict[species_code], gene)
    if new is None:
        raise KeyError(
            f"Gene '{gene}' not found for species '{species}' (code={species_code}). "
            f"Full leaf: '{full_leaf}'"
        )
    return new

def _gene_tree_worker(args):
    chunk_file, out_file, SpeciesDict, SequenceIDsDict, used_X, OGLeafMap = args
    with open(chunk_file) as infile, open(out_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            og_name, tree = line.split(":", 1)
            # read the tree with ete4, so we only ever change leaf names
            # (this replaces the old regex, which missed characters like + and
            # skipped any leaf starting with "n")
            gene_tree = ete4.Tree(tree.strip(), parser=1)
            # Fetch per-OG map once (only used in -X mode)
            og_map = None
            if used_X:
                og_map = OGLeafMap.get(og_name)
                if og_map is None:
                    raise KeyError(
                        f"OG '{og_name}' not found in Orthogroups.tsv mapping. "
                        f"Cannot convert -X gene tree."
                    )
            for leaf in gene_tree.leaves():
                if used_X:
                    new = FindGene(og_map, leaf.name)
                    if new is None:
                        raise KeyError(
                            f"Leaf '{leaf.name}' not found in OG map for {og_name} "
                            f"(OrthoFinder was run with -X)."
                        )
                else:
                    new = convert_leaf(leaf.name, SpeciesDict, SequenceIDsDict)
                    if new == leaf.name:
                        raise KeyError(f"Leaf '{leaf.name}' in {og_name} does not start with a known species name.")
                leaf.name = new
            outfile.write(f"{og_name}: {gene_tree.write(parser=1, format_root_node=True)}\n")

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

    # stop if the species tree has a polytomy
    CheckBifurcating(tree)

    # Replace leaf names using SpeciesDict
    # (leaf names are matched as written, or cleaned, e.g. dots -> _)
    for leaf in tree.leaves():
        leaf.name = FindSpecies(SpeciesDict, leaf.name)   # numeric code ("0", "1", "2", ...)

    # Ensure folder exists
    os.makedirs(os.path.dirname(st_out), exist_ok=True)

    # Write the numeric version — internal node labels remain untouched
    tree.write(outfile=st_out, parser=1)

# check for -X flag
_X_FLAG_RE = re.compile(r'(^|\s)-X(\s|$)')

def orthofinder_used_X(ortho_folder_path):
    log_path = os.path.join(ortho_folder_path, "Log.txt")
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                if line.startswith("Command Line:"):
                    return bool(_X_FLAG_RE.search(line[len("Command Line:"):]))
    print("Note: could not find the OrthoFinder command line in Log.txt, assuming OrthoFinder was run without -X")
    return False


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
