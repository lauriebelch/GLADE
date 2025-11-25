#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import ast
import shutil

# 1. Load reverse dictionaries

def load_reverse_maps(wd):
    species_file = os.path.join(wd, "SpeciesIDs.txt")
    seq_file     = os.path.join(wd, "SequenceIDs.txt")

    code_to_species = {}
    with open(species_file) as f:
        for line in f:
            if ": " not in line:
                continue
            code, _, name = line.partition(": ")
            name = os.path.splitext(name.strip())[0]
            code_to_species[code.strip()] = name

    code_to_gene = {}
    with open(seq_file) as f:
        for line in f:
            if ":" not in line:
                continue
            coded, _, rest = line.partition(":")
            coded = coded.strip()                 # e.g. "1_350"
            original = rest.strip().split()[0]   # e.g. "gi|2848..."
            code_to_gene[coded] = original

    return code_to_species, code_to_gene

# 2. Utility functions

def parse_list_field(value):
    """Convert '["1","3"]'→['1','3'], '1'→['1']."""
    if value is None:
        return []
    value = str(value)
    if value == "":
        return []
    try:
        parsed = ast.literal_eval(value)
        if isinstance(parsed, list):
            return [str(x) for x in parsed]
    except Exception:
        pass
    return [value]


def convert_species_code(val, code_to_species):
    """Convert '1'→Mycoplasma_..., leave 'N1' or 'root' unchanged."""
    val = str(val)
    if val.isdigit():
        return code_to_species.get(val, val)
    return val


def convert_gene(val, code_to_species, code_to_gene):
    """Convert '1_203'→'Mycoplasma_xxx_gi|...|' if known."""
    val = str(val)
    if "_" not in val:
        return val

    sp, gid = val.split("_", 1)
    coded = f"{sp}_{gid}"

    if coded not in code_to_gene:
        return val

    species = code_to_species.get(sp, sp)
    original = code_to_gene[coded]
    return f"{species}_{original}"

# 3. Generic table converter

def convert_table(infile, outfile, species_cols, gene_cols, code_to_species, code_to_gene):
    if not os.path.exists(infile):
        return

    with open(infile) as f_in, open(outfile, "w") as f_out:
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            new_row = row.copy()

            # species-coded columns (numeric or already names)
            for col in species_cols:
                if col not in row:
                    continue
                items = parse_list_field(row[col])
                conv = [
                    convert_species_code(x, code_to_species)
                    for x in items
                    if x not in ("", "None")
                ]
                new_row[col] = ",".join(conv)

            # gene-coded columns (numeric gene IDs)
            for col in gene_cols:
                if col not in row:
                    continue
                items = parse_list_field(row[col])
                conv = [
                    convert_gene(x, code_to_species, code_to_gene)
                    for x in items
                    if x not in ("", "None")
                ]
                new_row[col] = ",".join(conv)

            writer.writerow(new_row)

# 4. Convert ancestral FASTAs

def convert_ancestral_fasta(indir, outdir, code_to_species, code_to_gene):
    if not os.path.isdir(indir):
        return

    os.makedirs(outdir, exist_ok=True)

    for fname in os.listdir(indir):
        if not fname.endswith(".fasta"):
            continue

        inf = os.path.join(indir, fname)
        outf = os.path.join(outdir, fname)

        with open(inf) as fin, open(outf, "w") as fout:
            for line in fin:
                if not line.startswith(">"):
                    fout.write(line)
                    continue

                header = line[1:].strip()
                parts = header.split("_")
                if len(parts) < 3:
                    fout.write(line)
                    continue

                sp = parts[-2]
                gid = parts[-1]
                coded = f"{sp}_{gid}"

                if coded in code_to_gene:
                    species = code_to_species.get(sp, sp)
                    original = code_to_gene[coded]
                    parts[-2] = species
                    parts[-1] = original

                fout.write(">" + "_".join(parts) + "\n")

def convert_branch_statistics(infile, outfile, code_to_species):
    if not os.path.exists(infile):
        return

    with open(infile) as fin, open(outfile, "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            branch = row["branch"]

            if "___" in branch:
                parent, child = branch.split("___", 1)

                # convert species-coded nodes ONLY if they are pure digits
                if parent.isdigit():
                    parent = code_to_species.get(parent, parent)
                if child.isdigit():
                    child = code_to_species.get(child, child)

                row["branch"] = f"{parent}___{child}"

            writer.writerow(row)


# 5. MAIN

def main(ortho_folder):
    wd    = os.path.join(ortho_folder, "WorkingDirectory")
    glade = os.path.join(wd, "GladeWD")

    out_gld = os.path.join(ortho_folder, "GainsLossDuplication")
    out_anc = os.path.join(ortho_folder, "AncestralGenomes")

    os.makedirs(out_gld, exist_ok=True)
    os.makedirs(out_anc, exist_ok=True)

    code_to_species, code_to_gene = load_reverse_maps(wd)

    # Gains
    convert_table(
        infile=os.path.join(glade, "GainsLossDuplication", "Gains.tsv"),
        outfile=os.path.join(out_gld, "Gains.tsv"),
        species_cols=["Gain Node", "Parent Node"],
        gene_cols=[],
        code_to_species=code_to_species,
        code_to_gene=code_to_gene,
    )

    # Loss
    convert_table(
        infile=os.path.join(glade, "GainsLossDuplication", "Loss_speciation.tsv"),
        outfile=os.path.join(out_gld, "Loss_speciation.tsv"),
        species_cols=["Species", "Child Node"],
        gene_cols=[],
        code_to_species=code_to_species,
        code_to_gene=code_to_gene,
    )

    # Duplications
    convert_table(
        infile=os.path.join(glade, "GainsLossDuplication", "Duplications.tsv"),
        outfile=os.path.join(out_gld, "Duplications.tsv"),
        species_cols=["speciestree_node"],
        gene_cols=["leaves1", "leaves2"],
        code_to_species=code_to_species,
        code_to_gene=code_to_gene,
    )

    # Post-duplication losses
    convert_table(
        infile=os.path.join(glade, "GainsLossDuplication", "Loss_postduplication.tsv"),
        outfile=os.path.join(out_gld, "Loss_postduplication.tsv"),
        species_cols=["Lost Species"],
        gene_cols=["Child Node"],
        code_to_species=code_to_species,
        code_to_gene=code_to_gene,
    )

    # Convert ancestral FASTAs
    convert_ancestral_fasta(
        indir=os.path.join(glade, "AncestralGenomes"),
        outdir=out_anc,
        code_to_species=code_to_species,
        code_to_gene=code_to_gene,
    )

    # Copy ancestral summary tables
    for fname in ["AncestralGenomes.txt", "Ancestral_HOG_counts.csv"]:
        src = os.path.join(glade, "AncestralGenomes", fname)
        dst = os.path.join(out_anc, fname)
        if os.path.exists(src):
            shutil.copy(src, dst)

    # Copy by-branch & stats files unchanged
    for fname in [
        "Gains_bybranch.tsv",
        "Loss_speciation_bybranch.tsv",
        "Duplications_bybranch.tsv",
        "Loss_postduplication_bybranch.tsv",
        "Branch_statistics.tsv",
    ]:
        src = os.path.join(glade, "GainsLossDuplication", fname)
        dst = os.path.join(out_gld, fname)
        if os.path.exists(src):
            shutil.copy(src, dst)

    bs_in  = os.path.join(glade, "GainsLossDuplication", "Branch_statistics.tsv")
    bs_out = os.path.join(out_gld, "Branch_statistics.tsv")
    convert_branch_statistics(bs_in, bs_out, code_to_species)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", help="Path to Orthofinder results folder")
    args = parser.parse_args()
    main(args.folder)
