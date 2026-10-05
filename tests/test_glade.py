# -*- coding: utf-8 -*-
"""
Tests for GLADE, run on the ExampleData OrthoFinder results.

The edge-case tests copy ExampleData and edit it to recreate reported problems:
  - species names containing dots
  - gene IDs containing ( ) and +   (OrthoFinder writes ( ) as _ in trees/Orthogroups.tsv)
  - species names starting with "n"
  - a polytomy in the species tree
  - gene-tree polytomies (compared against OrthoFinder's own duplication output)
  - non-deterministic ancestral gene sets

Run from the repository root with:   pytest -v tests/
"""

import csv
import os
import shutil
import subprocess
import sys
import zipfile

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GLADE = os.path.join(REPO, "scripts", "GLADE.py")
RESULTS = "ExampleData/OrthoFinder/Results_ExampleDataGLADE"

csv.field_size_limit(sys.maxsize)


#### helper functions ########################################################

def make_example(tmp_path):
    # unzip a fresh copy of ExampleData and return the OrthoFinder results folder
    with zipfile.ZipFile(os.path.join(REPO, "ExampleData.zip")) as z:
        z.extractall(tmp_path)
    return os.path.join(tmp_path, RESULTS)


def run_glade(folder, *extra):
    cmd = [sys.executable, GLADE, "-f", folder, "-t", "2"] + list(extra)
    return subprocess.run(cmd, capture_output=True, text=True)


def replace_in_files(folder, old, new, files):
    # replace text in some OrthoFinder output files (like editing them with sed)
    for f in files:
        path = os.path.join(folder, f)
        with open(path) as fh:
            text = fh.read()
        with open(path, "w") as fh:
            fh.write(text.replace(old, new))


# the OrthoFinder files that GLADE reads
OF_FILES = ["WorkingDirectory/SpeciesIDs.txt",
            "Orthogroups/Orthogroups.tsv",
            "Species_Tree/SpeciesTree_rooted_node_labels.txt",
            "Resolved_Gene_Trees/Resolved_Gene_Trees.txt"]


def rename_species(folder, old, new):
    replace_in_files(folder, old, new, OF_FILES)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def count_rows(folder, name):
    return len(read_tsv(os.path.join(folder, "GainsLossDuplication", name)))


def all_leaves_converted(folder):
    # every leaf in GLADE's converted gene trees should look like "<species code>_<gene code>"
    import ete4
    path = os.path.join(folder, "WorkingDirectory", "GladeWD", "Resolved_Gene_Trees.txt")
    with open(path) as fh:
        for line in fh:
            og, newick = line.split(":", 1)
            for leaf in ete4.Tree(newick.strip(), parser=1).leaves():
                parts = leaf.name.split("_")
                if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
                    return False
    return True


#### tests ####################################################################

def test_example_data(tmp_path):
    # ExampleData: 338 duplications and 29 speciation losses
    folder = make_example(tmp_path)
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr
    assert count_rows(folder, "Duplications.tsv") == 338
    assert count_rows(folder, "Loss_speciation.tsv") == 29
    assert count_rows(folder, "Gains.tsv") == 332


def test_duplications_match_orthofinder(tmp_path):
    # GLADE duplication gene sets should match OrthoFinder's own
    # (Gene_Duplication_Events/Duplications.tsv). They reported 330/338 matching,
    # because GLADE only used the first two child clades at gene-tree polytomies.
    folder = make_example(tmp_path)
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr

    of = {}
    for row in read_tsv(os.path.join(folder, "Gene_Duplication_Events", "Duplications.tsv")):
        genes = set(row["Genes 1"].split(", ")) | set(row["Genes 2"].split(", "))
        genes.discard("")   # at polytomies OrthoFinder puts all genes in "Genes 1" and leaves "Genes 2" empty
        of[(row["Orthogroup"], row["Gene Tree Node"])] = genes

    glade = {}
    for row in read_tsv(os.path.join(folder, "GainsLossDuplication", "Duplications.tsv")):
        genes = set(row["leaves1"].split(",")) | set(row["leaves2"].split(","))
        glade[(row["Orthogroup"], row["genetree_node"])] = genes

    same = sum(1 for key in of if glade.get(key) == of[key])
    print(f"duplication gene sets identical to OrthoFinder: {same}/{len(of)}")
    assert same == len(of)


def test_species_name_with_dots(tmp_path):
    # species names with dots (e.g. Rozella.allomycis)
    folder = make_example(tmp_path)
    rename_species(folder, "Mycoplasma_agalactiae", "Mycoplasma.agalactiae")
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr
    assert all_leaves_converted(folder)
    assert count_rows(folder, "Duplications.tsv") == 338
    assert count_rows(folder, "Loss_speciation.tsv") == 29


def test_species_name_starting_with_n(tmp_path):
    # species names starting with a lower-case "n" (e.g. nicotiana_...) must still be converted
    folder = make_example(tmp_path)
    rename_species(folder, "Mycoplasma_agalactiae", "nicotiana_agalactiae")
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr
    assert all_leaves_converted(folder)
    assert count_rows(folder, "Loss_speciation.tsv") == 29


def test_gene_ids_with_special_characters(tmp_path):
    # OrthoFinder writes ( and ) as _ in Orthogroups.tsv and the gene trees,
    # but not in SequenceIDs.txt. Also '+' in gene IDs (strand-annotated IDs).
    folder = make_example(tmp_path)
    changes = [("gi|31541247|gb|AAP56549.1|", "gi|31541247|gb|AAP56549.1|(minus)", "gi|31541247|gb|AAP56549.1|_minus_"),
               ("gi|284811961|gb|ADB96864.1|", "gi|284811961|gb|ADB96864.1|+", "gi|284811961|gb|ADB96864.1|+")]
    for old, in_seqids, in_trees in changes:
        replace_in_files(folder, old, in_seqids, ["WorkingDirectory/SequenceIDs.txt"])
        replace_in_files(folder, old, in_trees, ["Orthogroups/Orthogroups.tsv", "Orthogroups/Orthogroups.txt",
                                                 "Resolved_Gene_Trees/Resolved_Gene_Trees.txt"])
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr
    assert all_leaves_converted(folder)
    assert count_rows(folder, "Duplications.tsv") == 338


def test_species_tree_polytomy_stops_with_error(tmp_path):
    # a polytomy in the species tree made GLADE silently lose most
    # speciation losses (29 -> 6). GLADE should now stop with a clear message.
    folder = make_example(tmp_path)
    tree_file = os.path.join(folder, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    with open(tree_file) as fh:
        tree = fh.read()
    # collapse node N1 into the root: ((A,B)N1,(C,D)N2)N0; -> (A,B,(C,D)N2)N0;
    import ete4
    t = ete4.Tree(tree, parser=1)
    n1 = next(t.search_nodes(name="N1"))
    n1.delete()
    with open(tree_file, "w") as fh:
        fh.write(t.write(parser=1, format_root_node=True))
    out = run_glade(folder)
    assert out.returncode != 0
    assert "bifurcating" in (out.stdout + out.stderr)


def test_same_seed_same_output(tmp_path):
    # ancestral gene sets were different on every run.
    # Same seed -> identical output, also with a different number of threads.
    folders = []
    for i, threads in enumerate(["2", "2", "4"]):
        folder = make_example(os.path.join(tmp_path, str(i)))
        out = subprocess.run([sys.executable, GLADE, "-f", folder, "-t", threads, "--seed", "7"],
                             capture_output=True, text=True)
        assert out.returncode == 0, out.stdout + out.stderr
        folders.append(folder)
    for name in ["N0.fasta", "N1.fasta", "N2.fasta", "AncestralGenomes.txt"]:
        files = [open(os.path.join(f, "AncestralGenomes", name)).read() for f in folders]
        assert files[0] == files[1] == files[2], name
    for name in ["OrthogroupBranchChange.tsv", "Duplications.tsv", "Loss_speciation.tsv"]:
        files = [open(os.path.join(f, "GainsLossDuplication", name)).read() for f in folders]
        assert files[0] == files[1] == files[2], name


def test_output_files_where_readme_says(tmp_path):
    # some outputs stayed in WorkingDirectory/GladeWD/ instead of the
    # locations given in the README
    folder = make_example(tmp_path)
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr
    for name in ["Gains.tsv", "Loss_speciation.tsv", "Loss_postduplication.tsv", "Duplications.tsv",
                 "Branch_statistics.tsv", "Gains_bybranch.tsv", "Loss_speciation_bybranch.tsv",
                 "Duplications_bybranch.tsv", "Loss_postduplication_bybranch.tsv",
                 "extant_OG_counts.tsv", "OrthogroupBranchChange.tsv"]:
        assert os.path.exists(os.path.join(folder, "GainsLossDuplication", name)), name
    for name in ["AncestralGenomes.txt", "Ancestral_HOG_counts.csv", "N0.fasta"]:
        assert os.path.exists(os.path.join(folder, "AncestralGenomes", name)), name
    # species names (not numeric codes) in the extant counts and the by-branch files
    with open(os.path.join(folder, "GainsLossDuplication", "extant_OG_counts.tsv")) as fh:
        assert "Mycoplasma_agalactiae" in fh.readline()
    for name in ["Gains_bybranch.tsv", "Duplications_bybranch.tsv", "OrthogroupBranchChange.tsv"]:
        for row in read_tsv(os.path.join(folder, "GainsLossDuplication", name)):
            assert not row["Branch"].split("___")[1].isdigit(), (name, row["Branch"])


def test_version_flag():
    out = subprocess.run([sys.executable, GLADE, "--version"], capture_output=True, text=True)
    assert out.returncode == 0
    assert "GLADE" in out.stdout and "1.0.0" in out.stdout


def test_gene_tree_polytomy_small_example():
    # a small hand-made example: species tree ((0,1)N1,2)N0 and a gene tree with a
    # 3-way polytomy at n1 containing two genes from species 0
    sys.path.insert(0, os.path.join(REPO, "scripts"))
    import ete4
    from common_functions import FindDuplications
    import GainAndLossAndDuplication
    species_tree = ete4.Tree("((0:1,1:1)N1:1,2:1)N0;", parser=1)
    gene_tree = ete4.Tree("((0_1:1,0_2:1,1_1:1)n1:1,2_1:1)n0;", parser=1)
    dupes = FindDuplications(gene_tree, species_tree, ["0", "1", "2"])
    assert len(dupes) == 1
    d = dupes[0]
    # all three genes under n1 are reported (the old code dropped the third child)
    assert sorted(d["leaves1"] + d["leaves2"]) == ["0_1", "0_2", "1_1"]
    # at a polytomy the order of duplications and losses is unknown, so no losses after duplication are inferred
    losses = GainAndLossAndDuplication.FindLossesAfterDuplications(gene_tree, species_tree, ["0", "1", "2"])
    assert losses == []


def test_postduplication_loss_branches():
    # species tree ((0,1)N1,(2,3)N2)N0. Each row of Loss_postduplication.tsv is one species that lost one copy;
    # species that form a clade are one loss, on the branch leading to that clade
    sys.path.insert(0, os.path.join(REPO, "scripts"))
    import ete4
    import BranchGainLossDuplication
    species_tree = ete4.Tree("((0:1,1:1)N1:1,(2:1,3:1)N2:1)N0;", parser=1)
    node_leaves, node_parent_list, _ = BranchGainLossDuplication.SpeciesTreeTraverse(species_tree)
    node_parent_dict = dict(node_parent_list)

    def branches(dupe, lost):
        # dupe: (Orthogroup, gene-tree node, genes in copy 1, genes in copy 2, species node); lost: {child node: [species]}
        og, focal, genes1, genes2, sn = dupe
        dupes = [{"Orthogroup": og, "genetree_node": focal, "leaves1": str(genes1), "leaves2": str(genes2)}]
        rows = [{"Orthogroup": og, "Focal Node": focal, "Lost Species": s, "Child Node": child, "Species Node": sn}
                for child, species in lost.items() for s in species]
        events = BranchGainLossDuplication.FindLossBranch(rows, node_leaves, node_parent_dict, dupes)
        return sorted(e["Branch_name"] for e in events)

    # duplication at N2, one copy lost in species 3 -> the branch to species 3 (not the branch to species 2)
    assert branches(("OG1", "n1", ["2_0", "3_0"], ["2_1"], "N2"), {"2_1": ["3"]}) == ["N2___3"]
    # duplication at N0, one copy lost in species 2 and 3 -> one loss on N0___N2 (not two losses on N0___N1)
    assert branches(("OG2", "n0", ["0_1", "1_0", "2_4", "3_3"], ["0_2", "1_1"], "N0"), {"n4": ["2", "3"]}) == ["N0___N2"]
    # as above but species 3 has neither copy: still one loss on N0___N2
    assert branches(("OG3", "n0", ["0_1", "1_0", "2_4"], ["0_2", "1_1"], "N0"), {"n4": ["2"]}) == ["N0___N2"]
    # each copy lost in a different species -> two losses, one per copy
    assert branches(("OG4", "n0", ["0_1", "2_1"], ["1_1", "2_2"], "N0"), {"n1": ["1"], "n2": ["0"]}) == ["N1___0", "N1___1"]


def test_outputs_agree_with_orthofinder_and_each_other(tmp_path):
    # cross-check the GLADE output tables against OrthoFinder's own files and each other
    import ete4
    folder = make_example(tmp_path)
    out = run_glade(folder)
    assert out.returncode == 0, out.stdout + out.stderr
    gld = os.path.join(folder, "GainsLossDuplication")

    species_tree = ete4.Tree(open(os.path.join(folder, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")).read(), parser=1)
    species_tree.name = "N0"
    species = list(species_tree.leaf_names())
    orthogroups = {row["Orthogroup"]: row for row in read_tsv(os.path.join(folder, "Orthogroups", "Orthogroups.tsv"))}

    # extant counts = OrthoFinder's Orthogroups.GeneCount.tsv
    gene_count = {row["Orthogroup"]: row for row in read_tsv(os.path.join(folder, "Orthogroups", "Orthogroups.GeneCount.tsv"))}
    extant = {row["Family_ID"]: row for row in read_tsv(os.path.join(gld, "extant_OG_counts.tsv"))}
    for og in gene_count:
        for sp in species:
            assert int(extant[og][sp]) == int(gene_count[og][sp])

    # gain node = most recent common ancestor of the species that have the orthogroup
    gains = read_tsv(os.path.join(gld, "Gains.tsv"))
    for row in gains:
        present = [sp for sp in species if orthogroups[row["Orthogroup"]][sp]]
        mrca = present[0] if len(present) == 1 else species_tree.common_ancestor(present).name
        assert row["Gain Node"] == mrca, row["Orthogroup"]

    # branch statistics add up to the event tables
    branch_stats = read_tsv(os.path.join(gld, "Branch_statistics.tsv"))
    dups = [row for row in read_tsv(os.path.join(gld, "Duplications.tsv")) if float(row["support"]) >= 0.5]
    assert sum(int(row["N_gains"]) for row in branch_stats) == len(gains)
    assert sum(int(row["N_speciation_losses"]) for row in branch_stats) == count_rows(folder, "Loss_speciation.tsv")
    assert sum(int(row["N_duplications"]) for row in branch_stats) == len(dups)

    # duplications per species-tree node = OrthoFinder's
    of_dups = {}
    for row in read_tsv(os.path.join(folder, "Gene_Duplication_Events", "Duplications.tsv")):
        if float(row["Support"]) >= 0.5:
            of_dups[row["Species Tree Node"]] = of_dups.get(row["Species Tree Node"], 0) + 1
    glade_dups = {}
    for row in dups:
        glade_dups[row["speciestree_node"]] = glade_dups.get(row["speciestree_node"], 0) + 1
    assert glade_dups == of_dups

    # orthogroup sizes on each branch agree with the extant and ancestral counts
    # (this caught a bug where the first orthogroup's ancestral sizes were read as 0)
    with open(os.path.join(folder, "AncestralGenomes", "Ancestral_HOG_counts.csv")) as fh:
        ancestral = {row["Orthogroup"]: row for row in csv.DictReader(fh)}

    def size(og, node):
        if node in species:
            return int(extant[og][node])
        return int(ancestral.get(og, {}).get(node, 0))

    for row in read_tsv(os.path.join(gld, "OrthogroupBranchChange.tsv")):
        assert int(row["family_size"]) == size(row["Orthogroup"], row["Focal_Node"]), row["Orthogroup_Branch"]
        assert int(row["parent_size"]) == size(row["Orthogroup"], row["Parent_Node"]), row["Orthogroup_Branch"]
        assert int(row["change"]) == int(row["family_size"]) - int(row["parent_size"])
