# -*- coding: utf-8 -*-
"""
Integration test with a fresh OrthoFinder run.

Makes 6 proteomes from OrthoFinder's ExampleData with awkward names:
  - a species name with dots        (M.agalactiae.PG2.faa)
  - a species name starting with n  (nicotiana_haemocanis.fa)
  - gene IDs with ( )               (every 3rd gene in M. genitalium, e.g. gi|...|(chr1))
  - gene IDs with +                 (every 3rd gene in M. gallisepticum, e.g. gi|...|+)
runs OrthoFinder on them, then runs GLADE and checks:
  1. GLADE runs and every gene-tree leaf is converted
  2. duplication gene sets match OrthoFinder's Gene_Duplication_Events/Duplications.tsv
  3. a species-tree polytomy stops GLADE with a clear message
  4. two runs (1 and 8 threads) give identical ancestral gene sets

Needs OrthoFinder v3 on the PATH. Usage:
  python tests/fresh_orthofinder_run.py path/to/OrthoFinder/ExampleData work_folder
"""

import csv
import glob
import os
import shutil
import subprocess
import sys

import ete4

GLADE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "GLADE.py")
csv.field_size_limit(sys.maxsize)


def edit_fasta(src, dst, suffix):
    # add suffix to every 3rd gene ID
    count = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                count += 1
                if count % 3 == 0:
                    parts = line[1:].rstrip("\n").split(" ", 1)
                    line = ">" + parts[0] + suffix + (" " + parts[1] if len(parts) > 1 else "") + "\n"
            fout.write(line)


def run_glade(folder, threads="4"):
    return subprocess.run([sys.executable, GLADE, "-f", folder, "-t", threads], capture_output=True, text=True)


def main(of_example, work):
    prot = os.path.join(work, "proteomes")
    os.makedirs(prot)
    shutil.copy(os.path.join(of_example, "Mycoplasma_agalactiae.faa"), os.path.join(prot, "M.agalactiae.PG2.faa"))
    shutil.copy(os.path.join(of_example, "Mycoplasma_hyopneumoniae.faa"), prot)
    shutil.copy(os.path.join(of_example, "AdditionalSpecies", "M_arthritidis.fa"), prot)
    shutil.copy(os.path.join(of_example, "AdditionalSpecies", "M_haemocanis.fa"), os.path.join(prot, "nicotiana_haemocanis.fa"))
    edit_fasta(os.path.join(of_example, "Mycoplasma_genitalium.faa"), os.path.join(prot, "Mycoplasma_genitalium.faa"), "(chr1)")
    edit_fasta(os.path.join(of_example, "Mycoplasma_gallisepticum.faa"), os.path.join(prot, "Mycoplasma_gallisepticum.faa"), "+")

    subprocess.run(["orthofinder", "-f", prot, "-t", "8", "-a", "4", "-n", "B"], check=True, capture_output=True)
    results = os.path.join(prot, "OrthoFinder", "Results_B")

    # 1. GLADE runs, all leaves converted
    folder = os.path.join(work, "run1")
    shutil.copytree(results, folder)
    out = run_glade(folder, "1")
    print("1. GLADE on the new OrthoFinder run:", "OK" if out.returncode == 0 else "FAILED\n" + out.stdout + out.stderr)
    bad = 0
    with open(os.path.join(folder, "WorkingDirectory", "GladeWD", "Resolved_Gene_Trees.txt")) as fh:
        for line in fh:
            for leaf in ete4.Tree(line.split(":", 1)[1].strip(), parser=1).leaves():
                if not leaf.name.replace("_", "").isdigit():
                    bad += 1
    print("   unconverted gene-tree leaves:", bad)

    # 2. duplications vs OrthoFinder (compare gene names with OrthoFinder's cleaning applied)
    def clean(name):
        for char in [".", "(", ")", ":", ","]:
            name = name.replace(char, "_")
        return name
    of = {}
    for row in csv.DictReader(open(os.path.join(folder, "Gene_Duplication_Events", "Duplications.tsv")), delimiter="\t"):
        genes = {clean(g) for g in row["Genes 1"].split(", ") + row["Genes 2"].split(", ") if g}
        of[(row["Orthogroup"], row["Gene Tree Node"])] = genes
    glade = {}
    for row in csv.DictReader(open(os.path.join(folder, "GainsLossDuplication", "Duplications.tsv")), delimiter="\t"):
        genes = {clean(g) for g in row["leaves1"].split(",") + row["leaves2"].split(",") if g}
        glade[(row["Orthogroup"], row["genetree_node"])] = genes
    same = sum(1 for key in of if glade.get(key) == of[key])
    print(f"2. duplication gene sets identical to OrthoFinder: {same}/{len(of)}")

    # 3. species-tree polytomy
    folder3 = os.path.join(work, "polytomy")
    shutil.copytree(results, folder3)
    tree_file = os.path.join(folder3, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    t = ete4.Tree(open(tree_file).read(), parser=1)
    next(t.search_nodes(name="N3")).delete()
    open(tree_file, "w").write(t.write(parser=1, format_root_node=True))
    out = run_glade(folder3)
    print("3. species-tree polytomy:", "stops with message" if out.returncode != 0 and "bifurcating" in out.stdout
          else "NOT STOPPED")

    # 4. same seed, different threads -> identical ancestral gene sets
    folder4 = os.path.join(work, "run2")
    shutil.copytree(results, folder4)
    run_glade(folder4, "8")
    differ = 0
    for f in glob.glob(os.path.join(folder, "AncestralGenomes", "*.fasta")):
        if open(f).read() != open(os.path.join(folder4, "AncestralGenomes", os.path.basename(f))).read():
            differ += 1
    print(f"4. ancestral FASTAs that differ between two runs: {differ}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
