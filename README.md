# GLADE

Gains, Losses, AncestralGenomes, Duplications, Evolution!

GLADE is a Python tool for reconstructing the full evolutionary history of orthogroups — including gene gains, losses, duplications, and ancestral gene repertoires — using only an OrthoFinder v3 results directory as input.

It maps every event onto the species tree and produces rich output for comparative genomics.

<p align="center">
  <img src="https://github.com/lauriebelch/GLADE/blob/main/GLADE.png" width="50%">
</p>

GLADE is a work in progress - proceed with caution!!!

## Table of contents
- [Installation](#What-is-GLADE)
- [Installation](#Installation)
- [How-to-use](#Simple-usage)
- [Output files](#Output-files)
- [Example-data](#Example-data)
- [Citation](#Citation)

## What is GLADE?

GLADE reconstructs the history of orthogroups — defined as sets of genes descended from a single gene in the most recent common ancestor — across a species tree.

Given a complete OrthoFinder v3 run, GLADE:
-identifies where each orthogroup first appeared (gain)
-detects losses
-identifies gene duplication events
-reconstructs ancestral gene content at every internal node,
-quantifies orthogroup size changes along every branch,
-outputs complete evolutionary histories for all orthogroups.

<p align="center">
<img src="glade_workflow.png" alt="workflow" width="700"/>
</p>

## Installation

GLADE requires Python 3.9 or later

I will add installation instructions for conda and manual

## Simple usage

```python GLADE.py -f path/to/orthofinder/results```

## Output files

GLADE produces a structured directory containing:

1. Gains, Losses, and Duplications (GainsLossDuplication/)
-Gains.tsv — where each orthogroup first appeared
-Loss_speciation.tsv — orthogroup losses due to speciation
-Loss_postduplication.tsv — losses after duplication events
-Duplications.tsv — duplication events with support values
-Branch_statistics.tsv — event counts per species-tree branch
-*_bybranch.tsv — expanded lists of orthogroups per event type
-extant_OG_counts.tsv — gene counts in extant species
-OrthogroupBranchChange.tsv — size changes per branch

2. AncestralGenomes/

One FASTA file per internal node containing reconstructed ancestral sequences
-AncestralGenomes.txt — summary statistics
-Ancestral_HOG_counts.csv — orthogroup copy numbers for all nodes

## Example data

## Citation

coming soon to biorxiv
Belcher L.J. et al. (2025) GLADE

