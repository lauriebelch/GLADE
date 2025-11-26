
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:08:45 2024

Author: OrthoLaurie
"""
import argparse
import os
import time
import sys
import common_functions
import ConvertFiles
import GainAndLossAndDuplication
import BranchGainLossDuplication
import AncestralGenome
import OrthoBranchChange
import ReStringFiles

def main():    
    
    parser = argparse.ArgumentParser(
        prog='GLADE.py',
        description="-----------------------------------------------\n"
        "-----------------------------------------------\n"
        "Welcome to GLADE\nGain, Loss, AncestralGenome, Duplication, Evolution!",
        epilog="Example usage:\n  GLADE.py -f /path/to/Orthofinder/results -t 8\n"
        "-----------------------------------------------\n"
        "-----------------------------------------------",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-f', '--folder', 
        type=str, 
        required=True, 
        help="Path to the Orthofinder results folder"
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=8,
        help="Number of threads to use for multiprocessing (default: 8)"
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    ortho_folder_path = args.folder
    n_threads = args.threads

    # print welcome messages
    print("---------------------------------------------------")
    print("Welcome to Glade\n")
    if os.path.isdir(ortho_folder_path):
        print(f"Orthofinder results folder: {ortho_folder_path}")
    else:
        print(f"The folder '{ortho_folder_path}' does not exist.")
        sys.exit(1) 
        
    # record start time
    start_time = time.time()
    # Get the current local time
    current_hour = time.localtime().tm_hour
    print("\n")
    # greet the customer
    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")

    # Execute the imported scripts' main functions sequentially
    print("---------------------------------------------------")
    print("Converting Files...")
    ConvertFiles.main(ortho_folder_path, n_threads)
    print("Files converted.")
    print("---------------------------------------------------")
    print("Finding Gains, Losses, Duplications...")
    GainAndLossAndDuplication.main(ortho_folder_path, n_threads)
    print("Gains, Losses, Duplications found.")
    print("---------------------------------------------------")
    print("Mapping events to branches...")
    BranchGainLossDuplication.main(ortho_folder_path, n_threads)
    print("Events mapped to branches.")
    print("---------------------------------------------------")
    print("Reconstructing Ancestral Genomes...")
    AncestralGenome.main(ortho_folder_path, n_threads)
    print("Ancestral Genomes reconstructed.")
    print("---------------------------------------------------")
    print("Calculating Branch statistics...")
    OrthoBranchChange.main(ortho_folder_path, n_threads)
    print("Branch statistics calculated.")
    print("Writing files...")
    ReStringFiles.main(ortho_folder_path, n_threads)
    print("Done.")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"---Finished! This run took {minutes} minutes {seconds} seconds.")
    print("---Files have landed in", os.path.join(ortho_folder_path, ""))
    print("---Thank you for choosing Orthofinder and GLADE.") 
    print("---GLADE: Belcher L. & Kelly S. (2025), bioRxiv")
    print("---OrthoFinder v3: Emms D.M., Liu Y., Belcher L., Holmes J. & Kelly S. (2025), bioRxiv https://doi.org/10.1101/2025.07.15.664860")

if __name__ == "__main__":
    main()