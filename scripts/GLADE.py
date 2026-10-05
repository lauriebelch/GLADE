
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
from version import __version__

def main():    
    
    parser = argparse.ArgumentParser(
        prog='GLADE.py',
        description="-----------------------------------------------\n"
        "-----------------------------------------------\n"
        f"Welcome to GLADE v{__version__}\nGain, Loss, Ancestral gene sets, Duplication, Evolution!",
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
    parser.add_argument(
        '-s', '--seed',
        type=int,
        default=1,
        help="Random seed for choosing genes in the ancestral gene sets (default: 1). Same seed = same output"
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f"GLADE v{__version__}"
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    ortho_folder_path = args.folder
    n_threads = args.threads
    seed = args.seed

    # print welcome messages
    print("---------------------------------------------------")
    print(f"Welcome to GLADE v{__version__}\n")
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
    try:
        ConvertFiles.main(ortho_folder_path, n_threads)
    except (KeyError, ValueError) as e:
        # problems with the input files (e.g. names that don't match, polytomy in the species tree)
        print("\nGLADE stopped because of a problem with the OrthoFinder input:")
        print(e)
        sys.exit(1)
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
    print("Reconstructing ancestral gene sets...")
    AncestralGenome.main(ortho_folder_path, n_threads, seed)
    print("Ancestral gene sets reconstructed.")
    print("---------------------------------------------------")
    print("Calculating Branch statistics...")
    OrthoBranchChange.main(ortho_folder_path, n_threads)
    print("Branch statistics calculated.")
    print("Writing files...")
    ReStringFiles.main(ortho_folder_path, n_threads)
    print("Done.")

    # save the version and settings used, for reproducibility
    with open(os.path.join(ortho_folder_path, "GLADE_run_info.txt"), "w") as f:
        f.write(f"GLADE version\t{__version__}\n")
        f.write(f"Command\t{' '.join(sys.argv)}\n")
        f.write(f"Seed\t{seed}\n")
        f.write(f"Threads\t{n_threads}\n")
        f.write(f"Date\t{time.strftime('%Y-%m-%d %H:%M')}\n")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"---Finished! This run took {minutes} minutes {seconds} seconds.")
    print("---Files have landed in", os.path.join(ortho_folder_path, ""))
    print(f"---Thank you for choosing Orthofinder and GLADE v{__version__}.")
    print("---GLADE: Belcher L. & Kelly S. (2026), bioRxiv https://doi.org/10.64898/2026.01.27.702036")
    print("---OrthoFinder v3: Emms D.M., Liu Y., Belcher L., Holmes J. & Kelly S. (2025), bioRxiv https://doi.org/10.1101/2025.07.15.664860")

if __name__ == "__main__":
    main()
