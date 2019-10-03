#!/bin/python

import os
import argparse
import sys
from rdkit import Chem
import py3Dmol
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from cresset import flare
import glob
from time import sleep

description_long="This will include a long description"

def print_progress_bar (percent, prefix = '', suffix = '', length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        percent   - Required  : percentage complete (Float)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = percent
    filledLength = int(percent)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %.2f%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if percent == 100.0: 
        print()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="falre_docking is a commandline interface to automatically dock ligands with Flare "
                                                 "free energies, PMFs and to view convergence.",
                                     epilog="flare_docking is built using pyflare ",
                                     prog="flare_docking")
    parser.add_argument('--description', action="store_true",
                        help="Print a complete description of this program.")

    parser.add_argument('--author', action="store_true",
                        help="Get information about the authors of this script.")

    parser.add_argument('--version', action="store_true",
                        help="Get version information about this script.")

    parser.add_argument('-l', '--ligands',
                        help="Directory which contains the ligands in either sdf or mol2 format.")

    parser.add_argument('-p', '--protein',
                         help="The file containing the pdb for the ligands to be docked to.")

    parser.add_argument('-o', '--flare_out', 
                        help="The filename of the Flare project to save the ligands and proteins docked.")

    parser.add_argument('--docking', action="store_true",
                         help="Actually run a docking calculation.")

    parser.add_argument('--max_poses', const=1, nargs= '?', type=int,
                        help="The number of poses to be evaluated with the docking. Default = 1")

    parser.add_argument('--grid',
                        help="This defines the docking grid, if not given then first loaded ligand in the ligands "
                        "directory will be used.")

    parser.add_argument('--docking_quality', const='ScoreOnly', nargs = '?', type=str,
                        help="The quality of the docking used. By default only the score of the current poses is evaluated.")

    #parsing all the arguemnts
    args, unknown = parser.parse_known_args()
    if len(unknown)>0:
        print('Some unkown arguemnts found')
        sys.exit(1)

    must_exit = False

    if args.description:
        print("%s\n" % description)
        must_exit = True

    if args.author:
        print("\nSome info about authorship")
        must_exit = True

    if args.version:
        print("Some info about versioning")
        must_exit = True

    if must_exit:
        sys.exit(0)

    #Defining the ligand list
    ligand_list = None
    if args.ligands:
        input_ligs = args.ligands

        #Check if the arguemnt provided is a directory
        if os.path.isdir(input_ligs):
            ligand_list = glob.glob(os.path.join(input_ligs,"*mol2"))
            if not ligand_list:
                ligand_list = glob.glob(os.path.join(input_ligs,"*smi"))
            if not ligand_list:
                print('No suitable .mol2 or .sdf ligand files found in the directory provided.')
                sys.exit(1)

        #If the argument povided is not a directory raise an error
        else:
            print("The argument -l %s is not a directory" %args.ligands)
            sys.exit(1)
            #glob.glob(os.path.join(input_ligs))

    #Protein information loaded
    protein_file = None
    if args.protein:
        print(args.protein)
        protein_file = args.protein
        if not os.path.isfile(protein_file):
            print('The file provided with the argument -p is not a file' %args.protein)
            sys.exit(1)

    flare_out = None
    if args.flare_out:
        flare_out = args.flare_out

    docking = False
    if args.docking:
        docking = True


    #Setting the Flare Project
    p = flare.Project()
    if (flare.main_window()):
        flare.main_window().project = p

    #loading in the protein and ligands
    for l in ligand_list:
        p.ligands.extend(flare.read_file(l))
    if protein_file is not None:
        p.proteins.extend(flare.read_file(protein_file))

    #now doing the docking
    if docking:

        #docking defaults:
        max_poses = None
        quality = 'ScoreOnly'
        grid = p.ligands[0]
        allowed_quality = ['ScoreOnly', 'Normal', 'VirtualScreen']

        if args.max_poses:
            max_poses = args.max_poses

        if args.docking_quality:
            q = args.docking_quality
            if q in allowed_quality:
                quality = q
            else:
                print('Docking quality provided with --docking_quality is not recognised \n')
                print('please chose one of the following list: \n')
                print(allowed_quality)

        if args.grid:
            print("This needs implementation still")

    #Initialise the docking and running the docking:
    print('Setting up the docking .....')
    dock = flare.Docking()
    dock.ligands = p.ligands
    dock.protein = p.proteins[0]
    dock.system.grid = grid
    #dock.max_poses = max_poses
    print('Running the docking ....')
    dock.start()
    while dock.is_running():

        status = dock.progress()
        if status[0][1] == 100.0:
            print_progress_bar(status[1][1], prefix = status[1][0], suffix = 'Complete', length = 100)
        else:
            print_progress_bar(status[0][1], prefix = status[0][0], suffix = 'Complete', length = 100)
        sleep(0.1)
        #print_progress_bar(status[1][1], prefix = status[1][0], suffix = 'Complete', length = 100)



    #saving the flare project
    p.save(flare_out)
