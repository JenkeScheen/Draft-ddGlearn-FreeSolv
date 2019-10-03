import oddt
from oddt import fingerprints

import numpy as np
np.set_printoptions(edgeitems=10)

import itertools
import csv
import os
import subprocess

morphs_targets_dict = {
# JM Lab datasets:
#    '../datasets/input/THROMBIN/morph.in': "THROMBIN",
#    '../datasets/input/HSP90_3/morph.in': "HSP90_3",
#    '../datasets/input/HSP90_2/morph.in': "HSP90_2",
#    '../datasets/input/HSP90_1/morph.in': "HSP90_1",
#    '../datasets/input/FXR_1/morph.in': "FXR_1",
#    '../datasets/input/FXR_2/morph.in': "FXR_2",
#    '../datasets/input/ACK1/morph.in': "ACK1",

# FEP+ datasets:
   '../datasets/input/BACE/morph_small.in': "BACE",
#    '../datasets/input/JNK1/morph.in' : "JNK1",
#    '../datasets/input/TYK2/morph.in' : "TYK2",
#    '../datasets/input/SCHR_THROMBIN/morph.in' : "SCHR_THROMBIN",
#   '../datasets/input/CDK2/morph.in' : "CDK2"
    }

def read_morph_file(MorphFilePath):
    # read in morphfile:
    with open(MorphFilePath, 'rt') as morph_file:

    # read morph_pairs for cleaning
        block = list(itertools.takewhile(lambda x: "[protein]" not in x,
            itertools.dropwhile(lambda x: "morph_pairs" not in x, morph_file)))

        morph_list = [w.replace("\n", "").replace("\t","").replace(",", ", ") for w in block]
        morph_pairs = "".join(morph_list)
        
    # clean data and return as nested list:
        try:
            first_cleaned = (morph_pairs.replace("morph_pairs","").replace("=","").replace(",","\n"))
        except:
            print("Error in reading morph file, check if the line \"morph_pairs = ...\" is ordered vertically. Exiting..")
            return
        second_cleaned = (first_cleaned.replace(" ", "").replace(">",", "))
        molecule_pairs = second_cleaned.split("\n")
        perturbation_list = []
        for pert in molecule_pairs:
            if len(pert) != 0: 
                perturbation_list.append(pert.split(", ", 1))
        # print("Total amount of perturbations is: ",len(perturbation_list))
        # print("#####################################")
    
    # replace ligand names with paths to their respective mol files:
    perturbations_paths = []
    for morph_pair in perturbation_list:
        member1_path = MorphFilePath.replace("morph_small.in", "")+"poses/" + str(morph_pair[0]) + "/ligand.sdf"
        member2_path = MorphFilePath.replace("morph_small.in", "")+"poses/" + str(morph_pair[1]) + "/ligand.sdf"
        perturbations_paths.append([member1_path, member2_path])

    return perturbations_paths



def calcCE(perturbations_paths, MorphFilePath, PathToProtein):
    CEs = []

    for pert in perturbations_paths:

        member1 = pert[0]
        member2 = pert[1]
        protein = PathToProtein

        subprocess.call(["/home/jscheen/Flare/cresset/Flare/bin/pyflare", 
                        "generateEC.py", "--ligand", member1, "--protein", protein])

    
for path, target in morphs_targets_dict.items():

    pdbs = read_morph_file(path)
    

    PathToProtein = path.replace("morph_small.in", "protein/"+target+"/protein.pdb")
    CEs = calcCE(pdbs, path, PathToProtein)


    


