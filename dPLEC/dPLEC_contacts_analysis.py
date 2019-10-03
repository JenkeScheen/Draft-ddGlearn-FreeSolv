import oddt
from oddt import fingerprints

import numpy as np
np.set_printoptions(edgeitems=10)

import itertools
import csv
import os
import pandas as pd

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
   '../datasets/input/BACE/morph.in': "BACE",
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
        print("Total amount of perturbations is: ",len(perturbation_list))
        print("#####################################")
    
    # replace ligand names with paths to their respective mol files:
    perturbations_paths = []
    for morph_pair in perturbation_list:
        member1_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[0]) + "/ligand.sdf"
        member2_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[1]) + "/ligand.sdf"
        perturbations_paths.append([member1_path, member2_path])

    return perturbations_paths



def pdbPaths_to_MOLs(perturbations_paths, MorphFilePath):
    newpaths = []

    # convert pdb files to mol because then rdkit can read in bond orders:
    # print("Converting .PDB ligands to .MOL..")
    # os.system("for f in "+MorphFilePath.replace("morph.in","")+"poses/*/ligand.pdb; do \
    #     molstring=$(echo $f | sed 's/.pdb/.mol/g'); \
    #     obabel -i pdb $f -O $molstring; \
    #     done")
    # print("finished converting")

    # load mol files, return nested list with molecule object pairs (NB: ODDT, not rdkit)
    for pert in perturbations_paths:

        member1 = next(oddt.toolkit.readfile("sdf", pert[0]))
        member2 = next(oddt.toolkit.readfile("sdf", pert[1]))

        

        newpaths.append([member1, member2])
    
    return newpaths

def calcFPs(perturbation_mols_list, path_to_protein, fpsize, pdbs):	
	# Generate protein ODDT mol object:
    try:
        print(path_to_protein)
        protein = next(oddt.toolkit.readfile("pdb", path_to_protein))
        
    except:
        print("Failed to generate protein object")
        return

    # Tell ODDT that the protein variable is a protein structure:
    protein.protein = True  

    # Make PLEC function:
    def PLEC_FP(ligand, protein):
        PLEC = fingerprints.PLEC(ligand, protein, 
        sparse=False,                           # equal bitsize for all processed ligands
        ignore_hoh=False,                       # incorporate H20 interactions in FP
        depth_ligand=1,                         # recommended ligand ECFP depth
        depth_protein=5,                        # recommended protein ECFP depth
        distance_cutoff=3.5,                    # maximum distance for atomic interaction
        size=fpsize)                            # 4096, 16384, 32768 or 65536 bits
        return np.array(PLEC, dtype="int64")    # change dtype from uint8 to int64 to allow subtraction


    # Call PLEC with set parameters for each ligand pair;
    # compute dPLEC by subtracting ligand A PLEC from ligand B PLEC:
    print("Computing dPLEC per perturbation..")
    dPLECs = []
    int_dfs = []
    for ligand_pair in perturbation_mols_list:

        for ligand in ligand_pair:
            

            ints = oddt.interactions.close_contacts(ligand.atom_dict, protein.atom_dict, 
                cutoff=3.5)
            pd.options.display.max_rows = 999
            
            idx_lig = [ atom_data[0] for atom_data in ints[0] ]
            atoms_lig = [ atom_data[5] for atom_data in ints[0] ]
            atoms_prot = [ atom_data[5] for atom_data in ints[1] ]
            residues_prot = [ atom_data[11] for atom_data in ints[1] ]
            df = pd.DataFrame(list(zip(idx_lig, atoms_lig, atoms_prot, residues_prot)), 
                columns=["Index_atom (ligand)","Ligand_atom", "Protein_atom", "Residue"])
            int_dfs.append(df)
            
        FP_Pair = [ PLEC_FP(ligand, protein) for ligand in ligand_pair ]
        deltaPLEC = FP_Pair[1] - FP_Pair[0]
        dPLECs.append(deltaPLEC.tolist())

    

    final_int_df = int_dfs[-1]
    final_int_df.to_csv("CloseContacts_control/sample_CloseContacts.csv", index=False)
    
    print(final_int_df)
    print(pdbs[-1][1])
    ###### Check which ligand PLEC is output.
    



def build_dPLECS(morphs_targets_dict):
    for path, target in morphs_targets_dict.items():
        print("#####################################")
        print("STARTING ON TARGET "+target)

        pdbs = read_morph_file(path)
        mol_paths = pdbPaths_to_MOLs(pdbs, path)

        PathToProtein = path.replace("morph.in", "protein/"+target+"/protein.pdb")
        
        deltaPLECs = calcFPs(mol_paths, PathToProtein, 16384, pdbs)

        
build_dPLECS(morphs_targets_dict)
