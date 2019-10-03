from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles, rdFMCS, rdChemReactions, rdMolDescriptors

import itertools
import numpy as np
import os
import csv

#####################

morphs_targets_dict = {
# JM Lab datasets:
#   '../datasets/input/THROMBIN/morph.in': "THROMBIN",
#    '../datasets/input/HSP90_3/morph.in': "HSP90_3",
#    '../datasets/input/HSP90_2/morph.in': "HSP90_2",
#    '../datasets/input/HSP90_1/morph.in': "HSP90_1",
#    '../datasets/input/FXR_1/morph.in': "FXR_1",
#    '../datasets/input/FXR_2/morph.in': "FXR_2",
#    '../datasets/input/ACK1/morph.in': "ACK1",

# # FEP+ datasets:
#   '../datasets/input/BACE/morph.in': "BACE",
#    '../datasets/input/JNK1/morph.in' : "JNK1",
#    '../datasets/input/TYK2/morph.in' : "TYK2",
#    '../datasets/input/SCHR_THROMBIN/morph.in' : "SCHR_THROMBIN",
#   '../datasets/input/CDK2/morph.in' : "CDK2",
  '../datasets/input/PTP1B/morph.in' : "PTP1B",
  '../datasets/input/MCL1/morph.in' : "MCL1"
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


    # replace ligand names with paths to their respective mol2 files:
    perturbations_paths = []
    for morph_pair in perturbation_list:
        member1_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[0]) + "/ligand.sdf"
        member2_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[1]) + "/ligand.sdf"
        perturbations_paths.append([member1_path, member2_path])

    return perturbations_paths


def build_reactions(perturbations_all_paths, MorphFilePath):
    # loop over each perturbation in the list and load the mol files:
    perturbation_reactions = []

    for perturbation_pair_path in perturbations_all_paths:
        path_to_mols = MorphFilePath.replace("morph.in","")+"poses/"

    # regenerate the perturbation (A>B):

        ligA = perturbation_pair_path[0].replace(path_to_mols, "").replace("/ligand.sdf","")
        ligB = perturbation_pair_path[1].replace(path_to_mols, "").replace("/ligand.sdf","")
        
        perturbation = str(ligA) + ">" + str(ligB)

    # read in mol files:
        mol1 = rdmolfiles.SDMolSupplier(perturbation_pair_path[0])[0]

        mol2 = rdmolfiles.SDMolSupplier(perturbation_pair_path[1])[0]
    
    # construct SMILES string from the two members:
#        member1_fullsmiles = Chem.MolToSmiles(mol1, allHsExplicit=True)
#        member2_fullsmiles = Chem.MolToSmiles(mol2, allHsExplicit=True)

#        reaction = str(member1_fullsmiles) + ">>" + str(member2_fullsmiles)
    # combine all results (SMILES):
        result = [perturbation, mol1, mol2]

        perturbation_reactions.append(result)
     
    return perturbation_reactions



def build_deltaFP(reactions):
    print("Building FPs and writing to CSV..")
    FP_column = np.arange(0, 256).tolist()
    FP_column = ["pfp" + str(item) for item in FP_column]

    PerturbationFingerprints = [
    "Perturbation", 
     
    "Member_Similarity (Dice)", 
    ]
    PerturbationFingerprints = [PerturbationFingerprints + FP_column]
    for reaction_members in reactions:
        pert = str(reaction_members[0])
    # deconstruct reaction smiles back into members:    
#        head, sep, tail = reaction_members[1].partition(">>")

    # take mol object from each member, retain hydrogens and override valency discrepancies
        member1 = reaction_members[1]
        member2 = reaction_members[2]
        member1.UpdatePropertyCache(strict=False)
        member2.UpdatePropertyCache(strict=False)

     # create bitstring of 256 bits for each member. 
        FP1 = (rdMolDescriptors.GetHashedAtomPairFingerprint(member1, 256))
        FP2 = (rdMolDescriptors.GetHashedAtomPairFingerprint(member2, 256))
        similarity = DataStructs.DiceSimilarity(FP1, FP2)

     # subtract and return reaction FP (=deltaFP) as list
        deltaFP = np.array(list(FP2)) - np.array(list(FP1))
        
     # join all the data together into one list and append to output:
        #reaction_members = pert

        result = [reaction_members[0]] + ([str(similarity)]) + deltaFP.tolist()

        PerturbationFingerprints.append(result)
    
    return PerturbationFingerprints





def collect_pFPs(morphs_targets_dict):
    # run all functions looping over the specified targets:
    for path, target in morphs_targets_dict.items():
        print("Working on", target+"..")

        
        perturbations_all_paths = read_morph_file(path)

        reactions = build_reactions(perturbations_all_paths, path)
        results_with_FPs = build_deltaFP(reactions)

        if not os.path.exists("./dFP_output"):
            os.makedirs("./dFP_output")

        with open('./dFP_output/perts_APFPs_'+str(target)+'.csv', 'w') as csvfile:
            for row in results_with_FPs:
                writer = csv.writer(csvfile)
                writer.writerow(row)
        print("Success; wrote to file.")
        print("############################")


collect_pFPs(morphs_targets_dict)
