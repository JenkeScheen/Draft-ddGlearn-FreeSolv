from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles, rdFMCS, rdChemReactions, rdMolDescriptors

import itertools
import numpy as np
import os
import csv

#####################

morphs_targets_dict = {
## JM Lab datasets:
#   '../datasets/input/THROMBIN/morph.in': "THROMBIN",
#    '../datasets/input/HSP90_3/morph.in': "HSP90_3",
#    '../datasets/input/HSP90_2/morph.in': "HSP90_2",
#    '../datasets/input/HSP90_1/morph.in': "HSP90_1",
#    '../datasets/input/FXR_1/morph.in': "FXR_1",
#    '../datasets/input/FXR_2/morph.in': "FXR_2",
#   '../datasets/input/ACK1/morph.in': "ACK1",

# ## FEP+ datasets:
#   '../datasets/input/BACE/morph.in': "BACE",
#    '../datasets/input/JNK1/morph.in' : "JNK1",
#    '../datasets/input/TYK2/morph.in' : "TYK2",
#    '../datasets/input/SCHR_THROMBIN/morph.in' : "SCHR_THROMBIN",
#   '../datasets/input/CDK2/morph.in' : "CDK2",

   '../datasets/input/CDK2/morph.in' : "MCL1",
#   '../datasets/input/CDK2/morph.in' : "PTP1B",
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

    # replace ligand names with paths to their respective mol2 files:
    perturbations_paths = []
    for morph_pair in perturbation_list:
        member1_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[0]) + "/ligand.sdf"
        member2_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[1]) + "/ligand.sdf"
        perturbations_paths.append([member1_path, member2_path])

    return perturbations_paths


def DeleteSubstructs_unique(mol, submol):
    # this function is called in the rare case that the MCS substructure fits 'twice' in 
    # the larger perturbation member, deleting the whole structure resulting in a .. >> ..
    # this function defines a list of possible matches and takes only the first:
    matches = mol.GetSubstructMatches(submol)
    res = []
    for match in matches:
        match = [m for m in match if mol.GetAtomWithIdx(m).GetAtomicNum() > 1]
        exp_hs_to_add = []
        indices_to_remove = set()
        bonds_to_remove = set()
        mol_copy = Chem.Mol(mol)
        for b in mol_copy.GetBonds():
            is_ba_in_match = (b.GetBeginAtomIdx() in match)
            is_ea_in_match = (b.GetEndAtomIdx() in match)
            if (is_ba_in_match or is_ea_in_match):
                bonds_to_remove.add((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            if ((b.GetBeginAtom().GetAtomicNum() == 1 and is_ea_in_match)
                or (b.GetEndAtom().GetAtomicNum() == 1 and is_ba_in_match)):
                indices_to_remove.add(b.GetBeginAtomIdx() if is_ea_in_match
                                      else b.GetEndAtomIdx())
                continue
            if ((is_ba_in_match and (not is_ea_in_match))
                or (is_ea_in_match and (not is_ba_in_match))):
                if (is_ba_in_match):
                    a = b.GetEndAtom() if is_ba_in_match else b.GetBeginAtom()
                try:
                    exp_h_add = a.GetIntProp('__exp_h_add')
                except KeyError:
                    exp_h_add = 0
                exp_h_add += 1
                a.SetIntProp('__exp_h_add', exp_h_add)
        indices_to_remove_sorted = sorted(indices_to_remove.union(match),
                                          reverse=True)
        rwmol = Chem.RWMol(mol_copy)
        [rwmol.RemoveBond(ba, ea) for (ba, ea) in bonds_to_remove]
        [rwmol.RemoveAtom(i) for i in indices_to_remove_sorted]
        for a in rwmol.GetAtoms():
            try:
                exp_h_add = a.GetIntProp('__exp_h_add')
            except KeyError:
                continue
            a.SetNumExplicitHs(a.GetNumExplicitHs() + exp_h_add)
        mol_copy = Chem.AddHs(rwmol,
                              addCoords=(mol.GetNumConformers() > 0),
                              explicitOnly=True)
        Chem.SanitizeMol(mol_copy)
        mol_copy.ClearComputedProps()
        mol_copy.UpdatePropertyCache()
        res.append(mol_copy)
    return res


def build_reactions(perturbations_all_paths, MorphFilePath):
    # loop over each perturbation in the list and load the mol files:
    perturbation_reactions = []

    # print("Converting .PDB ligands to .MOL2..")
    # os.system("for f in "+MorphFilePath.replace("morph.in","")+"poses/*/ligand.pdb; do \
    #     molstring=$(echo $f | sed 's/.pdb/.mol2/g'); \
    #     obabel -i pdb $f -O $molstring; \
    #     done")
    # print("finished converting")

    for perturbation_pair_path in perturbations_all_paths:
        path_to_mols = MorphFilePath.replace("morph.in","")+"poses/"

    # regenerate the perturbation (A>B):

        ligA = perturbation_pair_path[0].replace(path_to_mols, "").replace("/ligand.sdf","")
        ligB = perturbation_pair_path[1].replace(path_to_mols, "").replace("/ligand.sdf","")
        perturbation = str(ligA) + ">" + str(ligB)

    # read in mol files:
        mol1 = rdmolfiles.SDMolSupplier(perturbation_pair_path[0])
        mol2 = rdmolfiles.SDMolSupplier(perturbation_pair_path[1])
        perturbation_pair = [mol1[0], mol2[0]]

    # generate MCS (taking into account substitutions in ring structures)
        print("Generating MCS for perturbation " + str(perturbation) + "..")
        MCS_object = rdFMCS.FindMCS(perturbation_pair, completeRingsOnly=True)
        MCS_SMARTS = Chem.MolFromSmarts(MCS_object.smartsString)

        if MCS_SMARTS == None:
            print("Could not generate MCS pattern")
            return

    # use SMARTS pattern to isolate unique patterns in each pair member
    # if multiple unique patterns exist in one molecule they are written as:
    # pattern1.pattern2 ('.' signifies a non-bonded connection)
        member1 = perturbation_pair[0]
        member2 = perturbation_pair[1]
        member1_stripped = AllChem.DeleteSubstructs(member1, MCS_SMARTS)
        member2_stripped = AllChem.DeleteSubstructs(member2, MCS_SMARTS)
        member1_stripped_smiles = Chem.MolToSmiles(member1_stripped, allHsExplicit=True)
        member2_stripped_smiles = Chem.MolToSmiles(member2_stripped, allHsExplicit=True)

    # when regular method creates a .. >> .., call alternative function:
        if len(member1_stripped_smiles) == 0 and len(member2_stripped_smiles) == 0:
            member1_stripped = DeleteSubstructs_unique(member1, MCS_SMARTS)
            member2_stripped = DeleteSubstructs_unique(member2, MCS_SMARTS)
            member1_stripped_smiles = Chem.MolToSmiles(member1_stripped[0])
            member2_stripped_smiles = Chem.MolToSmiles(member2_stripped[0])

    # if either member turns out empty, place a hydrogen for clarity (doesn't influence bits):            
        if len(member1_stripped_smiles) == 0:
            member1_stripped_smiles = "[H]"
        if len(member2_stripped_smiles) == 0:
            member2_stripped_smiles = "[H]"
    
    # if either member contains only a CH4, make it C-CH4 so the AP FP activates a C-C bond:
        if member1_stripped_smiles == "[CH4]":
            member1_stripped_smiles = "[C][CH4]"
        if member2_stripped_smiles == "[CH4]":
            member2_stripped_smiles = "[C][CH4]"
    # construct SMILES string from the two members (stripped and full):
        reaction = str(member1_stripped_smiles) + ">>" + str(member2_stripped_smiles)
        member1 = str(member1_stripped_smiles)
        member2 = str(member2_stripped_smiles)

        member1_fullsmiles = Chem.MolToSmiles(perturbation_pair[0])
        member2_fullsmiles = Chem.MolToSmiles(perturbation_pair[1])

    # combine all results (SMILES):
        result = [perturbation, reaction, member1_fullsmiles, member2_fullsmiles]

        perturbation_reactions.append(result)
     
    return perturbation_reactions


def build_deltaFP(reactions, target):
    #print("Building FPs and writing to CSV..")
    FP_column = np.arange(0, 256).tolist()
    FP_column = ["pfp" + str(item) for item in FP_column]

    PerturbationFingerprints = [
    "Perturbation", 
    "Reaction_SMILES", 
    "fullmember1",
    "fullmember2",
    "Member_Similarity (Dice)", 
    ]
    PerturbationFingerprints = [PerturbationFingerprints + FP_column]

     # create list of bitstring sizes to loop over:
    fp_sizes = list(np.arange(2, 500, 1))

    collected_dim_data = []
    for size in fp_sizes:
        print("size:", size)
        bitstrings = []
        for reaction_members in reactions:

            pert = str(reaction_members[0])
        # deconstruct reaction smiles back into members:    
            head, sep, tail = reaction_members[1].partition(">>")

        # take mol object from each member, retain hydrogens and override valency discrepancies
            # member1 = Chem.MolFromSmiles(head, sanitize=False)
            # member2 = Chem.MolFromSmiles(tail, sanitize=False)
            # member1.UpdatePropertyCache(strict=False)
            # member2.UpdatePropertyCache(strict=False)
            # instead of the MCS subtracted, take the whole ligands to FP on:
            member1 = Chem.MolFromSmiles(reaction_members[2], sanitize=False)
            member2 = Chem.MolFromSmiles(reaction_members[3], sanitize=False)
            member1.UpdatePropertyCache(strict=False)
            member2.UpdatePropertyCache(strict=False)

            FP1 = (rdMolDescriptors.GetHashedAtomPairFingerprint(member1, int(size)))
            FP2 = (rdMolDescriptors.GetHashedAtomPairFingerprint(member2, int(size)))
            similarity = DataStructs.DiceSimilarity(FP1, FP2)

         # subtract and return reaction FP (=deltaFP) as list
            deltaFP = np.array(list(FP2)) - np.array(list(FP1))
            
            bitstrings.append(deltaFP)
            
        # DF of all perts' bitstrings for this specific fp-size
        df_bitstrings = pd.DataFrame(bitstrings)

        
        sparsity = len(df_bitstrings.columns[(df_bitstrings == 0).all()])

        #PCA reduction:
        PCA.__init__
        pca = PCA(n_components=0.95)

        # Fit to df, reduce:          
        dim_size = len(pd.DataFrame(pca.fit_transform(df_bitstrings)).columns.values)
        
        print("null columns:", sparsity)
        print("dims after reduction:", dim_size) #?

        tmp_dict = {}
        tmp_dict["Number of bits in atom-pair fingerprint"] = size
        tmp_dict["Null-columns"] = sparsity
        tmp_dict["Dimensions after PCA reduction"] = dim_size
        collected_dim_data.append(tmp_dict)


        print("@@@@@@@@@@@@@@")
    collected_dim_data_df = (pd.DataFrame(collected_dim_data))
    collected_dim_data_df.to_csv("hashing_loss_check/sparsity_df_"+target+".csv")
    #sns.set()
    plt.rcParams.update({'font.size': 18})
    ax = sns.lineplot(x="Number of bits in atom-pair fingerprint", y="Dimensions after PCA reduction", data=collected_dim_data_df, color="skyblue")
    plt.ylabel("")
    plt.xlabel("")
    ax2 = plt.twinx()
    ax2 = sns.lineplot(x="Number of bits in atom-pair fingerprint", y="Null-columns", data=collected_dim_data_df, color="crimson")
    plt.ylabel("")

    plt.title(target)
    plt.axvline(x=256, linestyle="--", color="black")
    plt.tight_layout()
    plt.savefig("hashing_loss_check/sparsity_plot_"+target+".png", dpi=300)
    plt.show()
       

        
#        print("##########################################################################")
    return PerturbationFingerprints


import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
import seaborn as sns 


def collect_pFPs(morphs_targets_dict):
    # run all functions looping over the specified targets:
    for path, target in morphs_targets_dict.items():
        print(target, "###########################################################################################")
        
        perturbations_all_paths = read_morph_file(path)

        reactions = build_reactions(perturbations_all_paths, path)
        results_with_FPs = build_deltaFP(reactions, target)

        # if not os.path.exists("./dFP_output"):
        #     os.makedirs("./dFP_output")

        # with open('./dFP_output/perts_APFPs_'+str(target)+'.csv', 'w') as csvfile:
        #     for row in results_with_FPs:
        #         writer = csv.writer(csvfile)
        #         writer.writerow(row)



collect_pFPs(morphs_targets_dict)
