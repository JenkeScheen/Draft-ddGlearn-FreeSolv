import pandas as pd 
from mordred import Calculator, descriptors
from rdkit import ML, Chem
from rdkit.Chem import rdmolfiles
import time

import itertools
import csv
import os
import subprocess

# open user-specified descriptors and create list:
# consult http://mordred-descriptor.github.io/documentation/master/descriptors.html
descriptors_raw = open("./descriptors/used_descriptors.txt", "r")
descriptors_raw_list = [line.split("\n") for line in descriptors_raw.readlines()]
descriptors_list = [desc[0] for desc in descriptors_raw_list]
print("Amount of descriptors: " + str(len(descriptors_list)))

morphs_targets_dict = {
# JM Lab datasets:
#     '../datasets/input/THROMBIN/morph.in': "THROMBIN",
#     '../datasets/input/HSP90_3/morph.in': "HSP90_3",
#     '../datasets/input/HSP90_2/morph.in': "HSP90_2",
#     '../datasets/input/HSP90_1/morph.in': "HSP90_1",
#     '../datasets/input/FXR_1/morph.in': "FXR_1",
#     '../datasets/input/FXR_2/morph.in': "FXR_2",
#     '../datasets/input/ACK1/morph.in': "ACK1",

# # FEP+ datasets:
#     '../datasets/input/BACE/morph.in': "BACE",
#     '../datasets/input/JNK1/morph.in' : "JNK1",
#     '../datasets/input/TYK2/morph.in' : "TYK2",
#     '../datasets/input/SCHR_THROMBIN/morph.in' : "SCHR_THROMBIN",
#     '../datasets/input/CDK2/morph.in' : "CDK2",
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
        print("#####################################")
    
    # replace ligand names with paths to their respective mol2 files:
    perturbations_paths = []
    for morph_pair in perturbation_list:
        member1_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[0]) + "/ligand.sdf"
        member2_path = MorphFilePath.replace("morph.in", "")+"poses/" + str(morph_pair[1]) + "/ligand.sdf"
        perturbations_paths.append([member1_path, member2_path])

    return perturbations_paths



def pdbPaths_to_MOLs(perturbations_paths, MorphFilePath):
    newpaths = []

    # convert pdb files to mol2 because then rdkit can read in bond orders:
#    print("Converting .PDB ligands to .MOL..")
#    os.system("for f in "+MorphFilePath.replace("morph.in","")+"poses/*/ligand.pdb; do \
#        molstring=$(echo $f | sed 's/.pdb/.mol/g'); \
#        obabel -i pdb $f -O $molstring; \
#        done")
#    print("finished converting")

    # load mol2 files, return nested list with molecule object pairs
    for pert in perturbations_paths:

        member1 = rdmolfiles.SDMolSupplier(pert[0])

        member2 = rdmolfiles.SDMolSupplier(pert[1])

        newpaths.append([member1[0], member2[0]])
    
    return newpaths


def ligFeaturesFromMols(mol_files, pert_paths, target, path, protein):
    
    # set up feature calculator, run per perturbation pair and calculate the feature difference 
    # i.e. subtract each member2 value from each member1 value:
    print("Generating deltaFeatures..")
    calc = Calculator(descriptors, ignore_3D=False)

    subtraction_values = []
    for pert in mol_files:
        # mordred descriptors:
        featured_members = calc.pandas(pert)
        featured_members.dropna()
        print(featured_members)
        featured_diff = featured_members_picked.diff(periods=1)

        # generate CEs using a pyflare subprocess:
        member1 = str(pert_paths[0][0])
      
        member2 = str(pert_paths[0][1])
        protein = str(protein.replace("morph.in", "protein/"+target+"/protein.pdb"))
        output1 = subprocess.check_output(["/home/jscheen/Flare/cresset/Flare/bin/pyflare", 
                    "generateEC.py", "--ligand", member1, "--protein", protein],)
        
        output1 = str(output1).replace("b'Complementarity\\n[", "").replace("]\\n'", "")
        output1 = output1.split(", ")
        output1 = [float(i) for i in output1]
        output2 = subprocess.check_output(["/home/jscheen/Flare/cresset/Flare/bin/pyflare", 
                    "generateEC.py", "--ligand", member1, "--protein", protein],)

        output2 = str(output2).replace("b'Complementarity\\n[", "").replace("]\\n'", "")
        output2 = output2.split(", ")
        output2 = [float(i) for i in output2]

        # calculate difference in CE metrics; make into series and add to features:
        diff_CE = [ round(output2[i]-output1[i], 3) for i in [0,1,2]]
        mordred = featured_diff.iloc[[1]]
        CE_dict = {"CE" : diff_CE[0], "CE_r" : diff_CE[1], "CE_rho" : diff_CE[2]}
        CE = pd.DataFrame(CE_dict, index=[1])
        result = pd.concat([mordred, CE], axis=1)
        subtraction_values.append(result)

    deltaFeatures = pd.concat(subtraction_values)

    # regenerate perturbation names:
    pert_names = []
    base_path = path.replace("morph.in", "poses/")
    for pert in pert_paths:
        member1 = pert[0].replace(base_path, "").replace("/ligand.sdf","")
        member2 = pert[1].replace(base_path, "").replace("/ligand.sdf","")
        pert_names.append(str(member1) + ">" + str(member2))
    
    # gather data, merge with perturbation names and output as CSV
    results_csv = [["Perturbation"] + deltaFeatures.columns.tolist()]
    deltaFeatures = deltaFeatures.values.tolist()
    
    # merge perturbation names with the corresponding deltaFeature data:
    for i in range(len(list(pert_names))):
        results_csv.append([pert_names[i]] + deltaFeatures[i])

    
    return results_csv

def write_to_file(data, target):
    # write to csv file 
    if not os.path.exists("./dFeatures_output"):
        os.makedirs("./dFeatures_output")

    with open('./dFeatures_output/deltaFeatures_'+target+'.csv', 'w') as csvfile:
        for row in data:
            writer = csv.writer(csvfile)
            writer.writerow(row)
    print("Success, wrote file to './dFeatures_output/deltaFeatures_"+target+".csv'")


def build_feats_on_dict(morphs_targets_dict):
    print(time.ctime())
    for path, target in morphs_targets_dict.items():
        print("#####################################")
        print("STARTING ON TARGET "+target)
        PathToProtein = path.replace("morph_small.in", "protein/"+target+"/protein.pdb")
        perturbation_paths = read_morph_file(path)

        mol2_files = pdbPaths_to_MOLs(perturbation_paths, path)

        data = ligFeaturesFromMols(mol2_files, perturbation_paths, target, path, PathToProtein)


        write_to_file(data, target)
    print(time.ctime())

build_feats_on_dict(morphs_targets_dict)

