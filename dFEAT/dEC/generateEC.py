from cresset import flare
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-l', '--ligand')
parser.add_argument('-p', '--protein')


args = parser.parse_args()

input_lig = args.ligand
input_prot = args.protein


print(input_lig)
print(input_prot)



# #Setting the Flare Project
# p = flare.Project()
# if (flare.main_window()):
#     flare.main_window().project = p

# #loading in the protein and ligands
# for l in ligand_list:
#     p.ligands.extend(flare.read_file(l))
# if protein_file is not None:
#     p.proteins.extend(flare.read_file(protein_file))