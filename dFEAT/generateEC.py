from cresset import flare
import argparse

# Extracting ligand and file paths
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--ligand')
parser.add_argument('-p', '--protein')
args = parser.parse_args()

input_lig = args.ligand
input_prot = args.protein


# Making the project
p = flare.Project()
if (flare.main_window()):
     flare.main_window().project = p

# Loading
p.ligands.extend(flare.read_file(input_lig))
p.proteins.extend(flare.read_file(input_prot))

# Remove waters
#print("Waters")
waters = p.proteins[0].residues.find("WAT")
p.proteins[0].residues.remove(waters)

# # Docking
# print("Docking")
# max_poses = 1
# grid = p.ligands[0]
# dock = flare.Docking()

# dock.ligands = p.ligands
# dock.protein = p.proteins[0]
# dock.system.grid = grid
# dock.system.quality = flare.LeadFinderSystem.Quality.ScoreOnly
# dock.max_poses = max_poses
# dock.start()
# dock.wait()

# associate protein to ligand:
p.ligands[0].protein = p.proteins[0]

# Complementarity
print("Complementarity")
elec = flare.ElectrostaticComplementarity()
elec.ligands = p.ligands
elec.start()
elec.wait()

elec_list = []
for ligand in elec.ligands:
    elec_list.append(ligand.properties['Complementarity'].value)
    elec_list.append(ligand.properties['Complementarity r'].value)
    elec_list.append(ligand.properties['Complementarity rho'].value)
print(elec_list)
