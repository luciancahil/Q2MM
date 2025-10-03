from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import rmsd
import numpy as np
from openff.toolkit import Molecule
import networkx as nx
from networkx.algorithms import isomorphism as iso
import argparse



def write_positions(raw_path, proccessed_path, mapping, name):
    raw_pos = os.path.join(raw_path, name)
    raw_pos = open(raw_pos)

    processed_pos_file = os.path.join(proccessed_path, name)
    processed_pos_file = open(processed_pos_file, mode='w')

    atom_initial_pos_dict = dict()
    for line in raw_pos:
        parts = line.strip().split(",")
        atom_initial_pos_dict[int(parts[0])] = parts[1:]

    processed_initial_positions = {mapping[k]:atom_initial_pos_dict[k] for k in atom_initial_pos_dict.keys()}
    for i in range(len(processed_initial_positions)):
        processed_pos_file.write("{}\n".format(",".join(processed_initial_positions[i])))

    processed_pos_file.close()

def generate_BO_files(smiles, num_bonds):
    BO_dir = os.path.join("BO_data", "q2mm-{}".format(smiles))
    os.makedirs(BO_dir, exist_ok=True)

    # the dimension files
    dim_file = os.path.join(BO_dir, "dimension.csv")
    dim_file = open(dim_file, mode='w')

    dim_file.write("{}\n".format(num_bonds))
    
    for _ in range(num_bonds):
        dim_file.write("50000,1000000\n")

    # the next_value file.
    next_file = (os.path.join(BO_dir, "initial.csv"))
    next_file = open(next_file, mode='w')
    next_file.write(",".join(["60000" for _ in range(num_bonds)]))

    print("Settings for Bayes Optimization has been stored in {}".format(BO_dir))

def write_modified_bonds(raw_path, dir_path, mapping, bonds):
    mod_bond_file = os.path.join(raw_path, "modified_bonds.csv")
    mod_bond_file = open(mod_bond_file)

    bond_value_dict = dict()
    for line in mod_bond_file:
        parts = line.split(",")
        bond = [mapping[int(parts[0])], mapping[int(parts[1])]]
        bond.sort()
        bond_value_dict[(bond[0], bond[1])] = float(parts[2])
    
    new_bond_file = os.path.join(dir_path, "new_bonds.csv")
    new_bond_file = open(new_bond_file, mode='w')

    modified_bond_file = os.path.join(dir_path, "modified_bonds.csv")
    modified_bond_file = open(modified_bond_file, mode='w')

    for key in bond_value_dict.keys():
        if key in bonds:
            # index, atom1, atom2, bond distance
            modified_bond_file.write("{},{},{},{}\n".format(bonds.index(key), key[0], key[1], bond_value_dict[key]))
        else:
            new_bond_file.write("{},{},{}\n".format(key[0], key[1], bond_value_dict[key]))

    return len(bond_value_dict)

def write_to_processed(smiles, mapping, bonds):
    dir_path = os.path.join("Processed_Molecule", smiles)
    raw_path = os.path.join("Raw_Molecule", smiles)
    os.makedirs(dir_path, exist_ok=True)


    write_positions(raw_path, dir_path, mapping, "initial_pos.csv")
    write_positions(raw_path, dir_path, mapping, "target_pos.csv")


    num_special_bonds = write_modified_bonds(raw_path, dir_path, mapping, bonds)
    print("Information aboout the molecule has been stored in {}".format(dir_path))


    generate_BO_files(smiles, num_special_bonds)

def process_raw_graph(smiles):
    dir = os.path.join("Raw_Molecule", smiles)


    # getting a map from index to atom type from the raw atom
    atom_file = os.path.join(dir, "Atoms.csv")
    atom_file = open(atom_file, mode='r')

    raw_atoms = dict()
    for line in atom_file:
        parts = line.split(",")
        raw_atoms[int(parts[0])] = int(parts[1])


    bonds_file = os.path.join(dir, "Bonds.csv")
    bonds_file = open(bonds_file, mode='r')
    
    raw_bonds = []
    for line in bonds_file:
        parts = line.split(",")
        raw_bonds.append((int(parts[0]), int(parts[1])))
    
    return raw_bonds, raw_atoms


"""
python get_graph.py --smiles "CC(C)=Cc1ccccc1c2ccccc2C=O"
"""
parser = argparse.ArgumentParser(description="A simple script with arguments.")

parser.add_argument("--smiles", type=str, help="name of function", required=True)

args = parser.parse_args()
smiles = args.smiles


molecule = Molecule.from_smiles(smiles)

#molecule = Molecule.from_smiles("C")

mol=molecule.to_rdkit()

bonds = mol.GetBonds()



mol = Molecule.from_rdkit(mol)
openff_bonds = set([(bond.to_dict()['atom1'], bond.to_dict()['atom2']) for bond in mol._bonds])
openff_atoms = [atom.atomic_number for atom in mol._atoms]


# Create the SMIRNOFF template generator with the default installed force field (openff-2.1.0)
from openmmforcefields.generators import (
    SMIRNOFFTemplateGenerator,
)

smirnoff = SMIRNOFFTemplateGenerator(molecules=mol)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField

forcefield = ForceField(
    "amber/protein.ff14SB.xml",
    "amber/tip3p_standard.xml",
    "amber/tip3p_HFE_multivalent.xml",
)
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(smirnoff.generator)

# create a System with the non-bonded settings of mainline OpenFF force fields
# (9 Angstrom cut-off, switching distance applied at 8 Angstrom)
system = forcefield.createSystem(
    topology=mol.to_topology().to_openmm(),
    nonbondedCutoff=0.9 * openmm.unit.nanometer,
    switchDistance=0.8 * openmm.unit.nanometer,
)

#system.getForces()[2].setParticleParameters(0,0,0,0)

open_mm_bonds = []
open_mm_atoms = {i:atom for i, atom in enumerate(openff_atoms)}
bonds = system.getForces()[0]
for i in range(system.getForces()[0].getNumBonds()):
    params = (bonds.getBondParameters(i))

    open_mm_bonds.append((params[0], params[1]))

# Build graphs
raw_bonds, raw_atoms = process_raw_graph(smiles)
G1 = nx.Graph()
G1.add_edges_from(raw_bonds)
nx.set_node_attributes(G1, raw_atoms, name="element")


G2 = nx.Graph()
G2.add_edges_from(open_mm_bonds)
nx.set_node_attributes(G2, open_mm_atoms, name="element")

# Only allow mappings where `element` is identical
node_match = iso.categorical_node_match("element", None)

GM = iso.GraphMatcher(G1, G2, node_match=node_match)
assert(GM.is_isomorphic())
print("Raw Atom Index : OpenMM Atom Index")
print(GM.mapping)

write_to_processed(smiles, GM.mapping, open_mm_bonds)