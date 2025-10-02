from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import rmsd
import numpy as np
from openff.toolkit import Molecule
import networkx as nx
from networkx.algorithms import isomorphism as iso

def process_raw_graph(smiles):
    dir = os.path.join("Raw_Molecule", smiles)


    # getting a map from index to atom type from the raw atom
    atom_file = os.path.join(dir, "Atoms.txt")
    atom_file = open(atom_file, mode='r')

    raw_atoms = dict()
    for line in atom_file:
        parts = line.split(",")
        raw_atoms[int(parts[0])] = int(parts[1])


    bonds_file = os.path.join(dir, "Bonds.txt")
    bonds_file = open(bonds_file, mode='r')
    
    raw_bonds = []
    for line in bonds_file:
        parts = line.split(",")
        raw_bonds.append((int(parts[0]), int(parts[1])))
    
    return raw_bonds, raw_atoms

smiles = "CC(C)=Cc1ccccc1c2ccccc2C=O"

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

for i, force in enumerate(system.getForces()):
    print(i, force.__class__.__name__)


#system.getForces()[2].setParticleParameters(0,0,0,0)

open_mm_bonds = []
open_mm_atoms = {i:atom for i, atom in enumerate(openff_atoms)}
bonds = system.getForces()[0]
for i in range(system.getForces()[0].getNumBonds()):
    params = (bonds.getBondParameters(i))

    print("{} {} {}".format(params[0], params[1], i))
    open_mm_bonds.append((params[0], params[1]))

# Build graphs
G1 = nx.Graph()
G1.add_edges_from(open_mm_bonds)
nx.set_node_attributes(G1, open_mm_atoms, name="element")

raw_bonds, raw_atoms = process_raw_graph(smiles)
G2 = nx.Graph()
G2.add_edges_from(raw_bonds)
nx.set_node_attributes(G2, raw_atoms, name="element")

# Only allow mappings where `element` is identical
node_match = iso.categorical_node_match("element", None)

GM = iso.GraphMatcher(G1, G2, node_match=node_match)
print(GM.is_isomorphic())   # -> False, because Cl cannot map to H