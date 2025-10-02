from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import rmsd
import numpy as np
from openff.toolkit import Molecule


molecule = Molecule.from_smiles("CC(C)=Cc1ccccc1c2ccccc2C=O")

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


bonds = system.getForces()[0]
for i in range(system.getForces()[0].getNumBonds()):
    params = (bonds.getBondParameters(i))

    print("{} {} {}".format(params[0], params[1], i))


print("\n\n\n\n\n")
print("Copy the above to the text box at https://csacademy.com/app/graph_editor/")
print("This will show you how all the atoms connect, with each number in the node representing the index.")