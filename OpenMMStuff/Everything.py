from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import rmsd
import numpy as np

# Dear Future me:
# If you're reading this, make it more adaptable by passing in both an array of modified bonds added bonds.
"""for atom in pdb.topology.atoms():
    print(atom)"""

# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit import Molecule

def get_pos(dir, file_name):
    file = os.path.join(dir, file_name)
    file = open(file)

    positions = []
    for line in file:
        parts = line.split(",")
        positions.append([float(part) for part in parts])

    return positions


def get_initial_pos(dir):
    return get_pos(dir, "initial_pos.csv")

def get_target_pos(dir):
    return get_pos(dir, "target_pos.csv")

def get_new_bond_list(dir):
    file = os.path.join(dir, "new_bonds.csv")
    file = open(file)

    new_bonds = []
    for line in file:
        parts = line.split(",")
        new_bonds.append((int(parts[0]), int(parts[1]), float(parts[2])))

    return new_bonds


def get_modified_bonds(dir):
    file = os.path.join(dir, "modified_bonds.csv")
    file = open(file)

    modified_bonds = []
    for line in file:
        parts = line.split(",")
        modified_bonds.append((int(parts[0]), int(parts[1]), int(parts[2]), float(parts[3])))

    return modified_bonds

def run_MM(smiles, x):
    molecule = Molecule.from_smiles(smiles)
    dir = os.path.join("Processed_Molecule", smiles)

    # each element is start_index, end_index, bond_length
    new_bonds = get_new_bond_list(dir)

    # each element is bond_index, start_atom, end_atom, bond_length
    modified_bonds = get_modified_bonds(dir)

    assert( len(modified_bonds) + len(new_bonds) == len(x))
    positions = get_initial_pos(dir)
    target_pos = get_target_pos(dir)
    breakpoint()

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


    """
    bonds.addBond(1, 17, 0.24100492249018482, new_CO_bond_strength)
    bonds.addBond(3, 16, 0.16490954889347675, new_CC_bond_strength)
    bonds.setBondParameters(2, 1, 3, 0.14809143719020354, replaced_CC_bond_strength)
    bonds.setBondParameters(16, 16, 17, 0.1345190260703295, replaced_CO_bond_strength)
    """



    # set positions (dummy)

    pos_file = open("./OpenMMStuff/InitialPositionsOpenFF.txt", mode='r')

    positions = []

    for line in pos_file:
        parts = line.split(",")

        positions.append([float(part) for part in parts[1:]])

    integrator = LangevinMiddleIntegrator(0*kelvin, 1/picosecond, 0.0004*picoseconds)
    simulation = Simulation(mol.to_topology().to_openmm(), system, integrator)

    print(positions)
    simulation.context.setPositions(positions)




    simulation.minimizeEnergy()
    simulation.reporters.append(DCDReporter('output.dcd', 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
    state = simulation.context.getState(positions=True, energy=True, velocities=True, parameters=True, integratorParameters=True)

    try:
        simulation.step(50000)
    except OpenMMException as e:
        print("Error:", e)
        breakpoint()
        return 100

    state = simulation.context.getState(positions=True, energy=True, forces=True, velocities=True, parameters=True, integratorParameters=True)


    predicted_pos = []
    for pos in state.getPositions(True):
        str_pos = [str(p).split(" ")[0] for p in pos]

        predicted_pos.append([float(n) for n in str_pos])


    print("Ended")

    target_file = open("./OpenMMStuff/target.txt", mode='r')

    num_molecules = int(target_file.readline())
    target_pos = []
    for _ in range(num_molecules):
        line = target_file.readline()
        target_pos.append([float(n) for n in line.split(",")])






    predicted_pos = np.array(predicted_pos)
    target_pos = np.array(target_pos)

    # Step 1: Center positions (forces are not translated)
    target_centered = target_pos - target_pos.mean(axis=0)
    predicted_centered = predicted_pos - predicted_pos.mean(axis=0)

    # Step 2: Get optimal rotation matrix
    U = rmsd.kabsch(predicted_centered, target_centered)

    # Step 3: Rotate predicted positions and forces
    predicted_pos_aligned = np.dot(predicted_centered, U)


    # remove all Hydrogens
    predicted_pos_aligned = predicted_pos_aligned[0:18] * 10
    target_centered = target_centered[0:18] * 10
    dist_err_square = np.sum((predicted_pos_aligned - target_centered)**2, axis=1)



    errors = [(i, err.item()) for i,err in enumerate(dist_err_square)]
    errors.sort(key=lambda x:x[1])
    for e in errors:
        print(e)

    # Step 4: Compute errors
    position_rmsd = np.sqrt(np.mean(dist_err_square))

    print(f"Position RMSD: {position_rmsd:.6f}")


    return position_rmsd


# f([0.09091178829241708, 124037.47090186493, 0.4427850843188363, 712228.4922465099, 0.47804226345288636, 689976.6946725135, 0.14856227286757806, 434534.95339649945]) = 
# 
