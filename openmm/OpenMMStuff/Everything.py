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



def run_MM(
    new_CO_bond_length=0.151,
    new_CO_bond_strength=317000,
    new_CC_bond_length=0.151,
    new_CC_bond_strength=317000,
    replaced_CO_bond_length=0.151,
    replaced_CO_bond_strength=317000,
    replaced_CC_bond_length=0.151,
    replaced_CC_bond_strength=317000,
    bond_group_coefs=(1.01, 1.02, 1.03, 1.04),
    angle_group_coefs=(1.01, 1.02, 1.03, 1.04,1.05),
    use_forces = True):


    bond_group_coefs = {i:val for i, val in enumerate(bond_group_coefs)}
    angle_group_coefs = {i: val for i, val in enumerate(angle_group_coefs)}
    bond_group_coefs[-1] = 1
    angle_group_coefs[-1] = 1

    print(bond_group_coefs)
    print(angle_group_coefs)


    bond_group_file = open("./OpenMMStuff/BondGroup.txt")
    angle_group_file = open("./OpenMMStuff/AngleGroup.txt")

    bond_coefs = dict()
    for line in bond_group_file:
        parts = line.split(",")
        bond_coefs[int(parts[0])] = bond_group_coefs[int(parts[1])]


    angle_coefs = dict()
    for line in angle_group_file:
        parts = line.split(",")
        angle_coefs[int(parts[0])] = angle_group_coefs[int(parts[1])]

    print(bond_coefs)
    print(angle_coefs)




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
    bonds.addBond(1, 17, new_CO_bond_length, new_CO_bond_strength)
    bonds.addBond(3, 16, new_CC_bond_length, new_CC_bond_strength)
    bonds.setBondParameters(2, 1, 3, replaced_CC_bond_length, replaced_CC_bond_strength)
    bonds.setBondParameters(16, 16, 17, replaced_CO_bond_length, replaced_CO_bond_strength)

    for i in range(system.getForces()[0].getNumBonds()):
        bond = system.getForces()[0].getBondParameters(i)
        new_len = bond[2] * bond_coefs[i]
        system.getForces()[0].setBondParameters(i, bond[0], bond[1], new_len, bond[3])



    for i in range(system.getForces()[4].getNumAngles()):
        angle = system.getForces()[4].getAngleParameters(i)
        new_len = angle[3] * angle_coefs[i]
        system.getForces()[4].setAngleParameters(i, angle[0], angle[1], angle[2], new_len, angle[4])


    print("Hello!")
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
    start_file = open('start.txt', mode='w')
    end_file = open("end.txt", mode='w')
    for pos in state.getPositions(True):
        start_file.write("{}\n".format(pos))

    try:
        simulation.step(50000)
    except(OpenMMException):
        print("Value Error")
        return 10**6
        sys.exit()

    state = simulation.context.getState(positions=True, energy=True, forces=True, velocities=True, parameters=True, integratorParameters=True)
    end_file.write("{}\n".format(len(state.getPositions())))


    predicted_pos = []
    for pos in state.getPositions(True):
        str_pos = [str(p).split(" ")[0] for p in pos]
        end_file.write("{}\n".format(",".join(str_pos)))

        predicted_pos.append([float(n) for n in str_pos])

    end_file.write("{}\n".format(len(state.getForces())))

    predicted_forces = []
    for force in state.getForces(True):
        str_force = [str(f).split(" ")[0] for f in force]

        end_file.write("{}\n".format(",".join(str_force)))
        predicted_forces.append([float(n) for n in str_force])

    end_file.write("1")
    end_file.write("{}".format(str(state.getPotentialEnergy()).split(" ")[0]))
    predicted_energy = float(str(state.getPotentialEnergy()).split(" ")[0])

    print("Ended")

    target_file = open("./OpenMMStuff/target.txt", mode='r')

    num_molecules = int(target_file.readline())
    target_pos = []
    for _ in range(num_molecules):
        line = target_file.readline()
        target_pos.append([float(n) for n in line.split(",")])


    target_file.readline()
    target_forces = []
    for _ in range(num_molecules):
        line = target_file.readline()
        target_forces.append([float(n) for n in line.split(",")])




    predicted_pos = np.array(predicted_pos)
    target_pos = np.array(target_pos)

    # Step 1: Center positions (forces are not translated)
    target_centered = target_pos - target_pos.mean(axis=0)
    predicted_centered = predicted_pos - predicted_pos.mean(axis=0)

    # Step 2: Get optimal rotation matrix
    U = rmsd.kabsch(predicted_centered, target_centered)

    # Step 3: Rotate predicted positions and forces
    predicted_pos_aligned = np.dot(predicted_centered, U)
    predicted_forces_rot = np.dot(predicted_forces, U)


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
    force_rmse = np.sqrt(np.mean(np.sum((predicted_forces_rot - target_forces)**2, axis=1)))

    print(f"Position RMSD: {position_rmsd:.6f}")
    print(f"Force RMSE:   {force_rmse:.6f}")

    if(use_forces):
        return position_rmsd + force_rmse
    else:
        return position_rmsd
