"""
Functions for running Monte Carlo Simulation.
"""

import random
import math


def calculate_LJ(r_ij):
    """
    The LJ interaction energy between two particles.

    Computes the pairwise Lennard Jones interaction energy based
    on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.

    Returns
    -------
    pairwise_energy : float
        The pairwise Lennard Jones interaction energy in reduced units.

    """

    r6_term = math.pow(1 / r_ij, 6)
    r12_term = math.pow(r6_term, 2)

    pairwise_energy = 4 * (r12_term - r6_term)

    return pairwise_energy


def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.

    Parameters
    ----------
    coord1, coord2: list
        The atomic coordinates

    Returns
    -------
    distance: float
        The distance between the two points.
    """

    distance = 0
    for i in range(3):
        dim_dist = coord1[i] - coord2[i]

        if box_length:
            dim_dist = dim_dist - box_length * round(dim_dist / box_length)

        dim_dist = dim_dist**2
        distance += dim_dist

    distance = math.sqrt(distance)
    return distance


def calculate_tail_correction(num_particles, box_length, cutoff):
    """
    Calculate the long range tail correction
    """

    const1 = (8 * math.pi * num_particles**2) / (3 * box_length**3)
    const2 = (1 / 3) * (1 / cutoff) ** 9 - (1 / cutoff) ** 3

    return const1 * const2


def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Calculate the total Lennard Jones energy of a system of particles.

    Parameters
    ----------
    coordinates : list
        Nested list containing particle coordinates.

    Returns
    -------
    total_energy : float
        The total pairwise Lennard Jones energy of the system of particles.
    """

    total_energy = 0

    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):

            dist_ij = calculate_distance(
                coordinates[i], coordinates[j], box_length=box_length
            )

            if dist_ij < cutoff:
                interaction_energy = calculate_LJ(dist_ij)
                total_energy += interaction_energy

    return total_energy


def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.

    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.

    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates
    """

    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()

    atomic_coordinates = []

    for atom in coordinates:
        split_atoms = atom.split()

        float_coords = []

        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))

        atomic_coordinates.append(float_coords)

    return atomic_coordinates, box_length


def calculate_virial(r_ij):
    """
    Calculate the pairwise virial term.
    """

    r6_term = math.pow(1 / r_ij, 6)
    r12_term = math.pow(r6_term, 2)

    pairwise_virial = 24 * (2 * r12_term - r6_term)

    return pairwise_virial


# def calculate_total_virial(coordinates, box_length, cutoff):
#    """
#    Calculate the total Lennard Jones energy of a system of particles.
#
#    Parameters
#    ----------
#    coordinates : list
#        Nested list containing particle coordinates.
#
#    Returns
#    -------
#    total_energy : float
#        The total pairwise Lennard Jones energy of the system of particles.
#    """
#
#    total_virial = 0
#
#    num_atoms = len(coordinates)
#
#    for i in range(num_atoms):
#        for j in range(i + 1, num_atoms):
#
#            dist_ij = calculate_distance(
#                coordinates[i], coordinates[j], box_length=box_length
#            )
#
#            if dist_ij < cutoff:
#                interaction_virial = calculate_virial(dist_ij)
#                total_virial += interaction_virial
#
#    return total_energy


def calculate_total_virial(coordinates, box_length, cutoff):
    """
    Calculate total_virial
    """

    total_virial = 0
    num_atoms = len(coordinates)

    for i in range(num_atoms):
        calc_coords = coordinates[i:]
        total_virial += calculate_pair_virial(calc_coords, 0, box_length, cutoff)

    return total_virial


def calculate_pair_virial(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment

    Parameters
    ----------
    coordinates : list
        The coordinates for all particles in the system

    i_particle : int
        The particle number for which to calculate the energy

    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated.

    Returns
    -------
    e_total : float
        The pairwise virial of the i_th particle with all other particles in the system.
    """

    e_total = 0

    i_position = coordinates[i_particle]

    num_atoms = len(coordinates)

    for j_particle in range(num_atoms):
        if i_particle != j_particle:
            j_position = coordinates[j_particle]
            rij = calculate_distance(i_position, j_position, box_length)

            if rij < cutoff:
                e_pair = calculate_virial(rij)
                e_total += e_pair

    return e_total


def calculate_total_pressure(coordinates, box_length, cutoff, reduced_temperature):
    """Calculate the pressure based on configuration"""

    vir = calculate_total_virial(
        coordinates=coordinates, box_length=box_length, cutoff=cutoff
    )
    volume = box_length**3
    num_particles = len(coordinates)

    pressure = 1 / (3 * volume) * (3 * num_particles * reduced_temperature + vir)

    return pressure


def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment.

    Parameters
    ----------
    coordinates : list
        The coordinates for all the particles in the system.

    i_particle : int
        The particle index for which to calculate the energy.

    box_length : float
        The length of the simulation box.

    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated

    Returns
    -------
    e_total : float
        The pairwise interaction energy of the i-th particle with all other particles in the system.

    """

    e_total = 0

    i_position = coordinates[i_particle]

    num_atoms = len(coordinates)

    for j_particle in range(num_atoms):
        if i_particle != j_particle:
            j_position = coordinates[j_particle]
            rij = calculate_distance(i_position, j_position, box_length)

            if rij < cutoff:
                e_pair = calculate_LJ(rij)
                e_total += e_pair

    return e_total
