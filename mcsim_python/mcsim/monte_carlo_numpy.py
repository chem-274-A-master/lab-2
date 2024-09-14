"""
Functions for running Monte Carlo Simulation.
"""

import random
import math

import numpy as np

from copy import deepcopy


def calculate_total_virial(coordinates, box_length, cutoff):
    """
    Numpy version
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
    coordinates : np.ndarray
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

    e_total = 0.0
    i_position = coordinates[i_particle]

    distance_array = calculate_distance(coordinates, i_position, box_length)

    # Just so we don't use it for calculation
    distance_array[i_particle] = cutoff * 2

    less_than_cutoff = distance_array[distance_array < cutoff]

    interaction_energies = calculate_virial(less_than_cutoff)

    e_total = np.sum(interaction_energies)

    return e_total


def calculate_virial(r_ij):
    """
    Calculate the pairwise virial term.
    """

    r6_term = np.power(1 / r_ij, 6)
    r12_term = np.power(r6_term, 2)

    pairwise_virial = 24 * (2 * r12_term - r6_term)

    return pairwise_virial


def calculate_total_pressure(coordinates, box_length, cutoff, reduced_temperature):
    """Calculate the pressure based on configuration"""

    vir = calculate_total_virial(
        coordinates=coordinates, box_length=box_length, cutoff=cutoff
    )
    volume = box_length**3
    num_particles = len(coordinates)

    pressure = 1 / (3 * volume) * (3 * num_particles * reduced_temperature + vir)

    return pressure


def calculate_pressure_correction(num_atoms, box_length, cutoff):
    """
    Calculate the pressure correction.
    """

    volume = math.pow(box_length, 3)
    sig_by_cutoff3 = math.pow((1 / cutoff), 3)
    sig_by_cutoff9 = math.pow(sig_by_cutoff3, 3)

    p_tail = (2 / 3) * sig_by_cutoff9 - sig_by_cutoff3
    p_tail *= 16 * math.pi * math.pow(num_atoms, 2)
    p_tail /= 3 * math.pow(volume, 2)

    return p_tail


def calculate_LJ(r_ij):
    """
    Calculate the Lennard Jones energy based on a separation distance.

    Parameters
    ----------
    r_ij : float
        The distance between the two particles.

    Returns
    -------
    pairwise_energy : float
        The Lennard Jones potential energy.
    """

    r6_term = np.power(1 / r_ij, 6)
    r12_term = np.power(r6_term, 2)

    pairwise_energy = 4 * (r12_term - r6_term)

    return pairwise_energy


def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.
    Parameters
    ----------
    coord1, coord2: np.ndarray
        The atomic coordinates
    box_length : float
        The box length to be used for minimum image distance.

    Returns
    -------
    distance: float
        The distance between the two points.
    """
    coord_dist = coord1 - coord2

    if box_length:
        coord_dist = coord_dist - box_length * np.round(coord_dist / box_length)

    if coord_dist.ndim < 2:
        coord_dist = coord_dist.reshape(1, -1)

    coord_dist = coord_dist**2

    coord_dist_sum = coord_dist.sum(axis=1)

    distance = np.sqrt(coord_dist_sum)

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
    Numpy version - rewrite (calculate_total_energy)
    """

    total_energy = 0
    num_atoms = len(coordinates)

    for i in range(num_atoms):
        calc_coords = coordinates[i:]
        total_energy += calculate_pair_energy(calc_coords, 0, box_length, cutoff)

    return total_energy


def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system) - rewrite

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
        The pairwise interaction energy of he i_th particle with all other particles in the system.
    """

    e_total = 0.0
    i_position = coordinates[i_particle]

    distance_array = calculate_distance(coordinates, i_position, box_length)

    # Just so we don't use it for calculation
    distance_array[i_particle] = cutoff * 2

    less_than_cutoff = distance_array[distance_array < cutoff]

    interaction_energies = calculate_LJ(less_than_cutoff)

    e_total = np.sum(interaction_energies)

    return e_total
