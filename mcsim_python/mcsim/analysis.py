"""
Analysis Functions
"""

import math

import numpy as np
from mcsim.monte_carlo_numpy import calculate_distance


def get_all_pairwise_distances(configuration, box_length):
    """Get all the pairwise distances for a configuration"""

    distances = np.array([])
    num_atoms = len(configuration)
    for i in range(num_atoms):
        atom_distances = calculate_distance(configuration[i], configuration, box_length)
        atom_distances = np.delete(atom_distances, i)

        distances = np.append(distances, atom_distances)
    return distances


def rdf(values, n_bins, max_value, num_particles, box_length):
    """
    Compute the RDF for a set of values

    Parameters
    -----------
    values : np.ndarray
        The distance values to compute the RDF for.

    n_bins : int
        the number of bins to use for the histogram

    max_value : float
        The maximum value for which to compute the RDF.

    num_particles : int
        The number of particles in the system.

    box_length : float
        The box length

    Returns
    -------
    bin_centers : np.ndarray
        An array of the bin centers

    rdf : np.ndarray
        An array of the rdf values.

    """

    histogram, bins = np.histogram(values, bins=n_bins, range=(0, max_value))

    bin_size = bins[1] - bins[0]

    bin_centers = bins + bin_size / 2
    bin_centers = bin_centers[:-1]

    rdf = []

    rdf = histogram / (
        4 * math.pi * bin_centers**2 * bin_size * num_particles**2 / box_length**3
    )

    return bin_centers, rdf


def time_average_rdf(configurations, box_length, num_atoms, max_value):
    """Compute the time averaged rdf"""

    rdfs = []
    for configuration in configurations:
        distances = get_all_pairwise_distances(configuration, box_length)

        bins, configuration_rdf = rdf(distances, 200, max_value, num_atoms, box_length)
        rdfs.append(configuration_rdf)

    rdfs = np.array(rdfs)

    avg_rdf = rdfs.mean(axis=0)

    return bins, avg_rdf, rdfs
