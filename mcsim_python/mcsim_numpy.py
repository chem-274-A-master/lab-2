"""
How to run a simulation using mcsim
"""

import mcsim
import random
import numpy as np

from copy import deepcopy

# Set number of steps
num_steps = 10000

# Generate random coordinates in a box
coordinates, box_length = mcsim.generate_random_coordinates(num_atoms=500, density=0.9)

# Make a copy
starting = deepcopy(coordinates)

random.seed(0)
np_coords = mcsim.run_simulation(starting, box_length, 3, 0.9, num_steps, freq=1000)
