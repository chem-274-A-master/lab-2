"""
How to run a simulation using mcsim
"""

import mcsim
import random
import time

from copy import deepcopy

# Set number of steps
num_steps = 10000

# Generate random coordinates in a box
coordinates, box_length = mcsim.generate_random_coordinates(num_atoms=500, density=0.9)

# Make a copy
starting = deepcopy(coordinates)

# Set the seed so we can compare our simulations

print("Running Python Standard Library ...")

start = time.time()
random.seed(0)
psl_coords = mcsim.run_simulation(
    coordinates, box_length, 3, 0.9, num_steps, freq=1000, use_numpy=False
)