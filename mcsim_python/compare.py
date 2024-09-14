"""
How to run a simulation using mcsim
"""

import mcsim
import random
import time

from copy import deepcopy

import numpy as np

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
end = time.time()

elapsed_psl = end - start

print("\n\nRunning NumPy ...")

start = time.time()
random.seed(0)
np_coords = mcsim.run_simulation(starting, box_length, 3, 0.9, num_steps, freq=1000)
end = time.time()

elapsed_np = end - start
try:
    assert np.array_equal(psl_coords, np_coords)
    print("NumPy and PSL calculated the same values.")
except:
    print("NumPy and PSL did not calculate the same values.")

print(f"PSL ran {num_steps} trials in {elapsed_psl} seconds.")
print(f"NP ran {num_steps} trials in {elapsed_np} seconds.")
print(f"NP was {elapsed_psl/elapsed_np} faster than PSL.")
