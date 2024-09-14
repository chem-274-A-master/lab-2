"""
Test running the simulation.
"""

import mcsim


def test_run_psl():
    """This just tests that simulation runs, not that any numbers are correct"""

    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]]

    mcsim.run_simulation(coords, 8, 3, 0.9, 10, use_numpy=False)


def test_run_numpy():
    """This just tests that simulation runs, not that any numbers are correct"""

    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]]

    mcsim.run_simulation(coords, 8, 3, 0.9, 10)
