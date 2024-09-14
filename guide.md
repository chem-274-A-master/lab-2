## Programming Questions

1. **What language is ASE written in?**
   ASE (Atomic Simulation Environment) is primarily written in Python (99.74% of the repository is Python according to repository analytics)

2. **Identify all of the classes we imported from ASE.**
   The following classes and functions were imported from ASE in the notebook:
   - `Atoms` (from `ase`)
   - `LennardJones` (from `ase.calculators.lj`)
   - `MaxwellBoltzmannDistribution` (from `ase.md.velocitydistribution`)
   - `VelocityVerlet` (from `ase.md.verlet`)
   - `Trajectory`, `write` (from `ase.io`)
   - `MDLogger` (from `ase.md`)

    kB and fs are also imported, but these are not classes that you instantiate. kB is a 

3. **ASE Lennard-Jones Calculator: Methods, Attributes, and Inheritance**

For commenting on methods and attributes, it can be an open discussion where students share what they found. Inheritance and encapsulation should be discussed for the class.
   - **Inheritance**: The Lennard-Jones calculator class inherits from `Calculator`, which is a generic base class for different types of calculators in ASE. Inheritance allows for code reuse and simplifies the extension of the base calculator class to other potential models. By inheriting from `Calculator`, Lennard-Jones gains access to common methods and attributes for calculators (it would be a good exercise to look at what is in calculator).
   - **Encapsulation/Modularity**: The Lennard-Jones calculator class encapsulates specific behavior for computing the Lennard-Jones potential, separating this functionality from other simulation logic. This makes it easy to switch between different potential models by simply changing the calculator without altering the rest of the code.
   
5. **Other ASE Class Example: VelocityVerlet**
This answers questions about VelocityVerlet: https://gitlab.com/ase/ase/-/blob/master/ase/md/verlet.py?ref_type=heads#L5
   - **Encapsulation/Modularity**: The `VelocityVerlet` class encapsulates the algorithm for time integration in molecular dynamics simulations. By separating this functionality into a class, ASE users can easily choose different integrators (like `VelocityVerlet` or other algorithms) by swapping out the integrator class.
   - **Inheritance**: `VelocityVerlet` inherits from `MolecularDynamics`, a base class that defines shared functionality for all integrators. This allows for consistent implementation of MD simulations while still enabling flexibility in the specific algorithm used. The base class handles common tasks like attaching loggers, while `VelocityVerlet` implements the specific integration method. **This concept will be used in Problem Set 2**

## Molecular Dynamics Exercises

1. **Reduced Density: Does this affect equilibration time?**
   - Yes, decreasing the reduced density to `0.009` tends to increase the equilibration time because the particles are more spread out, leading to fewer interactions and slower stabilization of energy. Starting with a random configuration can help with equilibration at lower densities, as it reflects the disordered state of a low-density system.

2. **Effect of Reduced Density on MD Properties (RDF, Diffusion Coefficient)**
   - At lower density, the RDF shows less pronounced peaks, more closely resembling the RDF of a gas. The diffusion coefficient increases.

3. **Simulating at Higher or Lower Temperature**
   - Increasing the temperature results in higher kinetic energy, leading to faster equilibration, broader RDF peaks, and a higher diffusion coefficient. Lower temperatures produce sharper RDF peaks due to reduced particle movement and stronger structural ordering, with a lower diffusion coefficient due to reduced thermal motion.

4. **Simulating a Different Substance with Lennard-Jones Parameters**
   - Students can discuss.
