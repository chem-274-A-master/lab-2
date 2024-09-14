# Chem 274A - Lab 2

In this lab, you will be building on what we have learned about molecular simulations so far to run a molecular dynamics simulation.

**In order to complete this lab, you should use JupyterLab (not VSCode) in order to see the simulation visualizations.**

## Exercises

### Section 1 - Environment Creation
1. Clone the repository.
2. Use the `makefile` to create an environment for this lab.
    ```
    make environment
    ```
3. Activate the environment created by the `makefile`
    ```
    conda activate chem274A_lab2
    ```
4. Use the `makefile` to install `mcsim`
    ```
    make install
    ```

### Section 2 - Runnning Simulations
**In order to complete this lab, you should use JupyterLab (not VSCode) in order to see the simulation visualizations.**

Use the notebook `lennard_jones_md.ipynb` to run molecular dynamics simulations simulations. Follow the instructions in the notebook.

You can start Jupyter Lab by typing the following in your activated environment

```bash
jupyter lab --no-browser
```

This will start the Jupyter server. Click one of the links in the message to open Jupyter Lab in your browser/