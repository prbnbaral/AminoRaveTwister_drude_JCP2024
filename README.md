# AminoRaveTwister Drude JCP2024

## Publication

**Grand canonical Monte Carlo and deep learning assisted enhanced sampling to characterize the distribution of Mg²⁺ and influence of the Drude polarizable force field on the stability of folded states of the twister ribozyme**

*J. Chem. Phys.* **161**, 225102 (2024)  
DOI: [10.1063/5.0241246](https://doi.org/10.1063/5.0241246)

## Authors

- **Prabin Baral**
- **Mert Y. Sengul** 
- **Alexander D. MacKerell Jr.**

## Overview

This repository contains a comprehensive molecular dynamics simulation workflow that combines OpenMM, metadynamics, and deep learning for iterative protein folding simulations, specifically developed for studying the twister ribozyme with Mg²⁺ ions using both additive and Drude polarizable force fields.

## Directory Structure

```
simulation_folder/
├── input/
│   ├── ml/
│   ├── md/
│   │   ├── system.psf
│   │   ├── system1_initial.crd
│   │   ├── system2_initial.crd
│   │   ├── system3_initial.crd
│   │   ├── system4_initial.crd
│   │   ├── system5_initial.crd
│   │   ├── system6_initial.crd
│   │   ├── system7_initial.crd
│   │   ├── system8_initial.crd
│   │   ├── system9_initial.crd
│   │   ├── system8_kcl_nomg_c36.psf
│   │   ├── system8_kcl_nomg_c36_initial.crd
│   │   ├── system8_kcl_nomg_drude.psf
│   │   ├── system8_kcl_nomg_drude_initial.crd
│   │   ├── step3_charmm2omm.str
│   │   ├── step4_equilibration.inp
│   │   ├── step5_production.inp
│   │   ├── toppar_drude.str
│   │   ├── drude_toppar_2020/
│   │   ├── toppar.str
│   │   └── toppar/
│   └── openmm/
│       ├── omm_vfswitch.py
│       ├── omm_restraints.py
│       ├── omm_barostat.py
│       ├── omm_readinputs.py
│       └── omm_readparams.py
├── output/
│   ├── meta/
│   ├── ml/
│   ├── md/
│   └── rave/
├── restraints/
│   ├── prot_ext.txt
│   └── prot_pos.txt
├── toppar/
├── drude_toppar_2020/
├── folding_package/
└── main_code_drude.py
```

## Directory Descriptions

### Input Files

**`toppar/`** - Contains all additive force field related files including `.rtf`, `.str`, and `.prm` files.

**`drude_toppar_2020/`** - Contains all Drude force field related `.str` files.

**`input/`** - Contains subfolders with input files organized by workflow component:
- **`openmm/`** - Files for running unbiased molecular dynamics simulations
- **`md/`** - Force field files (`toppar.str` for additive, `toppar_drude.str` for Drude), equilibration/production settings, and initial simulation configurations
- **`ml/`** - Initially empty; populated by the code with machine learning input files during execution

### Output Files

**`output/`** - Contains subfolders for output files (should be created empty):
- **`meta/`** - Metadynamics PLUMED input files generated during each iteration
- **`ml/`** - Machine learning related output files
- **`md/`** - Unbiased MD simulation results
- **`rave/`** - RAVE code output files

### Additional Components

**`restraints/`** - Contains restraint definition files:
- `prot_pos.txt` - Point restraint definitions
- `prot_ext.txt` - External restraining forces for secondary structure

**`folding_package/`** - Python package containing the complete iterative workflow code (imported by main script)

**`main_code_drude.py`** - Main execution script that coordinates deep learning, OpenMM, and metadynamics components

## Required Setup and Configuration

Before running the workflow, modify the following settings in `main_code_drude.py`:

### Path Configuration
```python
Path2script = "/home/username/python_codes"
```
Change this to the path where you've placed the complete workflow code.

### OpenMM Plugin Path
```python
Platform.loadPluginsFromDirectory('/home/username/modules/local/openmm/lib/plugins')
```
Update this path to point to your OpenMM PLUMED plugin location.

### Additional Settings
Customize metadynamics and deep learning parameters in the corresponding sections of `main_code_drude.py`.

## Dependencies

The following Python packages were used:

```
keras==2.8.0
MDAnalysis==2.2.0
numpy==1.23.1
pandas==1.4.2
parmed==3.4.3
OpenMM==7.7.0
```

## Usage

To run the simulation workflow:

```bash
python main_code_drude.py -i 0 1 2 3 4
```

This command executes:
- Iteration #0: Equilibration
- Iterations 1-4: Production runs

## Workflow Overview

This simulation workflow iteratively combines:
1. **OpenMM** for molecular dynamics simulations
2. **Metadynamics** for enhanced sampling
3. **Deep Learning** for analysis and optimization
4. **RAVE** for additional analysis

The system is designed for protein folding studies using both additive and Drude polarizable force fields.
