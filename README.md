# Numeric Simulation Laboratory - Davide Nicoli
My solutions to the assignments given in the 2021/2022 edition of the laboratory.

## Repository structure
The root directory of the project is organized in 5 directories plus some support files:

- `data/`: where configurations and dataset files are stored, for example the prime numbers lists for the random generator or the molecular dynamics settings;
- `results/`: where results from simulations and samplings are stored by default. Each exercise has its own subdirectory; some results are already cached, others will be computed on the fly;
- `solutions/`: it contains jupyter notebooks presenting the exercises' solutions;
- `src/`: it contains the C++ code 


## Project dependencies
The following software needs to be installed to build the project:
- a modern gcc compiler (supporting c++17);
- cmake >= 3.18;
- an anaconda distribution;
- an MPI distribution (both mpich and OpenMPI were tested).

Other build dependencies are specified in:
- "environment.yml" or "requirements.txt" (python dependencies; managed respectively by anaconda and pip);
- "conanfile.txt" (C++ dependencies; managed by conan).

Both these sets 

## How to build the project
After having installed the required build dependencies and cloned the current repository to a folder, the following steps must be followed:
### 1. Create the conda environment and install python dependencies
From the project root run:
```bash
conda env create -f environment.yml -n <Your environment name>
```
This command create a conda environment named \<Your environment name\> and installs the python dependencies (+ conan).
Now activate your environment:
```bash
conda activate <Your environment name>
```

### 2. Configure Conan and install C++ dependencies
C++ dependencies are managed by [conan](https://conan.io). Create a build directory and move inside it:
```bash
mkdir build
cd build
```
Conan has to be configured to properly install dependencies:






## How to display the exercises results
The jupyter notebooks that compute and diplay the solutions to the exercises are stored in `solutions/`. Every notebook should be already configured to run as intended without user configuration.

Many notebooks provide the option to use cached results setting the variable `USE_CACHED` to `True`: errors related to files not being found are usually related to the fact that computation has not already been performed. Try setting `USE_CACHED=False`.

