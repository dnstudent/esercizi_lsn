# Numeric Simulation Laboratory - Davide Nicoli
My solutions to the assignments given in the 2021/2022 edition of the laboratory.

## Repository structure
The root directory of the project is organized in 5 directories plus some support files:

- `data/`: where configurations and dataset files are stored, for example the prime numbers lists for the random generator or the molecular dynamics settings;
- `results/`: where results from simulations and samplings are stored by default. Each exercise has its own subdirectory; some results are already cached, others will be computed on the fly;
- `solutions/`: it contains jupyter notebooks presenting the exercises' solutions;
- `src/`: it contains the C++ code 

## How to display the exercises results
Launch a jupyter server from the `solutions/` directory:
```bash
cd solutions
jupyter notebook
```

The jupyter notebooks that display the solutions to the exercises are stored under `solutions/`. Provided that a jupyter server is available to the user, every notebook should be ready for consultation.

If a full recomputation of results is wanted instead, [the following section](#how-to-prepare-the-project-for-recomputation) presents which steps need to be taken in order to prepare the source tree and install python runtime dependencies.
Running the notebooks then triggers recomputation, which may take a while.
Many notebooks provide the option to use cached results setting the variable `USE_CACHED` to `True`: errors related to files not being found are usually related to the fact that solutions are not in cache. Try setting `USE_CACHED=False` if that is the case.

## How to prepare the project for recomputation
### Project dependencies
The following software needs to be present on the computer in order to build the project:
- a modern compiler supporting c++17 (the project was tested with gcc-10, gcc-12 and clang-14);
- cmake >= 3.18;
- an anaconda distribution;
- an MPI distribution (both mpich and OpenMPI were tested).

Build and runtime dependencies, that will instead be installed in the following sections, are specified in:
- "environment.yml" (for python dependencies; managed by anaconda);
- "conanfile.txt" (for C++ dependencies; managed by conan, which will be installed by conan).


### How to build the project
After having installed the required build dependencies and cloned the current repository to a folder, the following steps should be followed:
#### 1. Create a conda environment and install python dependencies
From the project root run:
```bash
conda env create -f environment.yml -n <YourEnvName>
```
This command creates a conda environment named \<YourEnvName\> and installs the dependencies listed in "environment.yml".
Now activate your environment:
```bash
conda activate <YourEnvName>
```

#### 2. Configure Conan and install C++ dependencies
C++ dependencies are managed by [conan](https://conan.io).
It was installed by conda during the previous step, and now it has to be configured:
```bash
conan profile new default --detect
```
the prompt may show a warning message starting with:
```
************************* WARNING: GCC OLD ABI COMPATIBILITY ***********************
```
In case, I suggest following the directions given in the warning and execute:
```bash
conan profile update settings.compiler.libcxx=libstdc++11 default
```
to make the software work properly with code compiled for recent versions of C++.

Dependencies can now be installed.
Start with creating a build directory:
```bash
mkdir build
cd build
```
where the cmake configuration files and the build tree will reside.
Now install C++ dependencies running:
```
conan install -pr:h=default -pr:b=default --build=missing -s compiler.cppstd=17 ..
```

#### 3. Build the project using cmake
Inside the `build/` directory created above run:
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../
cmake --build .
cmake --install .
```
This should build and install the executables used to solve the exercises under the `bin/` directory in the project's root, ready to be used by the jupyter notebooks.

