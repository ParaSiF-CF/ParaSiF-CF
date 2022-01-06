# ParaSiF_CF - Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation Framework

The Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation Framework employs Code_Saturne to solve the computational fluid dynamics (CFD), FEniCS to solve the computational structure mechanics (CSM) and MUI for data exchange. It offers a platform where users can carry out fluid-structure interaction studies using supercomputers.

The framework uses a partitioned approach. It takes advantage of several features of the MUI toolbox:

• by showing good scalability for large simulations;

• by keeping independent the development of the coupled software. It allows for easier independent testing of each software and avoids potential incompatibilities between them, i.e. if they use a given library X but with a different version;

• by coupling multiple solvers which have different programming language interfaces, e.g C, C++, FORTRAN and Python;

• by relying on a "Plug and play" strategy. One solver can be replaced by another one incrementally and without the need of recompiling if the MUI interface is used as a common adaptor;

• by using multiple software that could have incompatible licenses. It would exploit the dual licensing of the MUI library, as the solvers are never mixed source-wise or binary-wise.

**This framework is a beta-software at the moment and under active development**.

## Licensing

Copyright (C) 2021 Computational Engineering and Environment Group, Scientific Computing Department, Science and Technology Facilities Council, UK Research and Innovation. All rights reserved.

This code is licensed under the GNU General Public License version 3

## Acknowledgements
The Parallel Partitioned Multi-physical Simulation Framework is developed by the [Scientific Computing Department](https://www.scd.stfc.ac.uk/) of the [Science and Technology Facilities Council](https://stfc.ukri.org/).

## Publication


## Contact

Should you have any question please do not hesitate to contact the developers directly.

## Installation

A way on how to install Code_Saturne (V6.0.6), FEniCS (V2019.1.0) and MUI (V1.1.3) and their dependencies is provided under https://github.com/ParaSiF-CF/ParaSiF-CF_Archer2_install

## Source and export before run ParaSiF_CF cases

Assuming that $INSTALL_FOLDER has been defined as during the aforementioned installation, several settings have to be added:

```bash
source ${INSTALL_FOLDER}/FEniCS/V2019.1.0/dolfin/build/share/dolfin/dolfin.conf
export PATH=${INSTALL_FOLDER}/SATURNE/V6.0.6/code_saturne-6.0.6/arch/Linux/bin:$PATH
export PYTHONPATH=/path/to/ParaSiF_CF/src/CSM/FEniCS/V2019.1.0:$PYTHONPATH
```

## Demon

A demonstration case, e.g. the 3-D flow past a flexible beam is to be found in the ParaSiF_CF/demo folder (https://github.com/ParaSiF-CF/ParaSiF_CF/tree/main/demo/flexBeam3D).

To run this demonstration case:

• Go to the demo folder and extract the mesh file
```bash
cd MESH && tar -xf Flex_Beam.tar.gz && cd ../
```

• Go to the run folder
```bash
cd runFolder/
```

• Source and export before run ParaSiF_CF cases (Please refer to the previous section).

• Compile Code_Saturne and copy FEniCS run case to the RESU sub-folder.
```bash
./SCRIPTS/compileCS
```

• Run the ParaSiF_CF framework.
```bash
./SCRIPTS/runcase
```
