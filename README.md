# ParaSiF_CF - Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation Framework

Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation Framework employs Code_Saturne to solve the computational fluid dynamics (CFD), FEniCS to solve the computational structure mechanics (CSM) and MUI for data exchange. It offers a platform where users can carry out fluid-structure interaction studies using supercomputers.

The framework uses a partitioned approach. It takes several advantages of the MUI library:

• Good scalability on communications among domains for large simulations;

• Keeping the development of the coupled solvers decoupled. It allows for easier independent testing of each solver and avoids potential incompatibilities between two solvers (i.e if they both use a certain library X each one depends on a different version of it which are not compatible);

• Coupling of multiple solvers which have different programming language interfaces (e.g C, C++, FORTRAN and Python);

• "Plug and play" strategy. One solver can be replaced by another incrementally and without the need of recompiling if the MUI interface is used as a common adaptor;

• Use of multiple solvers which have two incompatible licenses exploiting the dual licensing of the MUI library (both solvers are never mixed source-wise or binary-wise).

**This framework is a beta software at the moment and under active development**.

## Licensing

Copyright (C) 2021 Engineering and Environment Group, Scientific Computing Department, Science and Technology Facilities Council, UK Research and Innovation. All rights reserved.

This code is licensed under the GNU General Public License version 3

## Acknowledgements
The Parallel Partitioned Multi-physical Simulation Framework is developed at the [Scientific Computing Department](https://www.scd.stfc.ac.uk/) of the [Science and Technology Facilities Council](https://stfc.ukri.org/). If you use this framework, please cite us:

*

## Publication


## Contact

Should you have any question please do not hesitate to contact the developers

## Installation

There is no need to install ParaSiF_CF itself, but the following should be done before ParaSiF_CF can be used:

• FEniCS v2019.0.1 should be installed;

• Code_Saturne v6.0.6 should be installed;

• MUI v1.1.2 (https://github.com/MxUI) should be obtained and its C & Python wrappers should be compiled;

• MUI_Utilities (https://github.com/MUI-Utilities) should be obtained;

• Compile the C & Python wrappers of the FSI Coupling Lab of the MUI Utility;

• Follow the Code_Saturne_MUI_Coupling library of the MUI_Utilities to establish the Code_Saturne - MUI coupling;

## Source and export before run parMupSiF cases

```bash
source /path/to/dolfin/dolfin.conf
source /path/to/Code_Saturne
export PYTHONPATH= /path/to/parMupSiF/src/CSM/FEniCS/V2019.1.0:$PYTHONPATH
```

## Demo

Demo of the single phase FSI case (3-D flow past a flexible beam) is in ParaSiF_CF/demo folder.
