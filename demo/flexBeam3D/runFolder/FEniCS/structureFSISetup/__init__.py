#---------------------------------------------------------------------------
#    Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation 
#    Framework (ParaSiF_CF) employs Code_Saturne to solve the computational 
#    fluid dynamics (CFD), FEniCS to solve the computational structure mechanics 
#    (CSM) and MUI for data exchange.
#    Copyright (C) 2021 Engineering and Environment Group, Scientific 
#    Computing Department, Science and Technology Facilities Council, 
#    UK Research and Innovation. All rights reserved.
#    This code is licensed under the GNU General Public License version 3
#    ** GNU General Public License, version 3 **
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------
#
# Filename 
#    __init__.py
#
# Created
#    25-September-2019
#
# Description
#    This is the dundant init file of the test case.
#
# Author
#    Wendi Liu
#
#---------------------------------------------------------------------------

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from structureFSISetup.structureSubDomain import Fixed
from structureFSISetup.structureSubDomain import Flex
from structureFSISetup.structureSubDomain import Symmetry
from structureFSISetup.structureBCS import boundaryConditions

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#