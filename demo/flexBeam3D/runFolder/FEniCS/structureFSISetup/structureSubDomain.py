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
#    structureSubDomain.py
#
# Created
#    25-September-2019
#
# Description
#    This is the sub-domain class of the test case.
#
# Author
#    Wendi Liu
#
#---------------------------------------------------------------------------

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

from dolfin import *
import configparser

#_________________________________________________________________________________________
#
#%% Import configure file
#_________________________________________________________________________________________

config = configparser.ConfigParser()
config.read('./structureFSISetup/structureInputPara.ini')

OBeamXtemp=float(config['GEOMETRY']['OBeamX'])
OBeamYtemp=float(config['GEOMETRY']['OBeamY'])
OBeamZtemp=float(config['GEOMETRY']['OBeamZ'])
XBeamtemp=float(config['GEOMETRY']['XBeam'])
YBeamtemp=float(config['GEOMETRY']['YBeam'])
ZBeamtemp=float(config['GEOMETRY']['ZBeam'])

#_________________________________________________________________________________________
#
#%% Define SubDomains classes
#%% for defining parts of the boundaries and the interior of the domain
#_________________________________________________________________________________________

class Fixed( SubDomain ):
    def inside (self , x, on_boundary ):
        tol = DOLFIN_EPS
        return near(x[1], (OBeamYtemp + tol))

class Flex( SubDomain ):
    def inside (self , x, on_boundary ):
        tol = DOLFIN_EPS
        return near(x[1], (OBeamYtemp + YBeamtemp - tol)) or near(x[0], (OBeamXtemp + XBeamtemp - tol)) or near(x[0], (OBeamXtemp + tol)) or near(x[2], (OBeamZtemp + tol)) 

class Symmetry( SubDomain ):
    def inside (self , x, on_boundary ):
        tol = DOLFIN_EPS
        return near(x[2], (OBeamZtemp + ZBeamtemp - tol))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#