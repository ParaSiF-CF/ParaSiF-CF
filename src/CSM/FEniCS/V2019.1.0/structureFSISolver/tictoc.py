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
#    tictoc.py
#
# Created
#    25-September-2019
#
# Description
#    This is the wall time class of the structure code.
#
# Author
#    Wendi Liu
#
#---------------------------------------------------------------------------

#_________________________________________________________________________________________
#
#%% Import packages
#_________________________________________________________________________________________

import time

#_________________________________________________________________________________________
#
#%% Wall clock class
#_________________________________________________________________________________________

class TicToc:
    # TicToc wall clock class
    def __init__(self):
        self.t0 = 0.0 # initial time
        self.t1 = 0.0 # final time

    # Records the beginning of a time interval
    def tic(self):
        self.t1 = time.time()
        self.t0 = self.t1

    # Records the ending of a time interval and return the time difference
    def toc(self):
        self.t1 = time.time()
        deltaT = self.t1 - self.t0 
        return deltaT

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILE END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#