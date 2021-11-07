/*---------------------------------------------------------------------------*\
    Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation 
	Framework (ParaSiF_CF) employs Code_Saturne to solve the computational 
	fluid dynamics (CFD), FEniCS to solve the computational structure mechanics 
	(CSM) and MUI for data exchange.
    Copyright (C) 2021 Engineering and Environment Group, Scientific 
    Computing Department, Science and Technology Facilities Council, 
    UK Research and Innovation. All rights reserved.
    This code is licensed under the GNU General Public License version 3
    ** GNU General Public License, version 3 **
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------

Filename 
    cs_MUI_Coupling.h

Created
    05 January 2021

Author
    Wendi Liu; Alex Skillen

\*---------------------------------------------------------------------------*/

#ifndef CS_MUI_COUPLING_H
#define CS_MUI_COUPLING_H

#include "wrappers/C/mui_c_wrapper_3d.h"

#include "bft_mem.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_selector.h"
#include "cs_field.h"
#include "cs_time_step.h"


typedef struct 
{
    char* URI;

    mui_uniface_3d* uniface;

    mui_sampler_pseudo_nearest2_linear_3d* spatial;
//    mui_sampler_nearest3d* spatial;
    mui_chrono_sampler_exact_3d* temporal;

} mui_coupling; 

extern mui_coupling cpl;

mui_uniface_3d* cs_get_uniface();
mui_sampler_pseudo_nearest2_linear_3d* cs_get_spatial();
mui_chrono_sampler_exact_3d* cs_get_temporal();

void setup_MUI_Coupling( const char* );

double cs_fetch_disp_MUI_Coupling(    const char*                fetch_name,
                                    double                     x_coord,
                                    double                     y_coord,
                                    double                     z_coord,
                                    int                        sub_iteration_numbers_MUI_Coupling,
                                    int                        current_sub_iteration_number,
                                    bool                       parallel_FSI_coupling);

int cs_get_MUI_rank();
int cs_get_MUI_size();
int cs_local_MPI_barrier();
void mui_announce_span_send(double coord_min_sendX, 
							double coord_min_sendY, 
							double coord_min_sendZ, 
							double coord_max_sendX, 
							double coord_max_sendY, 
							double coord_max_sendZ);
void mui_announce_span_rcv(double coord_min_rcvX, 
							double coord_min_rcvY, 
							double coord_min_rcvZ, 
							double coord_max_rcvX, 
							double coord_max_rcvY, 
							double coord_max_rcvZ);
void commit_MUI_Coupling(const int    sub_iteration_numbers_MUI_Coupling,
                         const int    current_sub_iteration_number);
void commit_MUI_Zero();

void barrier_MUICpl(const int    sub_iteration_numbers_MUI_Coupling,
                    const int    current_sub_iteration_number);

void cs_push_MUI_Coupling(const char*    field_name,
                          double         x_coord,
                          double         y_coord,
                          double         z_coord,
                          double         push_value);

void cs_push_field_MUI_Coupling(const char*    field_name,
                                char*          donate_type,
                                char*          donate);
#endif