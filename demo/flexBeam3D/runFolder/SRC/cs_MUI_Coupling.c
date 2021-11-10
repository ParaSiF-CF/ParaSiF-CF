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

#include "cs_MUI_Coupling.h"
#include "cs_halo.h"
#include "cs_parameters.h"
#include "cs_field_operator.h"
#include "cs_turbulence_model.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

mui_coupling cpl;

//-----

MPI_Comm* get_MPI_Comm()
{
    MPI_Comm* Pworld = &cs_glob_mpi_comm;
    
    return Pworld;
}

void setup_MUI_Coupling( const char* app_name )
{
    /// Declare MUI objects 
    const char* URI1 = "mpi://";
    const char* URI2 = app_name;
    const char* URI3 = "/threeDInterface0";
    float r    = 5.0;                      // search radius
    bft_printf("{cs} cs_MUI_Coupling.c\n");
    cpl.URI = malloc(strlen(URI1) + strlen(URI2) + strlen(URI3));
    
    strcpy(cpl.URI, URI1);
    strcat(cpl.URI, URI2);
    strcat(cpl.URI, URI3);
 
    cpl.uniface = mui_create_uniface_3d((const char*) cpl.URI );

    /// Define spatial and temporal samplers
    cpl.spatial = mui_create_sampler_pseudo_nearest2_linear_3d(r);
    cpl.temporal = mui_create_chrono_sampler_exact_3d(1e-37);

}

//-----
/// Declare MPI ranks
int cs_get_MUI_rank()
{

    int rank;
    MPI_Comm_rank( *get_MPI_Comm(), &rank );
    return rank;

}

//-----
/// Declare MPI ranks
int cs_get_MUI_size()
{

    int size;
    MPI_Comm_size( *get_MPI_Comm(), &size );
    return size;

}

//-----
/// Declare MPI ranks
int cs_local_MPI_barrier()
{

    int size;
    size = MPI_Barrier( *get_MPI_Comm());
    return size;

}

//-----

void mui_announce_span_send(double coord_min_sendX, 
                            double coord_min_sendY, 
                            double coord_min_sendZ, 
                            double coord_max_sendX, 
                            double coord_max_sendY, 
                            double coord_max_sendZ)
{

    if ((coord_min_sendX > coord_max_sendX)||(coord_min_sendY > coord_max_sendY)||(coord_min_sendZ > coord_max_sendZ)){

        mui_announce_send_disable_3d(cpl.uniface);

    } else {

        mui_announce_send_span_3d_box(cpl.uniface, coord_min_sendX,coord_min_sendY,coord_min_sendZ,
                                            coord_max_sendX,coord_max_sendY,coord_max_sendZ, 
                                            -1,(cs_glob_time_step->nt_max)*1e+5);
        
		int my_rank;
		my_rank=cs_get_MUI_rank();
		printf("{CS} rank contains beam: %d \n",my_rank);
        printf("{CS} coord_send: %f %f %f to %f %f %f from %d with time to %d \n",coord_min_sendX,coord_min_sendY,coord_min_sendZ,coord_max_sendX,coord_max_sendY,coord_max_sendZ, cs_get_MUI_rank(), (cs_glob_time_step->nt_max));

    }
    bft_printf("AFTER send creating box \n");

}
//-----

void mui_announce_span_rcv(double coord_min_rcvX, 
                            double coord_min_rcvY, 
                            double coord_min_rcvZ, 
                            double coord_max_rcvX, 
                            double coord_max_rcvY, 
                            double coord_max_rcvZ)
{

    if ((coord_min_rcvX > coord_max_rcvX)||(coord_min_rcvY > coord_max_rcvY)||(coord_min_rcvZ > coord_max_rcvZ)){

        mui_announce_recv_disable_3d(cpl.uniface);

    } else {

        mui_announce_recv_span_3d_box(cpl.uniface, coord_min_rcvX,coord_min_rcvY,coord_min_rcvZ,
                                            coord_max_rcvX,coord_max_rcvY,coord_max_rcvZ, 
                                            -1,(cs_glob_time_step->nt_max)*1e+5);
        printf("{CS} coord_rcv: %f %f %f to %f %f %f from %d with time to %d \n",coord_min_rcvX,coord_min_rcvY,coord_min_rcvZ,coord_max_rcvX,coord_max_rcvY,coord_max_rcvZ, cs_get_MUI_rank(), (cs_glob_time_step->nt_max));

    }
    bft_printf("AFTER rcv creating box \n");

}


//-----

mui_uniface_3d* cs_get_uniface()
{
 
    return cpl.uniface;
 
}

//-----

mui_sampler_pseudo_nearest2_linear_3d* cs_get_spatial()
{
 
    return cpl.spatial;
 
}

//-----

mui_chrono_sampler_exact_3d* cs_get_temporal()
{
 
    return cpl.temporal;
 
}

double cs_fetch_disp_MUI_Coupling(    const char*                fetch_name,
                                     double                     x_coord,
                                    double                     y_coord,
                                    double                     z_coord,
                                    int                   sub_iteration_numbers_MUI_Coupling,
                                    int                   current_sub_iteration_number,
                                    bool                   parallel_FSI_coupling) 
{
    char* name = malloc(strlen(fetch_name) + 1);
    strcpy(name, fetch_name);

    int iterations_current = -999;

    double disp_comp = 0.0;

    if (parallel_FSI_coupling)
    {

        iterations_current = ((cs_glob_time_step->nt_cur) * sub_iteration_numbers_MUI_Coupling  + current_sub_iteration_number )-1;

    }else
    {

        iterations_current = ((cs_glob_time_step->nt_cur) * sub_iteration_numbers_MUI_Coupling  + current_sub_iteration_number );

    }

    if (iterations_current > 0)
    {
        mui_point_3d fetch_point3d = { x_coord, y_coord, z_coord };
 
        disp_comp = mui_fetch_pseudo_nearest_neighbor2_linear_exact_3d( cpl.uniface,
                                                                    name,
                                                                    fetch_point3d,
                                                                    iterations_current,
                                                                    cpl.spatial,
                                                                    cpl.temporal );

    }else
    {
        disp_comp = 0.0;
    }
    
    return disp_comp;
}
//-----

void commit_MUI_Coupling(const int    sub_iteration_numbers_MUI_Coupling,
                         const int    current_sub_iteration_number)
{
    int iterations_current = ((cs_glob_time_step->nt_cur ) * sub_iteration_numbers_MUI_Coupling + current_sub_iteration_number);

    int iterations_forget = iterations_current -4;
    int cr = mui_commit_rank_3d( cpl.uniface, iterations_current );
    mui_forget_upper_3d( cpl.uniface, iterations_forget, 1 );
	mui_set_forget_length_3d( cpl.uniface, 4 );

    int my_rank;
		my_rank=cs_get_MUI_rank();
	printf("{CS} commit Step %d with commit peers %d at rank %d \n", iterations_current, cr, my_rank);
} 

//-----

void commit_MUI_Zero()
{
    int cr = mui_commit_rank_3d( cpl.uniface, 0 );
    bft_printf("{CS} commit number 0 with commit ranks %d \n", cr);
}

//----

void barrier_MUICpl(const int    sub_iteration_numbers_MUI_Coupling,
                         const int    current_sub_iteration_number)
{
    int iterations_current = ((cs_glob_time_step->nt_cur ) * sub_iteration_numbers_MUI_Coupling + current_sub_iteration_number);

    mui_barrier_3d( cpl.uniface, iterations_current );

}            

//-----

void cs_push_MUI_Coupling(const char*      field_name,
                          double           x_coord,
                          double           y_coord,
                          double           z_coord,
                          double           push_value)
{

    char* name = malloc(strlen(field_name) + 1);
    strcpy(name, field_name);

    mui_point_3d push_point3d = { x_coord, y_coord, z_coord };
    mui_push_3d(   cpl.uniface,
                name,
                push_point3d,
                push_value );

}

//-----

void cs_push_field_MUI_Coupling(const char*      field_name,
                                char*            donate_type,
                                char*            donate) 
{
    cs_field_t *f = cs_field_by_name(field_name);
   
    char* donate_boundary_faces;
    char* donate_cells;
    donate_boundary_faces = malloc( 100 );
    donate_cells = malloc( 100 );
    strcpy(donate_boundary_faces,"Boundary_Faces");
    strcpy(donate_cells,"Cells");

    if(strcmp(donate_type, donate_boundary_faces) == 0 ) {
    
        cs_lnum_t *face_list = NULL;
        cs_lnum_t n_faces=0;

        BFT_MALLOC(face_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

        char* selection_donate;
        selection_donate = malloc( 100 );
        strcpy( selection_donate, donate );

        cs_selector_get_b_face_list(selection_donate,
                                    &n_faces,
                                    face_list);
    
        for (cs_lnum_t j = 0; j < f->dim; j++) {
            char* name = malloc(strlen(field_name) + 1);
            strcpy(name, field_name);
            sprintf(&name[strlen(field_name)], "%i", j);

            for (cs_lnum_t i = 0; i < n_faces; i++) {
            
                cs_lnum_t face_id = face_list[i];
                cs_lnum_t cell_id = cs_glob_mesh->b_face_cells[face_id];

                mui_point_3d push_point3d = { cs_glob_mesh_quantities->b_face_cog[face_list[i]*3+0],
                            cs_glob_mesh_quantities->b_face_cog[face_list[i]*3+1],
                            cs_glob_mesh_quantities->b_face_cog[face_list[i]*3+2]};

                mui_push_3d(   cpl.uniface,
                            name,
                            push_point3d,
                            ( double ) f->val[cell_id]*f->dim+j );

                if(i == (n_faces-3)) {
                    bft_printf("{CS} %s is %lf \n", field_name, f->val[cell_id]*f->dim+j);
                    bft_printf("{CS} X face is %lf \n", cs_glob_mesh_quantities->b_face_cog[face_list[i]*3+0]);
                    bft_printf("{CS} Y face is %lf \n", cs_glob_mesh_quantities->b_face_cog[face_list[i]*3+1]);
                    bft_printf("{CS} Z face is %lf \n", cs_glob_mesh_quantities->b_face_cog[face_list[i]*3+2]);
                    bft_printf("{CS} X cell is %lf \n", cs_glob_mesh_quantities->cell_cen[cell_id*3+0]);
                    bft_printf("{CS} Y cell is %lf \n", cs_glob_mesh_quantities->cell_cen[cell_id*3+1]);
                    bft_printf("{CS} Z cell is %lf \n", cs_glob_mesh_quantities->cell_cen[cell_id*3+2]);
                 
                }

            }

            free( name );
        }
    
    BFT_FREE(face_list);
    } else if(strcmp(donate_type, donate_cells) == 0 ) {
    
        cs_lnum_t *elt_list = NULL;
        int n_elts;

        BFT_MALLOC(elt_list, cs_glob_mesh->n_cells, cs_lnum_t);

        char* selection_donate;
        selection_donate = malloc( 100 );
        strcpy( selection_donate, donate );

        cs_selector_get_cell_num_list(selection_donate,
                                        &n_elts,
                                        elt_list);
    
        for (cs_lnum_t j = 0; j < f->dim; j++) {
            char* name = malloc(strlen(field_name) + 1);
            strcpy(name, field_name);
            sprintf(&name[strlen(field_name)], "%i", j);

            for (cs_lnum_t i = 0; i < n_elts; i++) {
                mui_point_3d push_point3d = { cs_glob_mesh_quantities->cell_cen[elt_list[i]*3+0],
                            cs_glob_mesh_quantities->cell_cen[elt_list[i]*3+1],
                            cs_glob_mesh_quantities->cell_cen[elt_list[i]*3+2]};
                mui_push_3d(   cpl.uniface,
                            name,
                            push_point3d,
                            ( double ) f->val[elt_list[i]*f->dim+j] );

            }

            free( name );
        }
    
    BFT_FREE(elt_list);
    } else {

        bft_printf("Unrecognized donate type (%s) \n Please set the donate type as '%s' or '%s' \n", donate_type, donate_boundary_faces, donate_cells);
        exit(1);

    }
}

