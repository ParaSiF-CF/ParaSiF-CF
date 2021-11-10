/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* Code_Saturne version 6.0.5 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_geom.h"
#include "cs_interpolate.h"
#include "cs_mesh.h"
#include "cs_selector.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_probe.h"
#include "cs_time_plot.h"

#include "cs_field_pointer.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void)
{
 /* redefine profile output options such as format, frequency, etc... */

  cs_post_define_writer(CS_POST_WRITER_PROFILES,  /* writer_id */
                        "",                       /* writer name */
                        "profiles",
                        "plot",                   /* format name */
                        "dat",                    /* format options */
                        FVM_WRITER_TRANSIENT_COORDS,
                        true,                    /* output_at_start */
                        true,                     /* output at end */
                        400,                       /* time step frequency */
                        -1.0);                    /* time value frequency */
                        
                        
                        
/* Set time plot file writer flush behavior defaults. */

  /*! [post_set_tp_flush] */
  cs_time_plot_set_flush_default(18000, /* flush_wtime */
                                 -1);  /* n_buffer_steps */
  /*! [post_set_tp_flush] */

  /* Default output format and options */

  /* Redefine default writer */
  /* ----------------------- */

  /*! [post_define_writer_m1] */
  cs_post_define_writer(CS_POST_WRITER_DEFAULT,       /* writer_id */
                        "results",                    /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format_name */
                        "",                           /* format_options */
                        FVM_WRITER_TRANSIENT_COORDS,
                        true,                        /* output_at_start */
                        true,                         /* output_at_end */
                        400,                           /* frequency_n */
                        -1.0);                        /* frequency_t */
  /*! [post_define_writer_m1] */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define monitoring probes and profiles.
 *
 * Profiles are defined as sets of probes.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_probes(void)
{
/* A writer (id = CS_POST_WRITER_PROBES) using the format "time_plot" is
* associated by default to a set of monitoring probes.
* This is not the case for a profile. */
{
cs_probe_set_t *pset = cs_probe_set_create("Monitoring");
cs_probe_set_add_probe(pset, 0.45, 0.15, -0.15, "Point_A");
}
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, NULL otherwise
 * \param[in]       n_cells      local number of cells of post_mesh
 * \param[in]       n_i_faces    local number of interior faces of post_mesh
 * \param[in]       n_b_faces    local number of boundary faces of post_mesh
 * \param[in]       n_vertices   local number of vertices faces of post_mesh
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing
 *                               mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertices,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts)
{

  if (cat_id == CS_POST_MESH_VOLUME) {

    const cs_real_t  *volume = cs_glob_mesh_quantities->cell_vol;

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "Cell Volume",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      volume,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

  }
  /*< [postprocess_values_ex_1] */

}

/*----------------------------------------------------------------------------*/
/*!
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * \param  nt_max_abs  maximum time step number
 * \param  nt_cur_abs  current time step number
 * \param  t_cur_abs   absolute time at the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
