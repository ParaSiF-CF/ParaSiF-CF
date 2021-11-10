!-------------------------------------------------------------------------------

!                      Code_Saturne version 6.0.5
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_modules.f90
!>
!> \brief User-defined module: it allows to create any user array.
!>
!> See \subpage cs_user_modules for examples.
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!    Parallel Partitioned Fluid-Structure Interaction (FSI) Simulation 
!    Framework (ParaSiF_CF) employs Code_Saturne to solve the computational 
!    fluid dynamics (CFD), FEniCS to solve the computational structure mechanics 
!    (CSM) and MUI for data exchange.
!    Copyright (C) 2021 Engineering and Environment Group, Scientific 
!    Computing Department, Science and Technology Facilities Council, 
!    UK Research and Innovation. All rights reserved.
!    This code is licensed under the GNU General Public License version 3
!    ** GNU General Public License, version 3 **
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------

!> \Filename      cs_user_modules.f90
!> \Created       05 January 2019
!> \Author        Wendi Liu; Alex Skillen
!______________________________________________________________________________


module cs_user_mui_coupling

implicit none

    character (len = *), parameter :: push_MUI_name="8"
    character (len = *), parameter :: fetch_MUI_name="8"
    logical, parameter    ::    push_force_MUI = .true.
    logical, parameter    ::    push_in_C = .false.
    logical(kind=1), parameter    ::    parallel_FSI_coupling = .false.
    integer, parameter    ::    sub_iteration_MUI_Coupling = 15
    integer, parameter    ::    muiCouplingMethod = 2
    double precision, parameter    ::    init_und_relx_coupling = 0.8
    integer, parameter    ::    aitkenIterationN_IQNILS = 15
    double precision, parameter    ::    und_relx_coupling_Max = 0.5
    logical(kind=1)    ::    local_push = .true.
    logical(kind=1)    ::    local_fetch = .true.

end module cs_user_mui_coupling

!!========================================================================

module cs_commons_mui_coupling

implicit none

    double precision, dimension(:,:), allocatable :: fetchDataArray
    double precision, dimension(:), allocatable :: fetchIndexArray

end module cs_commons_mui_coupling

!!========================================================================

module cs_mui_coupling

    !===============================================================================
    !===============================================================================
    ! Module files
    !===============================================================================
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, &
      c_int, c_bool, c_double
    use entsor
    !===============================================================================

    implicit none
    type(c_ptr) :: C_cs_get_uniface, C_cs_get_spatial, C_cs_get_temporal
    type(c_ptr) :: cs_get_uniface, cs_get_spatial, cs_get_temporal
    type(c_ptr) :: C_get_MPI_Comm, get_MPI_Comm
    interface
        !-----
        !-----
        function C_get_MPI_Comm() bind(C, name="get_MPI_Comm") 
            use iso_c_binding 
        end function 

        function C_cs_get_uniface() &
          bind(C, name='cs_get_uniface')
          use, intrinsic :: iso_c_binding
          implicit none
        end function C_cs_get_uniface
        
        !-----
        
        function C_cs_get_spatial() &
          bind(C, name='cs_get_spatial')
          use, intrinsic :: iso_c_binding
          implicit none
        end function C_cs_get_spatial
        
        !-----
        
        function C_cs_get_temporal() &
          bind(C, name='cs_get_temporal')
          use, intrinsic :: iso_c_binding
          implicit none
        end function C_cs_get_temporal

        !-----

        integer(kind=c_int) function C_cs_get_MUI_rank() &
          bind(C, name='cs_get_MUI_rank')
          use, intrinsic :: iso_c_binding
          implicit none
        end function C_cs_get_MUI_rank

        !-----

        integer(kind=c_int) function C_cs_get_MUI_size() &
          bind(C, name='cs_get_MUI_size')
          use, intrinsic :: iso_c_binding
          implicit none
        end function C_cs_get_MUI_size

        !-----

        integer(kind=c_int) function C_cs_local_MPI_barrier() &
          bind(C, name='cs_local_MPI_barrier')
          use, intrinsic :: iso_c_binding
          implicit none
        end function C_cs_local_MPI_barrier

        
        real (kind=c_double) function C_cs_fetch_disp_MUI_Coupling(  fetch_name,  &
                                                                     x_coord,     &
                                                                     y_coord,     &
                                                                     z_coord,     &
                                                                     sub_iteration_numbers_MUI_Coupling, &
                                                                     current_sub_iteration_number, &
                                                                     parallel_FSI_coupling)       &
          bind(C, name='cs_fetch_disp_MUI_Coupling')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind=c_char, len=1), dimension(*), intent(in) :: fetch_name
          real (kind=c_double), VALUE :: x_coord
          real (kind=c_double), VALUE :: y_coord
          real (kind=c_double), VALUE :: z_coord
          integer(c_int), VALUE :: sub_iteration_numbers_MUI_Coupling
          integer(c_int), VALUE :: current_sub_iteration_number
          LOGICAL(KIND=C_BOOL), VALUE :: parallel_FSI_coupling
        end function C_cs_fetch_disp_MUI_Coupling

        subroutine C_mui_announce_span_send(coord_min_sendX,     &
                                            coord_min_sendY,     &
                                            coord_min_sendZ,     &
                                            coord_max_sendX,     &
                                            coord_max_sendY,     &
                                            coord_max_sendZ) &
          bind(C, name='mui_announce_span_send')
          use, intrinsic :: iso_c_binding
          implicit none
          real (kind=c_double), VALUE :: coord_min_sendX
          real (kind=c_double), VALUE :: coord_min_sendY
          real (kind=c_double), VALUE :: coord_min_sendZ
          real (kind=c_double), VALUE :: coord_max_sendX
          real (kind=c_double), VALUE :: coord_max_sendY
          real (kind=c_double), VALUE :: coord_max_sendZ
        end subroutine C_mui_announce_span_send
    
        subroutine C_mui_announce_span_rcv(coord_min_rcvX,     &
                                           coord_min_rcvY,     &
                                           coord_min_rcvZ,     &
                                           coord_max_rcvX,     &
                                           coord_max_rcvY,     &
                                           coord_max_rcvZ) &

          bind(C, name='mui_announce_span_rcv')
          use, intrinsic :: iso_c_binding
          implicit none
          real (kind=c_double), VALUE :: coord_min_rcvX
          real (kind=c_double), VALUE :: coord_min_rcvY
          real (kind=c_double), VALUE :: coord_min_rcvZ
          real (kind=c_double), VALUE :: coord_max_rcvX
          real (kind=c_double), VALUE :: coord_max_rcvY
          real (kind=c_double), VALUE :: coord_max_rcvZ
        end subroutine C_mui_announce_span_rcv

        subroutine C_commit_MUI_Coupling(sub_iteration_numbers_MUI_Coupling, current_sub_iteration_number) &
          bind(C, name='commit_MUI_Coupling')
          use, intrinsic :: iso_c_binding
          implicit none
          integer(c_int), VALUE :: sub_iteration_numbers_MUI_Coupling
          integer(c_int), VALUE :: current_sub_iteration_number
        end subroutine C_commit_MUI_Coupling

        subroutine C_commit_MUI_Zero() &
          bind(C, name='commit_MUI_Zero')
          use, intrinsic :: iso_c_binding
          implicit none
        end subroutine C_commit_MUI_Zero

        subroutine C_barrier_MUICpl(sub_iteration_numbers_MUI_Coupling, current_sub_iteration_number) &
          bind(C, name='barrier_MUICpl')
          use, intrinsic :: iso_c_binding
          implicit none
          integer(c_int), VALUE :: sub_iteration_numbers_MUI_Coupling
          integer(c_int), VALUE :: current_sub_iteration_number
        end subroutine C_barrier_MUICpl
        
        subroutine C_cs_push_MUI_Coupling(field_name, x_coord, y_coord, z_coord, push_value) &
          bind(C, name='cs_push_MUI_Coupling')
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char, len=1), dimension(*), intent(in) :: field_name
            real (kind=c_double), VALUE :: x_coord
            real (kind=c_double), VALUE :: y_coord
            real (kind=c_double), VALUE :: z_coord
            real (kind=c_double), VALUE :: push_value
        end subroutine C_cs_push_MUI_Coupling

        subroutine C_cs_push_field_MUI_Coupling(field_name, donate_type, donate) &
          bind(C, name='cs_push_field_MUI_Coupling')
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char, len=1), dimension(*), intent(in) :: field_name
            character(kind=c_char, len=1), dimension(*), intent(in) :: donate_type
            character(kind=c_char, len=1), dimension(*), intent(in) :: donate
        end subroutine C_cs_push_field_MUI_Coupling

    end interface

    contains
    !===============================================================================
    
    subroutine cs_push_MUI_Coupling(field_name, x_coord, y_coord, z_coord, push_value)
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        character(len=*), intent(in) :: field_name

        ! Local variables
        character(len=len_trim(field_name)+1, kind=c_char) :: c_field_name
        real (kind=c_double), VALUE :: x_coord
        real (kind=c_double), VALUE :: y_coord
        real (kind=c_double), VALUE :: z_coord
        real (kind=c_double), VALUE :: push_value

        c_field_name = trim(field_name)//c_null_char
        
        call C_cs_push_MUI_Coupling(c_field_name, x_coord, y_coord, z_coord, push_value)
    end subroutine cs_push_MUI_Coupling
    
    !===============================================================================

    function get_MPI_Comm()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        get_MPI_Comm = C_get_MPI_Comm()
        write(nfecra, *) "FORTRAN get_MPI_Comm"
        
        return
        
    end function get_MPI_Comm

    function cs_get_uniface()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        cs_get_uniface = C_cs_get_uniface()
        
        return
        
    end function cs_get_uniface

    function cs_get_spatial()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        cs_get_spatial = C_cs_get_spatial()
        
        return
        
    end function cs_get_spatial

    function cs_get_temporal()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        cs_get_temporal = C_cs_get_temporal()
        
        return
        
    end function cs_get_temporal

    integer(kind=c_int) function cs_get_MUI_rank()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        cs_get_MUI_rank = C_cs_get_MUI_rank()
        
        return
        
    end function cs_get_MUI_rank

    integer(kind=c_int) function cs_get_MUI_size()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        cs_get_MUI_size = C_cs_get_MUI_size()
        
        return
        
    end function cs_get_MUI_size

    integer(kind=c_int) function cs_local_MPI_barrier()
        
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        
        ! Local variables

        cs_local_MPI_barrier = C_cs_local_MPI_barrier()
        
        return
        
    end function cs_local_MPI_barrier

    real (kind=c_double) function cs_fetch_disp_MUI_Coupling( fetch_name,  &
                                                              x_coord,     &
                                                              y_coord,     &
                                                              z_coord,     &
                                                              sub_iteration_numbers_MUI_Coupling, &
                                                              current_sub_iteration_number)
        use, intrinsic :: iso_c_binding
        use cs_user_mui_coupling
        implicit none

        ! Arguments
        character(len=*), intent(in) :: fetch_name 
        real (kind=c_double), VALUE :: x_coord
        real (kind=c_double), VALUE :: y_coord
        real (kind=c_double), VALUE :: z_coord
        integer(c_int), VALUE :: sub_iteration_numbers_MUI_Coupling
        integer(c_int), VALUE :: current_sub_iteration_number

        ! Local variables
        character(len=len_trim(fetch_name)+1, kind=c_char) :: c_fetch_name

        c_fetch_name = trim(fetch_name)//c_null_char

        cs_fetch_disp_MUI_Coupling = C_cs_fetch_disp_MUI_Coupling( c_fetch_name, &
                                                                   x_coord,     &
                                                                   y_coord,     &
                                                                   z_coord,     &
                                                                   sub_iteration_numbers_MUI_Coupling, &
                                                                   current_sub_iteration_number, &
                                                                   parallel_FSI_coupling) 
        return
    end function cs_fetch_disp_MUI_Coupling

    !===============================================================================

    subroutine mui_announce_span_send(coord_min_sendX,     &
                                      coord_min_sendY,     &
                                      coord_min_sendZ,     &
                                      coord_max_sendX,     &
                                      coord_max_sendY,     &
                                      coord_max_sendZ)

        use, intrinsic :: iso_c_binding
        implicit none
        ! Arguments
        real (kind=c_double), VALUE :: coord_min_sendX
        real (kind=c_double), VALUE :: coord_min_sendY
        real (kind=c_double), VALUE :: coord_min_sendZ
        real (kind=c_double), VALUE :: coord_max_sendX
        real (kind=c_double), VALUE :: coord_max_sendY
        real (kind=c_double), VALUE :: coord_max_sendZ

        call C_mui_announce_span_send(coord_min_sendX,     &
                                      coord_min_sendY,     &
                                      coord_min_sendZ,     &
                                      coord_max_sendX,     &
                                      coord_max_sendY,     &
                                      coord_max_sendZ)

    end subroutine mui_announce_span_send    
    !===============================================================================

    subroutine mui_announce_span_rcv(coord_min_rcvX,     &
                                      coord_min_rcvY,     &
                                      coord_min_rcvZ,     &
                                      coord_max_rcvX,     &
                                      coord_max_rcvY,     &
                                      coord_max_rcvZ)

        use, intrinsic :: iso_c_binding
        implicit none
        ! Arguments
        real (kind=c_double), VALUE :: coord_min_rcvX
        real (kind=c_double), VALUE :: coord_min_rcvY
        real (kind=c_double), VALUE :: coord_min_rcvZ
        real (kind=c_double), VALUE :: coord_max_rcvX
        real (kind=c_double), VALUE :: coord_max_rcvY
        real (kind=c_double), VALUE :: coord_max_rcvZ

        call C_mui_announce_span_rcv(coord_min_rcvX,     &
                                      coord_min_rcvY,     &
                                      coord_min_rcvZ,     &
                                      coord_max_rcvX,     &
                                      coord_max_rcvY,     &
                                      coord_max_rcvZ)

    end subroutine mui_announce_span_rcv

    !===============================================================================
    subroutine commit_MUI_Coupling(sub_iteration_numbers_MUI_Coupling, current_sub_iteration_number)
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        integer(kind=c_int), VALUE :: sub_iteration_numbers_MUI_Coupling
        integer(kind=c_int), VALUE :: current_sub_iteration_number

        call C_commit_MUI_Coupling(sub_iteration_numbers_MUI_Coupling, current_sub_iteration_number)
    end subroutine commit_MUI_Coupling
    !===============================================================================
        
    subroutine commit_MUI_Zero()
        use, intrinsic :: iso_c_binding
        implicit none

        call C_commit_MUI_Zero()
    end subroutine commit_MUI_Zero
    !===============================================================================    
    subroutine cs_push_field_MUI_Coupling(field_name, donate_type, donate)
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        character(len=*), intent(in) :: field_name
        character(len=*), intent(in) :: donate_type
        character(len=*), intent(in) :: donate

        ! Local variables
        character(len=len_trim(field_name)+1, kind=c_char) :: c_field_name
        character(len=len_trim(donate_type)+1, kind=c_char) :: c_donate_type
        character(len=len_trim(donate)+1, kind=c_char) :: c_donate

        c_field_name = trim(field_name)//c_null_char
        c_donate_type = trim(donate_type)//c_null_char
        c_donate = trim(donate)//c_null_char

        call C_cs_push_field_MUI_Coupling(c_field_name, c_donate_type, c_donate)
    end subroutine cs_push_field_MUI_Coupling
    !===============================================================================
    subroutine barrier_MUICpl(sub_iteration_numbers_MUI_Coupling, current_sub_iteration_number)
        use, intrinsic :: iso_c_binding
        implicit none

        ! Arguments
        integer(kind=c_int), VALUE :: sub_iteration_numbers_MUI_Coupling
        integer(kind=c_int), VALUE :: current_sub_iteration_number

        call C_barrier_MUICpl(sub_iteration_numbers_MUI_Coupling, current_sub_iteration_number)
    end subroutine barrier_MUICpl
    !===============================================================================


end module cs_mui_coupling

!!========================================================================
module cs_push_force_MUI_Coupling_Module
!===============================================================================
! MUI Coupling files
use cs_mui_coupling
use cs_user_mui_coupling
! MUI Coupling files end
!===============================================================================

implicit none
    contains
    subroutine cs_push_force_MUI_Coupling ( )

        use numvar
        use optcal
        use parall
        use mesh
        use cs_c_bindings
        use cs_user_mui_coupling

        implicit none

        integer          ifac
        integer          ii
        integer          ilelt  , nlelt, iii
        
        double precision forces_f_cell(3), traction_f_cell(3), moment_f_obj(3), face_area, temp_x, temp_y, temp_z
        double precision, allocatable, dimension(:,:) :: cdgfbo0
        save cdgfbo0
        double precision, dimension(:,:), pointer :: bfprp_for

        double precision xnod(3)
        integer, allocatable, dimension(:) :: lstelt

        real(kind=8) :: zero_real_8 = 0.0, test_value_1, test_value_2, test_value_3

        double precision tempDelete

        !===============================================================================

        if (iforbr.ge.0) call field_get_val_v(iforbr, bfprp_for)

        ! Allocate a temporary array for cells or interior/boundary faces selection
        allocate(lstelt(max(ncel,nfac,nfabor)))


        ! store the coordinates of the centers of the boundary faces for the initial mesh at time zero
        if (ttcabs.eq.0.0d0) then

            allocate(cdgfbo0(3,max(ncel,nfac,nfabor)))

            ! get the faces corresponding to the structure surface
            call getfbr(push_MUI_name, nlelt, lstelt)
            !==========

            ! loop over the boundary faces to obtain the coordinates of the centers of the boundary faces
            do ilelt = 1, nlelt
                ifac = lstelt(ilelt)
            
                select case (ndim)

                case (3)
                    ! obtain and store the global coordinate of the centers of the present boundary face
                    cdgfbo0(1, ifac) = cdgfbo(1, ifac)
                    cdgfbo0(2, ifac) = cdgfbo(2, ifac)
                    cdgfbo0(3, ifac) = cdgfbo(3, ifac)

                case (2)
                    ! obtain and store the global coordinate of the centers of the present boundary face
                    cdgfbo0(1, ifac) = cdgfbo(1, ifac)
                    cdgfbo0(2, ifac) = cdgfbo(2, ifac)
                    cdgfbo0(3, ifac) = 0.0
                    
                case (1)
                    ! obtain and store the global coordinate of the centers of the present boundary face
                    cdgfbo0(1, ifac) = cdgfbo(1, ifac)
                    cdgfbo0(2, ifac) = 0.0
                    cdgfbo0(3, ifac) = 0.0

                case default
            
                    write(nfecra, *) "{CS} Error number of ndim: ", ndim, ". Please set it as '3', '2' or '1'."
                
                    return 1
            
                end select

            enddo
                
        else if (ttcabs.gt.0.0d0) then

            if (iforbr.ge.0) then

                ! set the force and moment components to zero
                do ii =  1, ndim
                    traction_f_cell(ii) = 0.d0
                    moment_f_obj(ii) = 0.d0
                enddo

                ! get the faces corresponding to the structure surface
                call getfbr(push_MUI_name, nlelt, lstelt)
                !==========
                iii = 0
                ! loop over the boundary faces to calculate the face area and force components
                do ilelt = 1, nlelt
                    ifac = lstelt(ilelt)
                    iii = iii + 1
                    select case (ndim)

                    case (3)
            
                        ! calculate the area of the present boundary face with the unit of [m^2]
                        face_area = sqrt(surfbo(1,ifac)**2+surfbo(2,ifac)**2+surfbo(3,ifac)**2)
            
                        ! obtain the force components of the present boundary face with the unit of [N]
                        forces_f_cell(1) = bfprp_for(1, ifac)
                        forces_f_cell(2) = bfprp_for(2, ifac)
                        forces_f_cell(3) = bfprp_for(3, ifac)

                        ! obtain the traction components of the present boundary face with the unit of [N/m^2]
                        traction_f_cell(1) = forces_f_cell(1)/face_area
                        traction_f_cell(2) = forces_f_cell(2)/face_area
                        traction_f_cell(3) = forces_f_cell(3)/face_area

                        ! obtain the global coordinate of the centers of the present boundary face for the moved mesh 
                        xnod(1) = cdgfbo(1, ifac)
                        xnod(2) = cdgfbo(2, ifac)
                        xnod(3) = cdgfbo(3, ifac)

                        !Global Moments
                        moment_f_obj(1)= moment_f_obj(1)+ (bfprp_for(3, ifac)* xnod(2) - bfprp_for(2, ifac)*xnod(3))
                        moment_f_obj(2)= moment_f_obj(2)+ (bfprp_for(1, ifac)* xnod(3) - bfprp_for(3, ifac)*xnod(1))
                        moment_f_obj(3)= moment_f_obj(3)+ (bfprp_for(2, ifac)* xnod(1) - bfprp_for(1, ifac)*xnod(2))


                        call cs_push_MUI_Coupling("forceX",      &
                                                  cdgfbo0(1, ifac),    &
                                                  cdgfbo0(2, ifac),    &
                                                  cdgfbo0(3, ifac),    &
                                                  forces_f_cell(1))

                        call cs_push_MUI_Coupling("forceY",      &
                                                  cdgfbo0(1, ifac),    &
                                                  cdgfbo0(2, ifac),    &
                                                  cdgfbo0(3, ifac),    &
                                                  forces_f_cell(2))
                        call cs_push_MUI_Coupling("forceZ",      &
                                                  cdgfbo0(1, ifac),    &
                                                  cdgfbo0(2, ifac),    &
                                                  cdgfbo0(3, ifac),    &
                                                  forces_f_cell(3))

                    case (2)
            
                        ! calculate the area of the present boundary face with the unit of [m^2]
                        face_area = sqrt(surfbo(1,ifac)**2+surfbo(2,ifac)**2)
            
                        ! obtain the force components of the present boundary face with the unit of [N]
                        forces_f_cell(1) = bfprp_for(1, ifac)
                        forces_f_cell(2) = bfprp_for(2, ifac)

                        ! obtain the traction components of the present boundary face with the unit of [N/m^2]
                        traction_f_cell(1) = forces_f_cell(1)/face_area
                        traction_f_cell(2) = forces_f_cell(2)/face_area

                        ! obtain the global coordinate of the centers of the present boundary face for the moved mesh
                        xnod(1) = cdgfbo(1, ifac)
                        xnod(2) = cdgfbo(2, ifac)

                        ! calculate the global Moment
                        moment_f_obj(1)= moment_f_obj(1)+ (bfprp_for(3, ifac)* xnod(2) - bfprp_for(2, ifac)*xnod(3))

                        call cs_push_MUI_Coupling("tractionX",      &
                                                  cdgfbo0(1, ifac),    &
                                                  cdgfbo0(2, ifac),    &
                                                  cdgfbo0(3, ifac),    &
                                                  traction_f_cell(1))

                    case (1)
            
                        ! calculate the area of the present boundary face with the unit of [m]
                        face_area = sqrt(surfbo(1,ifac)**2)
            
                        ! obtain the force components of the present boundary face with the unit of [N]
                        forces_f_cell(1) = bfprp_for(1, ifac)

                        ! obtain the traction components of the present boundary face with the unit of [N/m]
                        traction_f_cell(1) = forces_f_cell(1)/face_area

                        ! obtain the global coordinate of the centers of the present boundary face for the moved mesh
                        xnod(1) = cdgfbo(1, ifac)

                        call cs_push_MUI_Coupling("tractionX",      &
                                                  cdgfbo0(1, ifac),    &
                                                  cdgfbo0(2, ifac),    &
                                                  cdgfbo0(3, ifac),    &
                                                  traction_f_cell(1))

                    case default
            
                        write(nfecra, *) "{CS} Error number of ndim: ", ndim, ". Please set it as '3', '2' or '1'."
                
                        return 1
            
                    end select

                enddo

                write(nfecra, *) "total push size: ", iii
                print *, "total push size: ", iii


            endif

        else
            
            write(nfecra, *) "{CS} Error number of ttcabs: ", ttcabs
                
            return 1

        endif

        ! Deallocate the temporary array
        
        deallocate(lstelt)

        return

        end subroutine cs_push_force_MUI_Coupling

end module cs_push_force_MUI_Coupling_Module
!!========================================================================

module cs_mui_cu_fr

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_bool, &
                                       c_double, c_char, c_null_char

!===============================================================================

implicit none

!===============================================================================
! Interfaces
!===============================================================================

interface

    function create_muiCouplingFixedRelaxation_c(nArgs, pointSize, initUndRelxCpl, Fworld) &
        bind(C, name="create_muiCouplingFixedRelaxation")
        use iso_c_binding
        type(c_ptr) :: create_muiCouplingFixedRelaxation_c
        real(c_double), value :: initUndRelxCpl
        type(c_ptr), value :: Fworld
        integer(c_int), value :: pointSize
        integer(c_int), value :: nArgs
    end function create_muiCouplingFixedRelaxation_c

    subroutine delete_muiCouplingFixedRelaxation_c(muiCouplingFixedRelaxation) &
        bind(C, name="delete_muiCouplingFixedRelaxation")
        use iso_c_binding
        type(c_ptr), value :: muiCouplingFixedRelaxation
    end subroutine delete_muiCouplingFixedRelaxation_c

    function muiCouplingFixedRelaxation_undRelxCpl_c(muiCouplingFixedRelaxation)  &
        bind(C, name="muiCouplingFixedRelaxation_undRelxCpl")
        use iso_c_binding
        real(c_double) :: muiCouplingFixedRelaxation_undRelxCpl_c
        ! The const qualification is translated into an intent(in)
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
    end function muiCouplingFixedRelaxation_undRelxCpl_c

    function muiCouplingFixedRelaxation_getXDeltaDisp_c(muiCouplingFixedRelaxation, pointN) &
        bind(C, name="muiCouplingFixedRelaxation_getXDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingFixedRelaxation_getXDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
        integer(c_int), value :: pointN
    end function muiCouplingFixedRelaxation_getXDeltaDisp_c
    
    function muiCouplingFixedRelaxation_getYDeltaDisp_c(muiCouplingFixedRelaxation, pointN) &
        bind(C, name="muiCouplingFixedRelaxation_getYDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingFixedRelaxation_getYDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
        integer(c_int), value :: pointN
    end function muiCouplingFixedRelaxation_getYDeltaDisp_c
    
    function muiCouplingFixedRelaxation_getZDeltaDisp_c(muiCouplingFixedRelaxation, pointN) &
        bind(C, name="muiCouplingFixedRelaxation_getZDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingFixedRelaxation_getZDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
        integer(c_int), value :: pointN
    end function muiCouplingFixedRelaxation_getZDeltaDisp_c
    
    function muiCouplingFixedRelaxation_pointSize_c(muiCouplingFixedRelaxation) &
        bind(C, name="muiCouplingFixedRelaxation_pointSize")
        use iso_c_binding
        integer(c_int) :: muiCouplingFixedRelaxation_pointSize_c
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
    end function muiCouplingFixedRelaxation_pointSize_c

    ! void functions maps to subroutines
    subroutine muiCouplingFixedRelaxation_initialize_c(muiCouplingFixedRelaxation, nArgs, pointSize, initUndRelxCpl, Fworld) &
        bind(C, name="muiCouplingFixedRelaxation_initialize")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
        type(c_ptr), value :: Fworld
        integer(c_int), value :: pointSize
        real(c_double), value :: initUndRelxCpl
        integer(c_int), value :: nArgs
    end subroutine muiCouplingFixedRelaxation_initialize_c
    
    subroutine muiCouplingFixedRelaxation_collectResidual_c(muiCouplingFixedRelaxation, &
                                                            fetchMUIx,                  &
                                                            fetchMUIy,                  &
                                                            fetchMUIz,                  &
                                                            disPreX,                    &
                                                            disPreY,                    &
                                                            disPreZ,                    &
                                                            pointN)                     &
        bind(C, name="muiCouplingFixedRelaxation_collectResidual")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
        real(c_double), intent(in), value :: fetchMUIx
        real(c_double), intent(in), value :: fetchMUIy
        real(c_double), intent(in), value :: fetchMUIz
        real(c_double), intent(in), value :: disPreX
        real(c_double), intent(in), value :: disPreY
        real(c_double), intent(in), value :: disPreZ
        integer(c_int), value :: pointN
    end subroutine muiCouplingFixedRelaxation_collectResidual_c
    
    ! void functions maps to subroutines
    subroutine muiCouplingFixedRelaxation_process_c(muiCouplingFixedRelaxation) &
        bind(C, name="muiCouplingFixedRelaxation_process")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingFixedRelaxation
    end subroutine muiCouplingFixedRelaxation_process_c    
    

end interface


!===============================================================================
        
    private

    public :: muiCouplingFixedRelaxation


!===============================================================================


    type muiCouplingFixedRelaxation
        private
        type(c_ptr) :: ptr ! pointer to the muiCouplingFixedRelaxation class
    contains

        procedure :: delete => delete_muiCouplingFixedRelaxation_polymorph ! Destructor for gfortran

        ! Function member
        procedure :: initialize => muiCouplingFixedRelaxation_initialize
        procedure :: pointSize => muiCouplingFixedRelaxation_pointSize
        procedure :: collectResidual => muiCouplingFixedRelaxation_collectResidual
        procedure :: getXDeltaDisp => muiCouplingFixedRelaxation_getXDeltaDisp
        procedure :: getYDeltaDisp => muiCouplingFixedRelaxation_getYDeltaDisp
        procedure :: getZDeltaDisp => muiCouplingFixedRelaxation_getZDeltaDisp
        procedure :: process => muiCouplingFixedRelaxation_process
    end type muiCouplingFixedRelaxation


!===============================================================================

    ! This function will act as the constructor for muiCouplingFixedRelaxation type
    interface muiCouplingFixedRelaxation

        procedure create_muiCouplingFixedRelaxation

    end interface muiCouplingFixedRelaxation


!===============================================================================

contains ! Implementation of the functions. We just wrap the C function here.


    function create_muiCouplingFixedRelaxation(nArgs, pointSize, initUndRelxCpl, Fworld)
        type(muiCouplingFixedRelaxation) :: create_muiCouplingFixedRelaxation
        real(c_double), value :: initUndRelxCpl
        type(c_ptr), value :: Fworld
        integer(c_int), value :: pointSize
        integer(c_int), value :: nArgs
        create_muiCouplingFixedRelaxation%ptr = &
        create_muiCouplingFixedRelaxation_c(nArgs, pointSize, initUndRelxCpl, Fworld)

    end function create_muiCouplingFixedRelaxation

    subroutine delete_muiCouplingFixedRelaxation(this)
        type(muiCouplingFixedRelaxation) :: this
        call delete_muiCouplingFixedRelaxation_c(this%ptr)
    end subroutine delete_muiCouplingFixedRelaxation

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_muiCouplingFixedRelaxation_polymorph(this)
         
        class(muiCouplingFixedRelaxation) :: this
        call delete_muiCouplingFixedRelaxation_c(this%ptr)
    end subroutine delete_muiCouplingFixedRelaxation_polymorph

    real function muiCouplingFixedRelaxation_undRelxCpl(this)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        muiCouplingFixedRelaxation_undRelxCpl = muiCouplingFixedRelaxation_undRelxCpl_c(this%ptr)
    end function muiCouplingFixedRelaxation_undRelxCpl
    
    real function muiCouplingFixedRelaxation_getXDeltaDisp(this, pointN)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        integer, value :: pointN
        muiCouplingFixedRelaxation_getXDeltaDisp = muiCouplingFixedRelaxation_getXDeltaDisp_c(this%ptr, pointN)
    end function muiCouplingFixedRelaxation_getXDeltaDisp

    real function muiCouplingFixedRelaxation_getYDeltaDisp(this, pointN)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        integer, value :: pointN
        muiCouplingFixedRelaxation_getYDeltaDisp = muiCouplingFixedRelaxation_getYDeltaDisp_c(this%ptr, pointN)
    end function muiCouplingFixedRelaxation_getYDeltaDisp

    real function muiCouplingFixedRelaxation_getZDeltaDisp(this, pointN)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        integer, value :: pointN
        muiCouplingFixedRelaxation_getZDeltaDisp = muiCouplingFixedRelaxation_getZDeltaDisp_c(this%ptr, pointN)
    end function muiCouplingFixedRelaxation_getZDeltaDisp

    integer function muiCouplingFixedRelaxation_pointSize(this)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        muiCouplingFixedRelaxation_pointSize = muiCouplingFixedRelaxation_pointSize_c(this%ptr)
    end function muiCouplingFixedRelaxation_pointSize

    subroutine muiCouplingFixedRelaxation_initialize(this, nArgs, pointSize, initUndRelxCpl, Fworld)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        type(c_ptr), value :: Fworld
        integer(c_int), value :: pointSize
        real(c_double), value :: initUndRelxCpl
        integer(c_int), value :: nArgs
        call muiCouplingFixedRelaxation_initialize_c(this%ptr, nArgs, pointSize, initUndRelxCpl, Fworld)
    end subroutine muiCouplingFixedRelaxation_initialize

    subroutine muiCouplingFixedRelaxation_collectResidual(  this,       &
                                                            fetchMUIx,  &
                                                            fetchMUIy,  &
                                                            fetchMUIz,  &
                                                            disPreX,    &
                                                            disPreY,    &
                                                            disPreZ,    &
                                                            pointN)

        class(muiCouplingFixedRelaxation), intent(in) :: this
        real(c_double), intent(in), value :: fetchMUIx
        real(c_double), intent(in), value :: fetchMUIy
        real(c_double), intent(in), value :: fetchMUIz
        integer(c_int), value :: pointN
        real(c_double), intent(in), value :: disPreX
        real(c_double), intent(in), value :: disPreY
        real(c_double), intent(in), value :: disPreZ

        
        call muiCouplingFixedRelaxation_collectResidual_c(  this%ptr,  &
                                                            fetchMUIx, &
                                                            fetchMUIy, &
                                                            fetchMUIz, &
                                                            disPreX,   &
                                                            disPreY,   &
                                                            disPreZ,   &
                                                            pointN)
    end subroutine muiCouplingFixedRelaxation_collectResidual

    subroutine muiCouplingFixedRelaxation_process(this)
        class(muiCouplingFixedRelaxation), intent(in) :: this
        call muiCouplingFixedRelaxation_process_c(this%ptr)
    end subroutine muiCouplingFixedRelaxation_process

end module cs_mui_cu_fr


!!========================================================================


module cs_mui_cu_ak

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_bool, &
                                       c_double, c_char, c_null_char

!===============================================================================

implicit none

!===============================================================================
! Interfaces
!===============================================================================

interface

    function create_muiCouplingAitken_c(nArgs, pointSize, initUndRelxCpl, Fworld) &
        bind(C, name="create_muiCouplingAitken")
        use iso_c_binding
        type(c_ptr) :: create_muiCouplingAitken_c
        integer(c_int), value :: pointSize
        integer(c_int), value :: nArgs
        real(c_double), value :: initUndRelxCpl
        type(c_ptr), value :: Fworld
    end function create_muiCouplingAitken_c

    subroutine delete_muiCouplingAitken_c(muiCouplingAitken) &
        bind(C, name="delete_muiCouplingAitken")
        use iso_c_binding
        type(c_ptr), value :: muiCouplingAitken
    end subroutine delete_muiCouplingAitken_c

    function muiCouplingAitken_undRelxCpl_c(muiCouplingAitken)  &
        bind(C, name="muiCouplingAitken_undRelxCpl")
        use iso_c_binding
        real(c_double) :: muiCouplingAitken_undRelxCpl_c
        ! The const qualification is translated into an intent(in)
        type(c_ptr), intent(in), value :: muiCouplingAitken
    end function muiCouplingAitken_undRelxCpl_c

    function muiCouplingAitken_getXDeltaDisp_c(muiCouplingAitken, pointv) &
        bind(C, name="muiCouplingAitken_getXDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingAitken_getXDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingAitken
        integer(c_int), value :: pointv
    end function muiCouplingAitken_getXDeltaDisp_c
    
    function muiCouplingAitken_getYDeltaDisp_c(muiCouplingAitken, pointv) &
        bind(C, name="muiCouplingAitken_getYDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingAitken_getYDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingAitken
        integer(c_int), value :: pointv
    end function muiCouplingAitken_getYDeltaDisp_c
    
    function muiCouplingAitken_getZDeltaDisp_c(muiCouplingAitken, pointv) &
        bind(C, name="muiCouplingAitken_getZDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingAitken_getZDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingAitken
        integer(c_int), value :: pointv
    end function muiCouplingAitken_getZDeltaDisp_c
    
    function muiCouplingAitken_pointSize_c(muiCouplingAitken) &
        bind(C, name="muiCouplingAitken_pointSize")
        use iso_c_binding
        integer(c_int) :: muiCouplingAitken_pointSize_c
        type(c_ptr), intent(in), value :: muiCouplingAitken
    end function muiCouplingAitken_pointSize_c

    ! void functions maps to subroutines
    subroutine muiCouplingAitken_initialize_c(muiCouplingAitken, nArgs) &
        bind(C, name="muiCouplingAitken_initialize")
        use iso_c_binding
        integer(c_int), value :: nArgs
        type(c_ptr), intent(in), value :: muiCouplingAitken
    end subroutine muiCouplingAitken_initialize_c
    
    subroutine muiCouplingAitken_collectResidual_c( muiCouplingAitken,  &
                                                    fetchMUIx,          &
                                                    fetchMUIy,          &
                                                    fetchMUIz,          &
                                                    disPreX,            &
                                                    disPreY,            &
                                                    disPreZ,            &
                                                    pointv)             &
        bind(C, name="muiCouplingAitken_collectResidual")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingAitken
        real(c_double), intent(in), value :: fetchMUIx
        real(c_double), intent(in), value :: fetchMUIy
        real(c_double), intent(in), value :: fetchMUIz
        real(c_double), intent(in), value :: disPreX
        real(c_double), intent(in), value :: disPreY
        real(c_double), intent(in), value :: disPreZ
        integer(c_int), value :: pointv
    end subroutine muiCouplingAitken_collectResidual_c
    
    ! void functions maps to subroutines
    subroutine muiCouplingAitken_process_c( muiCouplingAitken,  &
                                            iterN,             &
                                            currentIterN)             &
        bind(C, name="muiCouplingAitken_process")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingAitken
        integer(c_int), value :: iterN
        integer(c_int), value :: currentIterN
    end subroutine muiCouplingAitken_process_c    
    

end interface


!===============================================================================
        
    private

    public :: muiCouplingAitken


!===============================================================================


    type muiCouplingAitken
        private
        type(c_ptr) :: ptr ! pointer to the muiCouplingAitken class
    contains

        procedure :: delete => delete_muiCouplingAitken_polymorph ! Destructor for gfortran

        ! Function member
        procedure :: initialize => muiCouplingAitken_initialize
        procedure :: pointSize => muiCouplingAitken_pointSize
        procedure :: collectResidual => muiCouplingAitken_collectResidual
        procedure :: getXDeltaDisp => muiCouplingAitken_getXDeltaDisp
        procedure :: getYDeltaDisp => muiCouplingAitken_getYDeltaDisp
        procedure :: getZDeltaDisp => muiCouplingAitken_getZDeltaDisp
        procedure :: process => muiCouplingAitken_process
    end type muiCouplingAitken


!===============================================================================

    ! This function will act as the constructor for muiCouplingAitken type
    interface muiCouplingAitken

        procedure create_muiCouplingAitken

    end interface muiCouplingAitken


!===============================================================================

contains ! Implementation of the functions. We just wrap the C function here.


    function create_muiCouplingAitken(nArgs, pointSize, initUndRelxCpl, Fworld)
        type(muiCouplingAitken) :: create_muiCouplingAitken
        integer(c_int), value :: pointSize
        integer(c_int), value :: nArgs
        real(c_double), value :: initUndRelxCpl
        type(c_ptr), value :: Fworld
        create_muiCouplingAitken%ptr = &
        create_muiCouplingAitken_c(nArgs, pointSize, initUndRelxCpl, Fworld)
    end function create_muiCouplingAitken

    subroutine delete_muiCouplingAitken(this)
        type(muiCouplingAitken) :: this
        call delete_muiCouplingAitken_c(this%ptr)
    end subroutine delete_muiCouplingAitken

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_muiCouplingAitken_polymorph(this)
         
        class(muiCouplingAitken) :: this
        call delete_muiCouplingAitken_c(this%ptr)
    end subroutine delete_muiCouplingAitken_polymorph

    real function muiCouplingAitken_undRelxCpl(this)
        class(muiCouplingAitken), intent(in) :: this
        muiCouplingAitken_undRelxCpl = muiCouplingAitken_undRelxCpl_c(this%ptr)
    end function muiCouplingAitken_undRelxCpl
    
    real function muiCouplingAitken_getXDeltaDisp(this, pointv)
        class(muiCouplingAitken), intent(in) :: this
        integer, value :: pointv
        muiCouplingAitken_getXDeltaDisp = muiCouplingAitken_getXDeltaDisp_c(this%ptr, pointv)
    end function muiCouplingAitken_getXDeltaDisp

    real function muiCouplingAitken_getYDeltaDisp(this, pointv)
        class(muiCouplingAitken), intent(in) :: this
        integer, value :: pointv
        muiCouplingAitken_getYDeltaDisp = muiCouplingAitken_getYDeltaDisp_c(this%ptr, pointv)
    end function muiCouplingAitken_getYDeltaDisp

    real function muiCouplingAitken_getZDeltaDisp(this, pointv)
        class(muiCouplingAitken), intent(in) :: this
        integer, value :: pointv
        muiCouplingAitken_getZDeltaDisp = muiCouplingAitken_getZDeltaDisp_c(this%ptr, pointv)
    end function muiCouplingAitken_getZDeltaDisp

    integer function muiCouplingAitken_pointSize(this)
        class(muiCouplingAitken), intent(in) :: this
        muiCouplingAitken_pointSize = muiCouplingAitken_pointSize_c(this%ptr)
    end function muiCouplingAitken_pointSize

    subroutine muiCouplingAitken_initialize(this, nArgs)
        class(muiCouplingAitken), intent(in) :: this
        integer(c_int), value :: nArgs
        call muiCouplingAitken_initialize_c(this%ptr, nArgs)
    end subroutine muiCouplingAitken_initialize

    subroutine muiCouplingAitken_collectResidual(   this,       &
                                                    fetchMUIx,  &
                                                    fetchMUIy,  &
                                                    fetchMUIz,  &
                                                    disPreX,    &
                                                    disPreY,    &
                                                    disPreZ,    &
                                                    pointv)

        class(muiCouplingAitken), intent(in) :: this
        real(c_double), intent(in), value :: fetchMUIx
        real(c_double), intent(in), value :: fetchMUIy
        real(c_double), intent(in), value :: fetchMUIz
        integer(c_int), value :: pointv
        real(c_double), intent(in), value :: disPreX
        real(c_double), intent(in), value :: disPreY
        real(c_double), intent(in), value :: disPreZ

        
        call muiCouplingAitken_collectResidual_c(   this%ptr,  &
                                                    fetchMUIx, &
                                                    fetchMUIy, &
                                                    fetchMUIz, &
                                                    disPreX,   &
                                                    disPreY,   &
                                                    disPreZ,   &
                                                    pointv)
    end subroutine muiCouplingAitken_collectResidual

    subroutine muiCouplingAitken_process(   this,   &
                                            iterN,   &
                                            currentIterN)
        class(muiCouplingAitken), intent(in) :: this
        integer(c_int), value :: iterN
        integer(c_int), value :: currentIterN
        call muiCouplingAitken_process_c(   this%ptr, &
                                            iterN, &
                                            currentIterN)
    end subroutine muiCouplingAitken_process

end module cs_mui_cu_ak


!!========================================================================


module cs_mui_cu_iqnils

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_bool, &
                                       c_double, c_char, c_null_char

!===============================================================================

implicit none

!===============================================================================
! Interfaces
!===============================================================================

interface

    function create_muiCouplingIQNILS_c(nArgs, pointSize, initUndRelxCpl, Fworld) &
        bind(C, name="create_muiCouplingIQNILS")
        use iso_c_binding
        type(c_ptr) :: create_muiCouplingIQNILS_c
        integer(c_int), value :: pointSize
        integer(c_int), value :: nArgs
        real(c_double), value :: initUndRelxCpl
        type(c_ptr), value :: Fworld
    end function create_muiCouplingIQNILS_c

    subroutine delete_muiCouplingIQNILS_c(muiCouplingIQNILS) &
        bind(C, name="delete_muiCouplingIQNILS")
        use iso_c_binding
        type(c_ptr), value :: muiCouplingIQNILS
    end subroutine delete_muiCouplingIQNILS_c

    function muiCouplingIQNILS_undRelxCpl_c(muiCouplingIQNILS)  &
        bind(C, name="muiCouplingIQNILS_undRelxCpl")
        use iso_c_binding
        real(c_double) :: muiCouplingIQNILS_undRelxCpl_c
        ! The const qualification is translated into an intent(in)
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
    end function muiCouplingIQNILS_undRelxCpl_c

    function muiCouplingIQNILS_getXDeltaDisp_c(muiCouplingIQNILS, pointv) &
        bind(C, name="muiCouplingIQNILS_getXDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingIQNILS_getXDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
        integer(c_int), value :: pointv
    end function muiCouplingIQNILS_getXDeltaDisp_c
    
    function muiCouplingIQNILS_getYDeltaDisp_c(muiCouplingIQNILS, pointv) &
        bind(C, name="muiCouplingIQNILS_getYDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingIQNILS_getYDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
        integer(c_int), value :: pointv
    end function muiCouplingIQNILS_getYDeltaDisp_c
    
    function muiCouplingIQNILS_getZDeltaDisp_c(muiCouplingIQNILS, pointv) &
        bind(C, name="muiCouplingIQNILS_getZDeltaDisp")
        use iso_c_binding
        real(c_double) :: muiCouplingIQNILS_getZDeltaDisp_c
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
        integer(c_int), value :: pointv
    end function muiCouplingIQNILS_getZDeltaDisp_c
    
    function muiCouplingIQNILS_pointSize_c(muiCouplingIQNILS) &
        bind(C, name="muiCouplingIQNILS_pointSize")
        use iso_c_binding
        integer(c_int) :: muiCouplingIQNILS_pointSize_c
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
    end function muiCouplingIQNILS_pointSize_c

    ! void functions maps to subroutines
    subroutine muiCouplingIQNILS_initialize_c(muiCouplingIQNILS,  &
                                              nArgs, pointSize, undRelxCplMax, &
                                              aitkenIterationN) &
        bind(C, name="muiCouplingIQNILS_initialize")
        use iso_c_binding
        integer(c_int), value :: pointSize
        integer(c_int), value :: aitkenIterationN
        integer(c_int), value :: nArgs
        real(c_double), value :: undRelxCplMax
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
    end subroutine muiCouplingIQNILS_initialize_c
    
    subroutine muiCouplingIQNILS_collectResidual_c( muiCouplingIQNILS,  &
                                                    fetchMUIx,          &
                                                    fetchMUIy,          &
                                                    fetchMUIz,          &
                                                    disPreX,            &
                                                    disPreY,            &
                                                    disPreZ,            &
                                                    pointv)             &
        bind(C, name="muiCouplingIQNILS_collectResidual")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
        real(c_double), intent(in), value :: fetchMUIx
        real(c_double), intent(in), value :: fetchMUIy
        real(c_double), intent(in), value :: fetchMUIz
        real(c_double), intent(in), value :: disPreX
        real(c_double), intent(in), value :: disPreY
        real(c_double), intent(in), value :: disPreZ
        integer(c_int), value :: pointv
    end subroutine muiCouplingIQNILS_collectResidual_c
    
    ! void functions maps to subroutines
    subroutine muiCouplingIQNILS_process_c( muiCouplingIQNILS,  &
                                            iterN)             &
        bind(C, name="muiCouplingIQNILS_process")
        use iso_c_binding
        type(c_ptr), intent(in), value :: muiCouplingIQNILS
        integer(c_int), value :: iterN
    end subroutine muiCouplingIQNILS_process_c    
    

end interface


!===============================================================================
        
    private

    public :: muiCouplingIQNILS


!===============================================================================


    type muiCouplingIQNILS
        private
        type(c_ptr) :: ptr ! pointer to the muiCouplingIQNILS class
    contains

        procedure :: delete => delete_muiCouplingIQNILS_polymorph ! Destructor for gfortran

        ! Function member
        procedure :: initialize => muiCouplingIQNILS_initialize
        procedure :: pointSize => muiCouplingIQNILS_pointSize
        procedure :: collectResidual => muiCouplingIQNILS_collectResidual
        procedure :: getXDeltaDisp => muiCouplingIQNILS_getXDeltaDisp
        procedure :: getYDeltaDisp => muiCouplingIQNILS_getYDeltaDisp
        procedure :: getZDeltaDisp => muiCouplingIQNILS_getZDeltaDisp
        procedure :: process => muiCouplingIQNILS_process
    end type muiCouplingIQNILS


!===============================================================================

    ! This function will act as the constructor for muiCouplingIQNILS type
    interface muiCouplingIQNILS

        procedure create_muiCouplingIQNILS

    end interface muiCouplingIQNILS


!===============================================================================

contains ! Implementation of the functions. We just wrap the C function here.


    function create_muiCouplingIQNILS(nArgs, pointSize, initUndRelxCpl, Fworld)
        type(muiCouplingIQNILS) :: create_muiCouplingIQNILS
        integer(c_int), value :: pointSize
        integer(c_int), value :: nArgs
        real(c_double), value :: initUndRelxCpl
        type(c_ptr), value :: Fworld
        create_muiCouplingIQNILS%ptr = &
        create_muiCouplingIQNILS_c(nArgs, pointSize, initUndRelxCpl, Fworld)
    end function create_muiCouplingIQNILS

    subroutine delete_muiCouplingIQNILS(this)
        type(muiCouplingIQNILS) :: this
        call delete_muiCouplingIQNILS_c(this%ptr)
    end subroutine delete_muiCouplingIQNILS

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_muiCouplingIQNILS_polymorph(this)
         
        class(muiCouplingIQNILS) :: this
        call delete_muiCouplingIQNILS_c(this%ptr)
    end subroutine delete_muiCouplingIQNILS_polymorph

    real function muiCouplingIQNILS_undRelxCpl(this)
        class(muiCouplingIQNILS), intent(in) :: this
        muiCouplingIQNILS_undRelxCpl = muiCouplingIQNILS_undRelxCpl_c(this%ptr)
    end function muiCouplingIQNILS_undRelxCpl
    
    real function muiCouplingIQNILS_getXDeltaDisp(this, pointv)
        class(muiCouplingIQNILS), intent(in) :: this
        integer, value :: pointv
        muiCouplingIQNILS_getXDeltaDisp = muiCouplingIQNILS_getXDeltaDisp_c(this%ptr, pointv)
    end function muiCouplingIQNILS_getXDeltaDisp

    real function muiCouplingIQNILS_getYDeltaDisp(this, pointv)
        class(muiCouplingIQNILS), intent(in) :: this
        integer, value :: pointv
        muiCouplingIQNILS_getYDeltaDisp = muiCouplingIQNILS_getYDeltaDisp_c(this%ptr, pointv)
    end function muiCouplingIQNILS_getYDeltaDisp

    real function muiCouplingIQNILS_getZDeltaDisp(this, pointv)
        class(muiCouplingIQNILS), intent(in) :: this
        integer, value :: pointv
        muiCouplingIQNILS_getZDeltaDisp = muiCouplingIQNILS_getZDeltaDisp_c(this%ptr, pointv)
    end function muiCouplingIQNILS_getZDeltaDisp

    integer function muiCouplingIQNILS_pointSize(this)
        class(muiCouplingIQNILS), intent(in) :: this
        muiCouplingIQNILS_pointSize = muiCouplingIQNILS_pointSize_c(this%ptr)
    end function muiCouplingIQNILS_pointSize

    subroutine muiCouplingIQNILS_initialize(this, &
                                              nArgs, pointSize, undRelxCplMax, &
                                              aitkenIterationN)
        class(muiCouplingIQNILS), intent(in) :: this
        integer(c_int), value :: nArgs
        integer(c_int), value :: pointSize
        integer(c_int), value :: aitkenIterationN
        real(c_double), value :: undRelxCplMax
        call muiCouplingIQNILS_initialize_c(this%ptr, &
                                              nArgs, pointSize, undRelxCplMax, &
                                              aitkenIterationN)
    end subroutine muiCouplingIQNILS_initialize

    subroutine muiCouplingIQNILS_collectResidual(   this,       &
                                                    fetchMUIx,  &
                                                    fetchMUIy,  &
                                                    fetchMUIz,  &
                                                    disPreX,    &
                                                    disPreY,    &
                                                    disPreZ,    &
                                                    pointv)

        class(muiCouplingIQNILS), intent(in) :: this
        real(c_double), intent(in), value :: fetchMUIx
        real(c_double), intent(in), value :: fetchMUIy
        real(c_double), intent(in), value :: fetchMUIz
        integer(c_int), value :: pointv
        real(c_double), intent(in), value :: disPreX
        real(c_double), intent(in), value :: disPreY
        real(c_double), intent(in), value :: disPreZ

        
        call muiCouplingIQNILS_collectResidual_c(   this%ptr,  &
                                                    fetchMUIx, &
                                                    fetchMUIy, &
                                                    fetchMUIz, &
                                                    disPreX,   &
                                                    disPreY,   &
                                                    disPreZ,   &
                                                    pointv)
    end subroutine muiCouplingIQNILS_collectResidual

    subroutine muiCouplingIQNILS_process(   this,   &
                                            iterN)
        class(muiCouplingIQNILS), intent(in) :: this
        integer(c_int), value :: iterN
        call muiCouplingIQNILS_process_c(   this%ptr, &
                                            iterN)
    end subroutine muiCouplingIQNILS_process

end module cs_mui_cu_iqnils


!!========================================================================
module cs_FSI_Coupling_Module

!===============================================================================
! MUI Coupling files
use cs_mui_coupling
use cs_user_mui_coupling
use cs_mui_cu_fr
use cs_mui_cu_ak
use cs_mui_cu_iqnils
! MUI Coupling files end
!===============================================================================

implicit none

    contains
    subroutine cs_FSI_Coupling_Fetch ( coordX, &
                                        coordY, &
                                        coordZ, &
                                        insubiter, &
                                        dispaleX, &
                                        dispaleY, &
                                        dispaleZ, &
                                        pointN, &
                                        fr, &
                                        ak, &
                                        IQNILS)


        use cs_c_bindings
        use cs_mui_coupling
        use cs_user_mui_coupling
        use cs_mui_cu_fr
        use cs_mui_cu_ak
        use cs_mui_cu_iqnils
        
        implicit none

        type(muiCouplingFixedRelaxation) :: fr
        type(muiCouplingAitken) :: ak
        type(muiCouplingIQNILS) :: IQNILS
        double precision dis_x, dis_y, dis_z
        double precision coordX, coordY, coordZ
        double precision dispaleX, dispaleY, dispaleZ
        integer insubiter
        integer pointN

        !===============================================================================
        if ((muiCouplingMethod .ne. 1) .and. (muiCouplingMethod .ne. 2)) then

          
            dis_x = cs_fetch_disp_MUI_Coupling ("dispX",              &
                                                coordX,       &
                                                coordY,       &
                                                coordZ,       &
                                                sub_iteration_MUI_Coupling, &
                                                insubiter)

            dis_y = cs_fetch_disp_MUI_Coupling ("dispY",              &
                                                coordX,       &
                                                coordY,       &
                                                coordZ,       &
                                                sub_iteration_MUI_Coupling, &
                                                insubiter)

            dis_z = cs_fetch_disp_MUI_Coupling ("dispZ",              &
                                                coordX,       &
                                                coordY,       &
                                                coordZ,       &
                                                sub_iteration_MUI_Coupling, &
                                                insubiter)

        endif
 !**************************************************************************************************
        
            select case (muiCouplingMethod)

                case (1)

                ! Loose Coupling do nothing

                case (2)

                    ! call fr%collectResidual(    dis_x,          &
                                                ! dis_y,          &
                                                ! dis_z,          &
                                                ! dispaleX, &
                                                ! dispaleY, &
                                                ! dispaleZ, &
                                                ! pointN)

                case (3)

                    call ak%collectResidual(    dis_x,          &
                                                dis_y,          &
                                                dis_z,          &
                                                dispaleX, &
                                                dispaleY, &
                                                dispaleZ, &
                                                pointN)

                case (4)

                    call IQNILS%collectResidual(    dis_x,          &
                                                    dis_y,          &
                                                    dis_z,          &
                                                    dispaleX, &
                                                    dispaleY, &
                                                    dispaleZ, &
                                                    pointN)

                case default

                    write(nfecra, *) "{CS} ERROR: MUI Coupling Method: ", muiCouplingMethod, " does not recognized!"
                    stop

        end select
            
 !************************************************************************************************** 

    end subroutine cs_FSI_Coupling_Fetch

    subroutine cs_FSI_Coupling_Process ( insubiter, &
                                         fr, &
                                         ak, &
                                         IQNILS)

        use cs_mui_cu_fr
        use cs_mui_cu_ak
        use cs_mui_cu_iqnils

        implicit none

        type(muiCouplingFixedRelaxation) :: fr
        type(muiCouplingAitken) :: ak
        type(muiCouplingIQNILS) :: IQNILS
        integer insubiter

        !===============================================================================
  

        select case (muiCouplingMethod)

            case (1)

                ! Loose Coupling do nothing

            case (2)

                ! Call bound procedures (member functions)
                !call fr%process()

            case (3)

                ! Call bound procedures (member functions)
                call ak%process(insubiter, insubiter)

            case (4)

                ! Call bound procedures (member functions)
                call IQNILS%process(insubiter)

            case default

                write(nfecra, *) "{CS} ERROR: MUI Coupling Method: ", muiCouplingMethod, " does not recognized!"
            stop

        end select

    end subroutine cs_FSI_Coupling_Process


    real function cs_FSI_Coupling_getXDeltaDisp(  coordX, &
                                                coordY, &
                                                coordZ, &
                                                insubiter, &
                                                dispaleX, &
                                                pointN, &
                                                fr, &
                                                ak, &
                                                IQNILS)
        use cs_c_bindings
        use cs_user_mui_coupling
        use cs_mui_cu_fr
        use cs_mui_cu_ak
        use cs_mui_cu_iqnils

        implicit none

        type(muiCouplingFixedRelaxation) :: fr
        type(muiCouplingAitken) :: ak
        type(muiCouplingIQNILS) :: IQNILS

        double precision dis_xx, dis_xx_temp
        double precision dispaleX
        double precision coordX, coordY, coordZ
        integer insubiter
        integer pointN

        !===============================================================================

            select case (muiCouplingMethod)

                case (1)

                    dis_xx = cs_fetch_disp_MUI_Coupling ("dispX",              &
                                                        coordX,       &
                                                        coordY,       &
                                                        coordZ,       &
                                                        sub_iteration_MUI_Coupling, &
                                                        insubiter)

                case (2)

                    !dis_xx = dispaleX + fr%getXDeltaDisp(pointN)   ! displacement /x

                    dis_xx_temp = cs_fetch_disp_MUI_Coupling ("dispX",              &
                                                        coordX,       &
                                                        coordY,       &
                                                        coordZ,       &
                                                        sub_iteration_MUI_Coupling, &
                                                        insubiter)

                    dis_xx = dispaleX + ((dis_xx_temp - dispaleX)*init_und_relx_coupling)



                case (3)

                    dis_xx = dispaleX + ak%getXDeltaDisp(pointN)   ! displacement /x

                case (4)

                    dis_xx = dispaleX + IQNILS%getXDeltaDisp(pointN)   ! displacement /x

                case default

                    write(nfecra, *) "{CS} ERROR: MUI Coupling Method: ", muiCouplingMethod, " does not recognized!"
                    stop

            end select

        cs_FSI_Coupling_getXDeltaDisp = dis_xx

    end function cs_FSI_Coupling_getXDeltaDisp

    real function cs_FSI_Coupling_getYDeltaDisp( coordX, &
                                                coordY, &
                                                coordZ, &
                                                insubiter, &
                                                dispaleY, &
                                                pointN, &
                                                fr, &
                                                ak, &
                                                IQNILS)

        use cs_c_bindings
        use cs_user_mui_coupling
        use cs_mui_cu_fr
        use cs_mui_cu_ak
        use cs_mui_cu_iqnils

        implicit none

        type(muiCouplingFixedRelaxation) :: fr
        type(muiCouplingAitken) :: ak
        type(muiCouplingIQNILS) :: IQNILS
        double precision dis_yy, dis_yy_temp
        double precision dispaleY
        double precision coordX, coordY, coordZ
        integer insubiter
        integer pointN

        !===============================================================================

            select case (muiCouplingMethod)

                case (1)

                    dis_yy = cs_fetch_disp_MUI_Coupling ("dispY",              &
                                                        coordX,       &
                                                        coordY,       &
                                                        coordZ,       &
                                                        sub_iteration_MUI_Coupling, &
                                                        insubiter)

                case (2)

                    !dis_yy = dispaleY + fr%getYDeltaDisp(pointN)   ! displacement /y

                    dis_yy_temp = cs_fetch_disp_MUI_Coupling ("dispY",              &
                                                        coordX,       &
                                                        coordY,       &
                                                        coordZ,       &
                                                        sub_iteration_MUI_Coupling, &
                                                        insubiter)

                    dis_yy = dispaleY + ((dis_yy_temp - dispaleY)*init_und_relx_coupling)

                case (3)

                    dis_yy = dispaleY + ak%getYDeltaDisp(pointN)   ! displacement /y

                case (4)

                    dis_yy = dispaleY + IQNILS%getYDeltaDisp(pointN)   ! displacement /y

                case default

                    write(nfecra, *) "{CS} ERROR: MUI Coupling Method: ", muiCouplingMethod, " does not recognized!"
                    stop

            end select

        cs_FSI_Coupling_getYDeltaDisp = dis_yy

    end function cs_FSI_Coupling_getYDeltaDisp
    
    real function cs_FSI_Coupling_getZDeltaDisp( coordX, &
                                                coordY, &
                                                coordZ, &
                                                insubiter, &
                                                dispaleZ, &
                                                pointN, &
                                                fr, &
                                                ak, &
                                                IQNILS)

        use cs_c_bindings
        use cs_user_mui_coupling
        use cs_mui_cu_fr
        use cs_mui_cu_ak
        use cs_mui_cu_iqnils

        implicit none

        type(muiCouplingFixedRelaxation) :: fr
        type(muiCouplingAitken) :: ak
        type(muiCouplingIQNILS) :: IQNILS
        double precision dis_zz, dis_zz_temp
        double precision dispaleZ
        double precision coordX, coordY, coordZ
        integer insubiter
        integer pointN

        !===============================================================================

            select case (muiCouplingMethod)

                case (1)

                    dis_zz = cs_fetch_disp_MUI_Coupling ("dispZ",              &
                                                        coordX,       &
                                                        coordY,       &
                                                        coordZ,       &
                                                        sub_iteration_MUI_Coupling, &
                                                        insubiter)

                case (2)

                    !dis_zz = dispaleZ + fr%getZDeltaDisp(pointN)   ! displacement /z

                    dis_zz_temp = cs_fetch_disp_MUI_Coupling ("dispZ",              &
                                                        coordX,       &
                                                        coordY,       &
                                                        coordZ,       &
                                                        sub_iteration_MUI_Coupling, &
                                                        insubiter)

                    dis_zz = dispaleZ + ((dis_zz_temp - dispaleZ)*init_und_relx_coupling)

                case (3)

                    dis_zz = dispaleZ + ak%getZDeltaDisp(pointN)   ! displacement /z

                case (4)

                    dis_zz = dispaleZ + IQNILS%getZDeltaDisp(pointN)   ! displacement /z

                case default

                    write(nfecra, *) "{CS} ERROR: MUI Coupling Method: ", muiCouplingMethod, " does not recognized!"
                    stop

            end select

        cs_FSI_Coupling_getZDeltaDisp = dis_zz

    end function cs_FSI_Coupling_getZDeltaDisp

end module cs_FSI_Coupling_Module
!!========================================================================
!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
