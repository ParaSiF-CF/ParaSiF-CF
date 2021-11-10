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

!> \file cs_user_extra_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

!--------
! Formats
!--------

! Local variables

!< [loc_var_dec]
integer          ifac
integer          ii
integer          ilelt  , nlelt

double precision xfor(3), tfor(3), sface
double precision, dimension(:,:), pointer :: bfprp_for
double precision, dimension(:), pointer :: cvar_pr

double precision xnod(3)
integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================
if (iforbr.ge.0) then 
    call field_get_val_v(iforbr, bfprp_for)
endif

call field_get_val_s(ivarfl(ipr), cvar_pr)

! Allocate a temporary array for cells or interior/boundary faces selection
allocate(lstelt(max(ncel,nfac,nfabor)))

!===============================================================================
! compute global efforts on a subset of faces
!===============================================================================

if(ntcabs.eq.1) then
    open(file="Forces.csv",unit=impusr(1))
    ! write headings to the dat file FX, FY, FZ are the x,y,z efforts
    write(impusr(1),"(5(a,1x))") " TIME STEP, TIME, FX, FY, FZ"
    close(unit=impusr(1))       
    open(file="Moments.csv",unit=impusr(2))
    ! write headings to the dat file TX, TY, TZ are the x,y,z efforts
    write(impusr(2),"(5(a,1x))") " TIME STEP, TIME, MX, MY,MZ"
    close(unit=impusr(2))
endif

! open the file
open(file="Forces.csv",unit=impusr(1),position="append")
open(file="Moments.csv",unit=impusr(2),position="append")

if (iforbr.ge.0) then
      
    ! set the force components to zero
    do ii =  1, ndim
        xfor(ii) = 0.d0
        tfor(ii) = 0.d0
    enddo

    ! get the cells corresponding to the cylinder surface
    call getfbr('8', nlelt, lstelt)
    !==========

    ! loop over the cells to integrate the force
    ! over the structure surface
    do ilelt = 1, nlelt
        ifac = lstelt(ilelt)
        
        ! Face center
        xnod(1) = cdgfbo(1, ifac)
        xnod(2) = cdgfbo(2, ifac)
        xnod(3) = cdgfbo(3, ifac)
  
        ! update the force
        do ii = 1, ndim
            xfor(ii) = xfor(ii) + bfprp_for(ii, ifac)
        enddo

    
        !Global Moments
  
        tfor(1)= tfor(1)+ (bfprp_for(3, ifac)* xnod(2) - bfprp_for(2, ifac)*xnod(3))
        tfor(2)= tfor(2)+ (bfprp_for(1, ifac)* xnod(3) - bfprp_for(3, ifac)*xnod(1))
        tfor(3)= tfor(3)+ (bfprp_for(2, ifac)* xnod(1) - bfprp_for(1, ifac)*xnod(2))  
 
    enddo

    ! if the calculation is parallel, add the data from the other processes
    if (irangp.ge.0) then
        call parrsm(ndim,xfor)
        call parrsm(ndim,tfor)
  
        ! print  efforts for current timestep
        write(impusr(1),*)ntcabs,"    ",ttcabs,"    ",xfor(1),"    ",xfor(2),"    ",xfor(3),"    "
        write(impusr(2),*)ntcabs,"    ",ttcabs,"    ",tfor(1),"    ",tfor(2),"    ",tfor(3),"    "
    
        ! close
        close(unit=impusr(1))
        close(unit=impusr(2))
  
    endif

endif

! Deallocate the temporary array
deallocate(lstelt)

!----
! End
!----

return
end subroutine cs_f_user_extra_operations
