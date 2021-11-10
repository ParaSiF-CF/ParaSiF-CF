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

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------

!===============================================================================

!> \brief User subroutine for selection of specific physics module

!> Define the use of a specific physics amongst the following:
!>   - combustion with gas / coal / heavy fuel oil
!>   - compressible flows
!>   - electric arcs
!>   - atmospheric modelling
!>   - radiative transfer
!>   - cooling towers modelling
!>
!>    Only one specific physics module can be activated at once.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixmlpu        indicates if the XML file from the GUI is used
!>                              (1 : yes, 0 : no)
!______________________________________________________________________________!

subroutine usppmo &
 ( ixmlpu )

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu
use coincl
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer ixmlpu

! --- compf: Compressible flows
! ==========

!        if = -1   module not activated
!        if = 0    module activated

if (ixmlpu.eq.0) then

  ippmod(icompf) = -1

endif




return
end subroutine usppmo


!===============================================================================

!> \brief User subroutine for input of model selection parameters.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      ixmlpu       indicates if the XML file from the GUI is used
!>                              used (1: yes, 0: no
!> \param[in, out] iturb        turbulence model
!> \param[in, out] itherm       thermal model
!> \param[in, out] iale         ALE module
!______________________________________________________________________________!

subroutine usipph &
 ( ixmlpu, iturb , itherm, iale )

!===============================================================================
! Module files
!===============================================================================

use entsor, only: nfecra ! No other module should appear here
use optcal, only: irijco ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale

! Local variables

!===============================================================================

!>    In this subroutine, only the parameters which already appear may
!>    be set, to the exclusion of any other.
!>
!>    If we are not using the Code_Saturne GUI:
!>    All the parameters which appear in this subroutine must be set.
!>
!>    If we are using the Code_Saturne GUI:
!>    parameters protected by a test of the form:
!>
!>      if (ixmlpu.eq.0) then
!>         ...
!>      endif
!>
!>    should already have been defined using the GUI, so only
!>    experts should consider removing the test and adapting them here.

!===============================================================================


!< [usipph]

if (ixmlpu.eq.0) then

! --- Turbulence
!       0: Laminar
!      10: Mixing length
!      20: k-epsilon
!      21: k-epsilon (linear production)
!      30: Rij-epsilon, (standard LRR)
!      31: Rij-epsilon (SSG)
!      32: Rij-epsilon (EBRSM)
!      40: LES (Smagorinsky)
!      41: LES (Dynamic)
!      42: LES (WALE)
!      50: v2f (phi-model)
!      51: v2f (BL-v2/k)
!      60: k-omega SST
!      70: Spalart Allmaras
!  For 10, contact the development team before use

  iturb = 0

! --- Thermal model
!      0: none
!      1: temperature
!      2: enthalpy
!      3: total energy (only for compressible module)
!
!  For temperature, the temperature scale may be set later using itpscl
!  (1 for Kelvin, 2 for Celsius).
!
!  Warning: When using specific physics, this value is
!           set automatically by the physics model.

  itherm = 0

! --- Activation of ALE (Arbitrary Lagrangian Eulerian) method
!    -  0: not activates the ALE module
!    -  1: activates the ALE module

  iale = 1

endif

!< [usipph]

!----
! Formats
!----

return
end subroutine usipph


!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

!===============================================================================

!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.

!===============================================================================
!< [usipsu]

! Calculation options (optcal)
! ============================

! In case of restart, read auxiliary restart file ileaux (= 1) or not (0).

! By default, this file is read, but it may be useful to deactivate
! its use when restarting after a preprocessing stage possibly leading
! to a different number of faces (such as simply joining meshes on
! a different architecture or optimization level or with different options).

! Writing of auxiliary restart files may also be deactivated using: iecaux = 0

ileaux = 1
!ntsuit = -2

! --- Time stepping  (0 : uniform and constant
!                     1 : variable in time, uniform in space
!                     2 : variable in time and space
!                    -1 : steady algorithm)

idtvar = 0


! --- Duration
!       ntmabs = absolute number of the last time step required
!         if we have already run 10 time steps and want to
!         run 10 more, ntmabs must be set to 10 + 10 = 20

ntmabs = 800 !!! Multiples of 80 (for t = 0.05d0)


! --- Reference time step
!     The example given below is probably not adapted to your case.


dtref  = 0.005d0

! --- Maximum time step: dtmax
!     Set a value base on characteristic values of your case.
!      otherwise, the code will use a multiple of dtref by default.
!     Example with
!        Ld: "dynamic" length (for example, the domain length)
!        Ud: characteristic flow velocity
!        Lt: thermal length (for example, the domain height gravity-wise)
!        Delta_rho/rho: relative density difference
!        g: gravity acceleration

!     dtmax = min(Ld/Ud, sqrt(Lt/(g.Delta_rho/rho)))

! --- Handling of hydrostatic pressure
!     iphydr = 0 : ignore hydrostatic pressure (by default)
!              1 : with hydrotatic pressure computation to handle the balance
!                  between the pressure gradient and source terms (gravity and
!                  head losses)
!              2 : with hydrostatic pressure computation to handle the imbalance
!                  between the pressure gradient and gravity source term


! = 1: least squares method based on the ﬁrst neighbour cells (cells which share

imrgra = 2

!indicates whether the source terms in transposed gradient and velocity divergence should be taken into account in the momentum equation.

ivisse = 0

ro0    = 1000.d0
viscl0 = 1.0d0
cp0    = 1219.d0

t0 = 20.d0 + 273.15d0
p0 = 0.0d0

xyzp0(1) = 0.d0
xyzp0(2) = 0.d0
xyzp0(3) = 0.d0

! --- Reference length scale in meters for initialization
!       of epsilon (and specific clipping of turbulence, but
!       this is not the default option)
!       Assign a value of the order of the largest dimension of the
!       physical domain in which the flow may develop.
!       If a negative value is set here, or no value set and the GUI not
!       used, the cubic root of the domain will be used.
!       (useful only for turbulence).

almax = 1.0d0

! integer, save iforbr
! field id of the stresses at boundary (if post-processed)
! set it to 0 if the drag and/or lift coefficient of the object need to be calculated.

iforbr = 0

! ALE (Arbitrary Lagrangian Eulerian) related options
!====================================================

! Number of iterations for fluid initialization. Contrary to ntmabs,
! nalinf is not an absolute iteration number, meaning that in case of
! restart calculation nalinf corresponds to the number of iterations
! for fuid initialization beginning from the first current iteration of
! the calculation restart. In general nalinf = 0 in that case.

nalinf = 0

! Maximum number of iterations in case of implicit Fluid Structure Coupling
! with structural calculations (internal and/or external
! (i.e. using Code_Aster)).
! nalimx = 1, in case of explicit FSI algorithm.

nalimx = 1


! Relative precision of sub-cycling Fluid Structure Coupling algorithm.

epalim = 1.d-5


!iteration (yes=1, no=0) to initialize ALE

italin = 1
!< [usipsu]

!----
! Formats
!----

return
end subroutine usipsu


!===============================================================================

!> \brief User subroutine for the input of additional user parameters for
!>        input/output.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp       number of active specific physics models
!______________________________________________________________________________!

subroutine usipes &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use field
use parall
use period
use ihmpre
use post
use ppppar
use ppthch
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii
integer f_id, ifllog

type(var_cal_opt) :: vcopt

!===============================================================================

!>     This subroutine allows setting parameters
!>     which do not already appear in the other subroutines of this file.
!>
!>     It is possible to add or remove parameters.
!>     The number of physical properties and variables is known here.

!===============================================================================

! ! Frequency of log output
ntlist = 1

!< [usipes_ex_02]
! IOPerformance ! do ii = 1, nvar
  ! IOPerformance ! call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  ! IOPerformance ! vcopt%iwarni = 1
  ! IOPerformance ! call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
! IOPerformance ! enddo

! IOPerformance ! call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
! IOPerformance ! vcopt%iwarni = 2
! IOPerformance ! call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! IOPerformance ! call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
! IOPerformance ! vcopt%iwarni = 2
! IOPerformance ! call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)
!< [usipes_ex_02]


! Logging of variables and properties
! (example for velocity: 1 to activate logging output, 0 to deactivate)

!< [usipes_ex_03]
f_id = ivarfl(iu)
ifllog = 1
call field_set_key_int(f_id, keyvis, ifllog)
f_id = ivarfl(iv)
ifllog = 1
call field_set_key_int(f_id, keyvis, ifllog)
f_id = ivarfl(iw)
ifllog = 1
call field_set_key_int(f_id, keyvis, ifllog)
f_id = ivarfl(ipr)
ifllog = 1
call field_set_key_int(f_id, keyvis, ifllog)
!< [usipes_ex_03]

! --- structures output step

!< [usipes_ex_05]
nthist = 1
frhist = -1.d0
!< [usipes_ex_05]

! Postprocessing of variables and properties
! (example for velocity: 1 to activate postprocessing output, 0 to deactivate)

!< [usipes_ex_07]
f_id = ivarfl(iu)
call field_set_key_int(f_id, keyvis, 1)
!< [usipes_ex_07]

! Probes for variables and properties
! (example for velocity)

!< [usipes_ex_08]
f_id = ivarfl(iu)
call field_set_key_int_bits(f_id, keyvis, POST_MONITOR)
!< [usipes_ex_08]

!< [usipes_ex_09]
f_id = ivarfl(iu)
call field_set_key_int_bits(f_id, keyvis, POST_BOUNDARY_NR)

f_id = ivarfl(ipr)
call field_set_key_int_bits(f_id, keyvis, POST_BOUNDARY_NR)
!< [usipes_ex_09]


return
end subroutine usipes


!===============================================================================


!> \brief Initialize non-standard calculation options for the atmospheric version.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine usati1

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use atincl
use atsoil
use atchem
use atimbr
use siream

!===============================================================================

implicit none

!===============================================================================

!----
! End
!----

return
end subroutine usati1


!===============================================================================
! Purpose:
! -------
!
!> 1. Additional Calculation Options
!>    a. Density Relaxation
!>
!> 2. Physical Constants
!>    a.Dynamic Diffusion Coefficient
!>    b.Constants of the chosen model (EBU, Libby-Williams, ...)
!
!> This routine is called:
!>
!>
!>  - Eddy Break Up pre-mixed flame
!>  - Diffusion flame in the framework of ``3 points'' rapid complete chemistry
!>  - Libby-Williams pre-mixed flame
!>  - Lagrangian module coupled with pulverized coal:
!>    Eulerian combustion of pulverized coal and
!>    Lagrangian transport of coal particles
!>  - Pulverised coal combustion
!>  - Fuel (oil) combustion
!
!===============================================================================

subroutine cs_user_combustion

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use ihmpre
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use radiat

!===============================================================================

implicit none

!----
! End
!----

return
end subroutine cs_user_combustion


!===============================================================================

!> \brief User subroutines for input of calculation parameters,
!>        and to initialize variables used for radiative transfer module

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine cs_f_user_radiative_transfer_param

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use ppincl
use ppppar
use radiat

!===============================================================================

implicit none

return

end subroutine cs_f_user_radiative_transfer_param


!===============================================================================

!> \brief User subroutine.

!> Initialize non standard options for the compressible flow scheme such
!> as the variability of the thermal conductivity and the volume viscosity.
!> Their values can be given in the subroutine \ref uscfx2 .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine uscfx1

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use field

!===============================================================================

implicit none

! Arguments

!===============================================================================

!===============================================================================


!----
! End
!----

return
end subroutine uscfx1


!===============================================================================


!> \brief User subroutine.
!>
!> Set values for the reference volumic viscosity, the reference
!> conductivity and the molar mass for compressible flow.
!>
!> Initialize non standard options for the compressible flow scheme such
!> as the hydrostatic equilibrium.
!>
!> In addition to options set in the user subroutine \ref uscfx1 (or in
!> the GUI): this subroutine allows to set a switch to indicate if the
!> molecular viscosity is constant, its values being given in the user
!> subroutine \ref usipsu .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine uscfx2

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

!===============================================================================


!----
! End
!----

return
end subroutine uscfx2


!===============================================================================

!> \brief Definition of cooling tower model and exchange zones

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine cs_user_cooling_towers

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl

!===============================================================================

implicit none

!===============================================================================

!===============================================================================


!----
! End
!----

return
end subroutine cs_user_cooling_towers


!===============================================================================

!> \brief User routine for definition of computation parameters dealing
!>        with Darcy module

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine user_darcy_ini1

!===============================================================================
! Module files
!===============================================================================

use ihmpre, only: iihmpr
use entsor
use darcy_module

!===============================================================================

implicit none

!===============================================================================


!----
! End
!----

return

end subroutine user_darcy_ini1

!===============================================================================
