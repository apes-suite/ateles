! Copyright (c) 2013, 2015-2016, 2018-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2013-2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !

! **************************************************************************** !
!> Module for the description of the equation system of the Ateles solver
!!
!! The equation to be solved is configured by the `equation` table.
!! To select the equation, its name needs to be provided. All further settings
!! depend on the specific equation.
!!
!!```lua
!!  equation = { name = 'equation to solve' }
!!```
!!
!! The `name` has to be string and the following equations are available:
!!
!! * `navier_stokes` Viscid compressible fluid flow 3D,
!!   see [[atl_eqn_nvrstk_module]]
!! * `navier_stokes_2d` Viscid compressible fluid flow 2D,
!!   see [[atl_eqn_nvrstk_module]]
!! * `euler` Inviscid compressible fluid flow 3D,
!!   see [[atl_eqn_euler_module]]
!! * `euler_2d` Inviscid compressible fluid flow 2D,
!!   see [[atl_eqn_euler_module]]
!! * `euler_1d` Inviscid compressible fluid flow 1D,
!!   see [[atl_eqn_euler_module]]
!! * `loclineuler` Locally linearized Euler equations 3D,
!!   see [[atl_eqn_euler_module]]
!! * `loclineuler_1d` Locally linearized Euler equations 1D,
!!   see [[atl_eqn_euler_module]]
!! * `lineareuler` Linearized Euler equations 3D,
!!   see [[atl_eqn_lineareuler_module]]
!! * `lineareuler_2d` Linearized Euler equations 2D,
!!   see [[atl_eqn_lineareuler_module]]
!! * `acoustic` Linear isothermal acoustic equation 3D,
!!   see [[atl_eqn_acoustic_module]]
!! * `acoustic_2d` Linear isothermal acoustic equation 2D,
!!   see [[atl_eqn_acoustic_module]]
!! * `maxwell` Electrodynamics with Maxwell's equations 3D,
!!   see [[atl_eqn_maxwell_module]]
!! * `maxwell_2d` Electrodynamics with Maxwell's equations 2D,
!!   see [[atl_eqn_maxwell_module]]
!! * `maxwelldivcorrection` Maxwell's equations with divergence correction,
!!   see [[atl_eqn_maxwell_module]]
!! * `advection_1d` Scalar linear 1D advection,
!!   see [[atl_eqn_advection_1d_module]]
!! * `heat` Heat diffusion 3D,
!!   see [[atl_eqn_heat_module]]
!! * `heat_2d` Heat diffusion 2D,
!!   see [[atl_eqn_heat_module]]
!! * `heat_1d` Heat diffusion 1D,
!!   see [[atl_eqn_heat_module]]
!! * `filtered_navier_stokes` Viscid compressible fluid flow with turbulence 3D,
!!   see [[atl_eqn_nvrstk_module]]
!! * `filtered_navier_stokes_2d` Viscid compressible fluid flow with turbulence
!!   2D, see [[atl_eqn_nvrstk_module]]
!! * `nernstplanck` The Nernst-Planck equations,
!!   see [[atl_eqn_nerplanck_module]]
!! * `bbmem` Black-Box Membrane,
!!   see [[atl_eqn_bbm_module]]
!!
!! A note for the lower dimensional systems: the treelm meshes are always of 3D
!! structure, but for lower dimensional equations the higher dimensions are
!! ignored. If there are multiple elements in the superfluous dimensions, the
!! equations will be solved there independently.
!!
!! This module defines the [[atl_Equations_type]] that provides all equation
!! specific data of the Ateles solver.
!!
module atl_equation_module

  ! include treelm modules
  use env_module,                     only: rk, labelLen
  use tem_bc_module,                  only: tem_bc_state_type
  use tem_coordinate_module,          only: coordRotation_type
  use tem_varSys_module,              only: tem_varSys_type
  use tem_spacetime_fun_module,       only: tem_st_fun_linkedList_type
  use tem_stringKeyValuePair_module,  only: grw_stringKeyValuePairArray_type
  use tem_varMap_module,              only: tem_varMap_type

  ! include aotus modules
  use aotus_module,                   only: flu_State
  use aot_path_module,                only: aot_path_type

  ! Ateles modules
  use atl_materialFun_module,         only: atl_materialFun_type

  ! include atl_eqn_* modules
  use atl_eqn_advection_1d_module,    only: atl_advection_1d_type
  use atl_eqn_heat_module,            only: atl_heat_type
  use atl_eqn_bbm_module,             only: atl_bbmem_type
  use atl_eqn_maxwell_module,         only: atl_maxwell_type
  use atl_eqn_euler_module,           only: atl_Euler_type
  use atl_eqn_LinearEuler_module,     only: atl_LinearEuler_type
  use atl_eqn_acoustic_module,        only: atl_acoustic_type
  use atl_eqn_nvrstk_module,          only: atl_NavierStokes_type, &
    &                                       atl_FiltNavierStokes_type
  use atl_eqn_nerplanck_module,       only: atl_NernstPlanck_type


  implicit none
  private

  public :: atl_equations_type
  public :: atl_eqn_var_trafo_type
  public :: atl_temp_flux_arrays_type
  public :: atl_eqn_load_bc


  ! ------------------------------------------------------------------------ !
  !> The number of temporary arrays required to evaluate the physical fluxes
  !! can be set from here.
  !! This is required so that they don't have to be decleared in the openmp
  !! parallel region where physical flux calculation takes place
  type atl_temp_flux_arrays_type

    !> Number of temp arrays of size of oversamp.
    integer :: overSamp = 0

    !> Number of temp arrays of size of modal elements.
    integer :: modal = 0

    !> Number of temp arrays of size of nodal elements.
    integer :: nodal = 0

    !> The size of the variables needed (can be more than equation%nScalars)
    integer :: nScal = 0

  end type atl_temp_flux_arrays_type
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Datatype representing the equation which is used for the simulation.
  !!
  !! The <tt>eq_kind</tt> component specifies which eqation system is simulated
  !! (e.g. the navier-stokes equations).
  !!
  !! For each possible equation system there is a component to store necessary
  !! parameters (e.g. ideal gas constant, viscosity...).
  !! The definition of these subtypes can be found in their own modules
  !! (atl_eqn_*_module).
  !!
  !! Provides information about the variables of the current equation system
  !! (variable names, derived quantities) and possible source terms.
  !!
  type atl_Equations_type
    !> The type of the equation
    !!
    !! possible values:
    !! * 'advection_1d'           :  1D advection advection
    !! * 'bbmem'                  :  blackbox membrane model equations
    !! * 'euler'                  :  (3D) Euler equations
    !! * 'euler_2d'               :  2D Euler equations
    !! * 'navier_stokes'          :  Navier-Stokes equations
    !! * 'filtered_navier_stokes' :  Navier-Stokes with turbulence modeling
    !! * 'maxwell'                :  Maxwell equations
    !! * 'maxwelldivcorrection'   :  Maxwell equations with divergence corection
    !! * 'nernstplanck'           :  Nernst-Planck equation
    !! * 'acoustic'               :  Linearized Gas Dynamics -isentropic
    !! * 'acoustic_2d'            :  2D Linearized Gas Dynamics -isentropic
    !! * 'heat_1d'                :  1D Heat Equation
    !! * 'heat_2d'                :  2D Heat Equation
    !! * 'heat'                   :  3D Heat Equation
    !! * 'linearEuler'            :  3d Linearized euler equation
    !! If you want to add an additional equation you should also add a
    !! new character representing the equation type here. DO NOT FORGET to
    !! add a comment above!
    character(len=labelLen)          :: eq_kind = ''
    !> Number of dimensions of the scheme
    integer                          :: nDimensions
    !> Flag for adaptive timestep calcualtion
    logical                          :: adaptive_timestep
    !> Euler equations parameters
    type(atl_Euler_type)             :: Euler
    !> Navier-Stokes equations parameters (additional to Euler)
    type(atl_NavierStokes_type)      :: NavierStokes
    !> filtered Navier-Stokes equation parameters
    !! (additional to NavierStokes and Euler)
    type(atl_FiltNavierStokes_type)  :: FiltNavierStokes
    !> Pure Maxwell equations parameters
    type(atl_maxwell_type)           :: Maxwell
    !> Advection-Diffusion equation parameters
    type(atl_advection_1d_type)      :: advection
    !> Membrane equations
    type (atl_BBMEM_type)            :: BBMEM
    !> Nernst-Planck equation
    type (atl_NernstPlanck_type)     :: NERPLANCK
    !> Lineraized Gas Dynamics/Acoustic equation
    type (atl_acoustic_type)         :: Acoustic
    !> Heat 1D,2D,3D equation parameters
    type(atl_heat_type)              :: heat
    !> Lineraized Euler equation, 3D and 2D
    type (atl_LinearEuler_type)      :: LinearEuler


    !> Variables of the equation system.
    !!
    !! The conservative variables are stored first.
    !! Additional variables for primitive and characteristic variable systems
    !! can be stored as derived quantities if needed.
    !!
    !!```
    !! +----------------------------+
    !| | consVar | derived vars     |
    !! +----------------------------+
    !!```
    !!
    !! Not all variable kinds have to be there.
    !! The entries of primitive variables can be accessed indirectly
    !! through varSys%Variable(primVar(i)).
    type(tem_varSys_type) :: varSys

    !> Contains all available space-time-functions.
    !!
    !! This list contains all space-time-function, regardless of their purpose.
    !! I.e. source terms, material parameters, etc. can point to this particular
    !! list to reference their space-time-functions.
    type(tem_st_fun_linkedList_type) :: stFunList

    !> The number of derivatives we have to use in our simulation.
    !! Zero means that we use only cell values. One means that we
    !! calculate first derivatives (for each spatial direction).
    integer :: nDerivatives = 0

    !> Does this equation type have primitive variables?
    logical :: hasPrimitiveVariables = .false.

    !> Index of primitive variables in varSys, not allocated if
    !! the current equation has no primitive variables.
    integer, allocatable :: primVar(:)

    !> Index of state variables in varSys
    integer, allocatable :: stateVar(:)

    !> Function pointer to the routine loading boundary conditions according
    !! to the equation system.
    procedure(atl_eqn_load_bc), pointer :: load_bc => NULL()

    !> Function pointer to transform conservative variables to primitive
    !! variables (if the equation has primitive variables)
    procedure(eqn_varElem_trafo), pointer, pass(equation) :: cons2prim => NULL()

    !> Function pointer to transform primitive variables to conservative
    !! variables (if the equation has primitive variables)
    procedure(eqn_varElem_trafo), pointer, pass(equation) :: prim2cons => NULL()

    !> Permutations for all variables of equation system to transform
    !! x,y,z axes to x axes aligned data.
    type(coordRotation_type) :: varRotation(3)

    !> Flag to indicate, whether the equation system requires estimates on
    !! polynomial deviations during computation.
    logical :: requiresDeviation = .false.

    !> Flag to indicate, whether the equation system requires estimates on
    !! gradients of the polynomials during computation.
    logical :: requires_gradmax = .false.

    !> Flag to determine if the given equation system is nonlinear
    logical :: isNonlinear = .false.

    !>block of temporary arrays that can be used within the flux computation for
    ! element local stuff but on multiple threads.
    type(atl_temp_flux_arrays_type) :: temp

    !> Common information about the material the equation system is capable to
    !! use.
    type(atl_materialFun_type) :: material

    !> Maps to reduction_transient operation variables in varSys
    type(tem_varMap_type) :: redTransVarMap

  end type atl_Equations_type
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Description of variable transformations from one system to another and back
  !! again.
  !!
  !! (Systems could be: conservative, primitive, characteristic, ...)
  !! Currently mainly used to transform variables for the common interface
  !! of the boundary conditions, as they might be given in different systems.
  type atl_eqn_var_trafo_type
    !> Flag to indicate if transformation is necessary at all.
    logical :: identity

    !> Transformation from a given system to the required one.
    procedure(eqn_var_trafo), nopass, pointer :: from => NULL()

    !> Inverse, back transformation again.
    procedure(eqn_var_trafo), nopass, pointer :: to => NULL()
  end type atl_eqn_var_trafo_type
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> interface descriptoin for a routine that loads boundary conditions
  !! (currently not really used yet)
  abstract interface
    subroutine atl_eqn_load_bc(equation, bc_state, bc_state_gradient, &
      &                        bc_varDict, bc_varDict_gradient,       &
      &                        bc_normal_vec, bc_normal_vec_gradient, &
      &                        bc_trafo, bc_trafo_gradient,           &
      &                        bc_label, bc_kind, thandle, conf       )

      import :: atl_equations_type, tem_bc_state_type, atl_eqn_var_trafo_type, &
        &       aot_path_type, flu_State, grw_stringKeyValuePairArray_type

      !> Contains everything
      class(atl_equations_type), intent(inout) :: equation
      !> boundary state variable definitions loaded from config file
      type(tem_bc_state_type), allocatable, intent(out) :: bc_state(:)
      !> boundary state gradient variable definitions loaded from config file
      type(tem_bc_state_type), allocatable, intent(out) :: bc_state_gradient(:)
      !> Dictionary of boundary variables in bc_state
      type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict
      !> Dictionary of boundary variables in bc_state_gradient
      type(grw_stringKeyValuePairArray_type), intent(out) :: bc_varDict_gradient
      logical, intent(out) :: bc_normal_vec
      logical, intent(out) :: bc_normal_vec_gradient
      type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo
      type(atl_eqn_var_trafo_type), intent(out) :: bc_trafo_gradient
      character(len=*), intent(in) :: bc_label
      character(len=*), intent(in) :: bc_kind
      integer, intent(in) :: thandle
      type(flu_State) :: conf
    end subroutine atl_eqn_load_bc
  end interface
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> interface description for the transormation of a set of state varibles to
  !! another set of state variables (e.g. conservative variables to primitive
  !! variable or vice versa)
  abstract interface
    subroutine eqn_var_trafo( equation, instate, outstate, material)
      import :: atl_equations_type, rk
      !> the eqations type that defines all necessary parameters
      class(atl_equations_type), intent(in) :: equation

      !> input state vector
      !!
      !! The array dimensions reflect the way the state vector is stored in
      !! the solver:
      !! dimension (nPnts nVars) with:
      !! * nPnts (used to store state at several points in each element)
      !! * nVars: number of necessary variables to define the state
      real(kind=rk), intent(inout)      :: instate(:,:)

      !> output transformed state vector
      !!
      !! The array dimensions reflect the way the state vector is stored in
      !! the solver:
      !! dimension (nPnts, nVars) with:
      !! * nPnts same as in instate
      !! * nVars: number of necessary variables to define the state
      real(kind=rk), optional, intent(out)    :: outstate(:,:)

      !> The material information for the state transformation.
      !!
      !! The array dimensionals are: (nDofs, nMaterials)
      !! * nDoFs same as in instate
      !! * nMaterials: number of material parameters.
      real(kind=rk), optional,  intent(in) :: material(:,:)

    end subroutine eqn_var_trafo
  end interface
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> interface description for the transormation of a set of state varibles to
  !! another set of state variables (e.g. conservative variables to primitive
  !! variable or vice versa)
  abstract interface
    subroutine eqn_varElem_trafo(equation, instate, outstate, nElems)
      import :: atl_equations_type, rk
      !> the eqations type that defines all necessary parameters
      class(atl_equations_type), intent(in) :: equation

      !> input state vector
      !!
      !! The array dimensions reflect the way the state vector is stored in
      !! the solver:
      !! dimension (nDoFs nVars) with:
      !!
      !! * nMaxInElems >= nElems,
      !! * nDoFs arbitrary (used to store state at several points in each
      !!   element)
      !! * nVars: number of necessary variables to define the state
      real(kind=rk), intent(inout)      :: instate(:,:,:)

      !> output transformed state vector
      !!
      !! The array dimensions reflect the way the state vector is stored in
      !! the solver:
      !! dimension (nElems, nDoFs, nVars) with:
      !!
      !! * nDoFs same as in instate
      !! * nVars: number of necessary variables to define the state
      real(kind=rk), optional, intent(out)    :: outstate(:,:,:)

      !> The number of elements.
      integer, intent(in) :: nElems

    end subroutine eqn_varElem_trafo
  end interface
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> interface description for the calculation of source terms; adds  necessary values
  !! to the right hand side of the equation system from source parameters (from a
  !! space-time function) and the state (from which you can possibly derive other quantities)
  abstract interface
    subroutine eqn_evaluate_source(equation, state, sourceParameter, RHS, &
                                  & background_material, nElems)
      import :: atl_equations_type, rk
      !> the eqations type that defines all necessary parameters
      class(atl_equations_type), intent(in) :: equation

      !> input state vector
      !!
      !! The array dimensions reflect the way the state vector is stored in
      !! the solver:
      !! dimension (1:nMaxInElems, nDoFs, nVars) with:
      !!
      !! * nMaxInElems >= nElems,
      !! * nDoFs arbitrary (used to store state at several points in each
      !!   element)
      !! * nVars: number of necessary variables to define the state
      real(kind=rk), intent(in)         :: state(:,:,:)

      !> vector with the parameters of the source
      !!
      !! The array dimensions reflect the way the state vector is stored in
      !! the solver:
      !! dimension (1:nElems, nDoFs, nSourceParameter) with:
      !!
      !! * nElems, number elements to calculate RHS for (< nMaxInElems!)
      !! * nDoFs arbitrary (used to store state at several points in each
      !!   element)
      !! * nSourceParameter: same as nComponents of the tem_source_type
      real(kind=rk), intent(in)         :: sourceParameter(:,:,:)

      !> right hand side vector
      !!
      !! attention: always updat right hand side with RHS = RHS + ...
      !!
      !! same dimensions as state (1:nElems, nDoFs, nVars) with:
      !! * nElems,
      !! * nDoFs arbitrary (used to store state at several points in each
      !!   element)
      !! * nVars: number of necessary variables to define the state
      real(kind=rk), intent(inout)       :: RHS(:,:,:)

      !> Description of the background material on the current level.
      real(kind=rk), intent(in) :: background_material(:)

      !> number of elems to calculate RHS for,
      !! if elemInd is not given, calculate RHS for (1:nElems,:,:)
      integer, intent(in)                :: nElems

    end subroutine eqn_evaluate_source
  end interface


end module atl_equation_module
! ******************************************************************************
