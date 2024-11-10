! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2015, 2017-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2020 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> The scheme describes the discretization to use in the simulation.
!!
!! There are two parts that need to be configured:
!!
!! * The `spatial` discretization
!! * And the `temporal` discretization
!!
!! Optionally a `stabilization` method may be defined for the scheme, see the
!! [[atl_stabilization_module]] for more details on that definition.
!!
!! The kind of spatial discretization is chosen via the `name` setting within
!! the `spatial` table. The following spatial schemes are available:
!!
!! * 'modg' Modal DG discretization in 3D
!! * 'modg_2d' Modal DG disretization in 2D
!! * 'modg_1d' Modal DG discretization in 1D
!!
!! The main configuration option for the spatial discretization is the polynomial
!! degree to use to represent the state in the DG elements.
!! This polynomial degree is set by the option `m` in the spatial table.
!! It may either be a simple scalar value, defining a single polynomial degree to
!! be used for all elements in the domain, or a table that provides polynomial
!! degrees based on the local refinement level of elements.
!! This can be achieved either by providing the polynomial degree for each level
!! individually or by choosing a predefined scheme to choose the polynomial
!! degree for the levels.
!!
!! Individual definitions take the following form:
!!
!!```lua
!!  m = {
!!    { level = 4, m = 4 },
!!    { level = 5, m = 6 }
!!  }
!!```
!!
!! A predefined scheme is offered by `'fixedfact'`, where the polynomial degree
!! on each level is computed by the following formula.
!!
!!```fortran
!!  m(iLevel) = nint( base_order * factor**(maxLevel-iLevel)) - 1
!!```
!!
!! Here the `base_order` and `factor` need to be defined in the configuration, where
!! the `base_order` sets the minimal polynomial degree (+1) that is to be used on the
!! finest level. `factor` defines the factor, by which the scheme order is to be
!! increased by each level. The polynomial degree definition for this case looks as
!! follows:
!!
!!```lua
!!  m = {
!!    predefined = 'fixedfact',
!!    base_order = 4,
!!    factor     = 1.5
!!  }
!!```
!!
!! The default factor for the 'fixedfact' scheme is \(\sqrt{2}\), which allows for
!! approximately the same time step restriction across the levels for hyperbolic
!! equations.
!!
!! Besides the polynomial degree `m` it is also possible to choose the polynomial
!! space to use for multidimensional representation.
!! This `modg_space` is either 'P' or 'Q'.
!! P indicates a multidimensional polynomial, where the sum of the mode indices
!! is at most equal to the configured polynomial degree `m`.
!! Q indicates that each index in the different dimensions itself may at most be
!! `m`. The 'Q' space is the default but requires more computational effort and memory
!! especially for 3D simulations.
!!
!! The explicit time integration is configured by the `temporal` table within `scheme`.
!! Following schemes are available:
!!
!! * 'explicitEuler', only for testing! (unstable)
!! * 'explicitRungeKutta', with `steps=4`
!! * 'imexRungeKutta', with `steps=4`
!! * 'explicitRungeKuttaTaylor', with arbitrary number of steps
!! * 'explicitSSPRungeKutta', with `steps=2`
!!
!! The 'imexRungeKutta' scheme should be used when penalization terms are to be used
!! in the flow simulations. The 'explicitRungeKuttaTaylor' is mainly intended for the
!! solution of linear equation systems.
!!
!! The time step width is controlled by a `control` subtable and the time step can
!! either be chosend adaptively according to the CFL condition, or set as a fixed
!! time step. See also [[atl_global_time_integration_module]].
!!
!! A complete definition of the scheme without the optional `stabilization` table
!! takes the following form:
!!
!!```lua
!!  scheme = {
!!    spatial = {
!!      name = 'modg',
!!      modg_space = 'Q',
!!      m = 11
!!    },
!!    temporal = {
!!      name = 'explicitRungeKutta',
!!      steps = 4,
!!      control = {
!!        name = 'cfl',
!!        cfl = 0.8
!!      }
!!    }
!!  }
!!```
!!
!! For details on the optional `stabilization` see the [[atl_stabilization_module]].
!!
module atl_scheme_module
  use aotus_module,                 only: flu_State, aot_get_val
  use aot_table_module,             only: aot_table_open, aot_table_close

  use env_module,                   only: rk, labelLen
  use tem_aux_module,               only: tem_abort
  use tem_tools_module,             only: upper_to_lower
  use tem_stencil_module,           only: tem_stencilHeader_type,             &
    &                                     tem_stencil_map_toTreelmDef,        &
    &                                     d3q6_cxDir, d2q4_cxDir, d1q2_cxDir, &
    &                                     init
  use tem_logging_module,           only: logUnit


  use ply_modg_basis_module,        only: ply_scalProdDualLeg,     &
    &                                     ply_scalProdDualLegDiff, &
    &                                     ply_modg_basis_type

  use atl_modg_scheme_module,       only: atl_modg_scheme_type, &
    &                                     atl_modg_scheme_init
  use atl_modg_2d_scheme_module,    only: atl_modg_2d_scheme_type,    &
    &                                     atl_modg_2d_scheme_init
  use atl_modg_1d_scheme_module,    only: atl_modg_1d_scheme_type,    &
    &                                     atl_modg_1d_scheme_init
  use atl_stabilization_module,     only: atl_stabilization_type, &
    &                                     atl_ini_stabilization

  implicit none

  private

  integer, parameter :: atl_modg_scheme_prp = 6
  integer, parameter :: atl_modg_2d_scheme_prp = 7
  integer, parameter :: atl_modg_1d_scheme_prp = 8

  !> Datatype to specify the timestepping method.
  type atl_local_timestep_type
    !> The local timestep.
    real(kind=rk) :: dt
  end type


  !> type to define a one dimensional stencil for reconstructions.
  type atl_oneDimStencil_type
    !> the 1D stencil in treelm coordinates.
    integer :: stencil

    !> the number of elements in the stencil, including the cell itself
    !! you reconstruct for.
    integer :: nElems

    !> relative position of the stencil elements to the current cell.
    !! Note, that this vector has length (nElems-1) since the current cell
    !! itself is not stored here.
    integer, allocatable :: elemPos(:)

    !> relative position of the stencil elements in negative direction to the
    !! current cell.
    !!
    !! Note, that this vector has length (nElems-1) since the current cell
    !! itself is not stored here. The entries start with the cell that is most
    !! far from the current cell away.
    integer, allocatable :: ngElemPos(:)

    !> for each element of the mesh we store the lowest and highest left shift
    !! that build correct stencils (i.e. correct means: not including
    !! any boundary element).
    !! The first dimension is the number of elements associated with this
    !! stencil. The second dimension is 2, the first is the lowest possible left
    !! shift index the second is the highest possible left shift index.
    integer,allocatable :: bnd(:,:)
  end type

  !> type specifying all informations about the stencil for the dimension
  !! by dimension reconstruction.
  type atl_dimbydimstencil_type
    !> the stencil in x direction
    type(atl_oneDimStencil_type) :: xStencil
    !> the stencil in y direction
    type(atl_oneDimStencil_type) :: yStencil
    !> the stencil in z direction
    type(atl_oneDimStencil_type) :: zStencil
  end type atl_dimbydimstencil_type


  !> type containing all the informations related to the scheme, e.g.:
  !! time and space discretization, scheme order, etc.
  type atl_scheme_type
    !> integer representing the current discretization scheme.
    integer             :: scheme

    !> the number of degrees of freedom for the selected scheme for a single cell
    !! and a single variable of the equation.
    !! For example we have: P1PM => nDofs=4, P2PM = 10). This number includes only the
    !! degrees of freedom which will be stored. We do not include the number of
    !! reconstructed degrees of freedom here!
    integer             :: nDoFs

    !> the number of reconstructed degrees of freedom for the selected scheme
    !! for a single cell and a single variable of the equation (including
    !! the reconstructed degrees of freedoms).
    integer             :: nDoFsRecons

    !> The number of dofs on the faces.
    integer             :: nFaceDofs

    !> variable to specify the space integration.
    ! ToDO VK: think we can delete this datatype, it was used in the weno scheme
    ! futher it is used flux/atl_hlleFlux_module to specify the number of surface points
    !type(space_quadrature_type) :: space_integration

    !> levelwise information of time discretization
    type(atl_local_timestep_type) :: time

    !! if you want to add another scheme you should add it here and give
    !! a unique code above atl_scheme_type%scheme! Please add a comment, too.
    !! usage is specified by atl_scheme_type%scheme.

    !> Parameters of the modal discontinuous Galerkin scheme if
    !! scheme is set to modg.
    type(atl_modg_scheme_type) :: modg

    !> Parameters of the modal discontinuous Galerkin scheme if
    !! scheme is set to modg 2d.
    type(atl_modg_2d_scheme_type) :: modg_2d

    !> Parameters of the modal discontinuous Galerkin scheme if
    !! scheme is set to modg 1d.
    type(atl_modg_1d_scheme_type) :: modg_1d

    !> Informations about the polynomial basis of a MODG scheme.
    type(ply_modg_basis_type) :: modg_basis

    !> The stabilization(s) for the scheme.
    !! Applied one after each other. Starting with index 1, then 2, ...
    type(atl_stabilization_type), allocatable :: stabilization(:)

    !> Precomputed Scalar Products
    real(kind=rk), allocatable :: dl_prod(:,:)
    real(kind=rk), allocatable :: dl_prodDiff(:,:)

    !> Temp Arrays needed for evaluation of physical fluxes
    real(kind=rk), allocatable :: temp_over(:,:,:)
    real(kind=rk), allocatable :: temp_modal(:,:,:)
    real(kind=rk), allocatable :: temp_nodal(:,:,:)

  end type

  public :: atl_scheme_type,          &
    &       atl_init_scheme,          &
    &       atl_dimbydimstencil_type, &
    &       atl_onedimstencil_type,   &
    &       atl_define_SchemeStencil, &
    &       atl_local_timestep_type,  &
    &       atl_modg_scheme_type,     &
    &       atl_modg_scheme_prp,      &
    &       atl_modg_2d_scheme_prp,   &
    &       atl_modg_1d_scheme_prp,   &
    &       atl_schemeID2ndim


contains


  ! ------------------------------------------------------------------------ !
  !> subroutine to intialize a scheme as specified by a given lua script file.
  subroutine atl_init_scheme(me, conf, minlevel, maxlevel)
    ! -------------------------------------------------------------------- !
    !> The global minimum level of the mesh
    integer, intent(in) :: minLevel
    !> The global maximum level of the mesh
    integer, intent(in) :: maxLevel
    !> the scheme you want to initialize.
    type(atl_scheme_type), intent(out) :: me(minlevel:maxlevel)
    !> flu binding to lua configuration file.
    type(flu_State), intent(in) :: conf
    ! -------------------------------------------------------------------- !
    integer :: scheme_table, spatial_table
    character(len=labelLen) :: scheme_name
    character(len=labelLen) :: sname
    integer :: iError, ilevel
    type(atl_stabilization_type), allocatable :: stabilization(:)
    ! -------------------------------------------------------------------- !

    ! open the scheme table
    call aot_table_open(L=conf, thandle=scheme_table, key='scheme')
    if(scheme_table.eq.0) then
      write(logUnit(1),*) 'ERROR in init_kernel_state: no scheme table in ' // &
        & 'lua configuration file found,stopping...'
      call tem_abort()
    end if

    ! Init the stabilzation (same for all the levels)
    call atl_ini_stabilization(conf = conf, parent_table = scheme_table, &
                              & filter = stabilization )
    do iLevel = minlevel, maxlevel
      allocate(me(iLevel)%stabilization(size(stabilization)))
      me(iLevel)%stabilization(:) = stabilization(:)
    end do

    ! open the spatial subtable
    call aot_table_open(L=conf, parent=scheme_table, &
      &                 thandle=spatial_table, key='spatial')

    ! get the name of the scheme
    call aot_get_val(L = conf, thandle = spatial_table, &
      &              key = 'name', &
      &              val = scheme_name, &
      &              ErrCode = iError)

    sname = upper_to_lower(scheme_name)
    sname = adjustl(sname)

    select case(trim(sname))

    case('modg')
      me(:)%scheme = atl_modg_scheme_prp

      write(logUnit(1),*) 'Init MODG scheme ...'
      do ilevel = minlevel, maxlevel
        call atl_modg_scheme_init( me           = me(ilevel)%modg,      &
          &                        nDofs        = me(ilevel)%nDofs,     &
          &                        nFaceDofs    = me(ilevel)%nFaceDofs, &
          &                        conf         = conf,                 &
          &                        thandle      = spatial_table,        &
          &                        currentLevel = iLevel,               &
          &                        maxLevel     = maxLevel              )

      !precompute and store the scalar product between ansatz and test function
      call compute_scalProd_DualLeg(me(iLevel)%dl_prod, me(iLevel)%dl_prodDiff,&
        &                             me(iLevel)%modg%maxpolyDegree)
      end do

    case('modg_2d')
      me(:)%scheme = atl_modg_2d_scheme_prp

      write(logUnit(1),*) 'Init 2D MODG scheme ...'
      do ilevel = minlevel, maxlevel
        call atl_modg_2d_scheme_init( me           = me(ilevel)%modg_2d,   &
          &                           nDofs        = me(ilevel)%nDofs,     &
          &                           nFaceDofs    = me(ilevel)%nFaceDofs, &
          &                           conf         = conf,                 &
          &                           thandle      = spatial_table,        &
          &                           currentLevel = iLevel,               &
          &                           maxLevel     = maxLevel              )
      ! precompute and store the scalar product between ansatz and test function
      call compute_scalProd_DualLeg(me(iLevel)%dl_prod, me(iLevel)%dl_prodDiff,&
        &                             me(iLevel)%modg_2d%maxpolyDegree)
      end do


    case('modg_1d')
      me(:)%scheme = atl_modg_1d_scheme_prp

      write(logUnit(1),*) 'Init 1D MODG scheme ...'
      do ilevel = minlevel, maxlevel
        call atl_modg_1d_scheme_init( me           = me(ilevel)%modg_1d,   &
          &                           nDofs        = me(ilevel)%nDofs,     &
          &                           nFaceDofs    = me(ilevel)%nFaceDofs, &
          &                           conf         = conf,                 &
          &                           thandle      = spatial_table,        &
          &                           currentLevel = iLevel,               &
          &                           maxLevel     = maxLevel              )
      end do


    case default
      write(logUnit(1),*) 'ERROR in init_kernel_state: unknown scheme name ' &
        &            // trim(scheme_name) // ' !'
      write(logUnit(1),*) 'Supported schemes are: '
      write(logUnit(1),*) '* modg'
      write(logUnit(1),*) '* modg_2d'
      write(logUnit(1),*) 'Stopping....'
      call tem_abort()
    end select

    call aot_table_close(L = conf, thandle = spatial_table)
    call aot_table_close(L = conf, thandle = scheme_table)

    ! Now, we init the timestepping scheme.
    do ilevel = minlevel, maxlevel
      call init_local_time_integration(me(ilevel)%time)
    end do



  end subroutine atl_init_scheme
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> precompute the scalar products of the anstaz and test function
  subroutine compute_scalProd_DualLeg (dl_prod, dl_prodDiff, maxPolyDegree)
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable, intent(out) :: dl_prod(:,:)
    real(kind=rk), allocatable, intent(out) :: dl_prodDiff(:,:)
    integer, intent(in) :: maxPolyDegree
    ! -------------------------------------------------------------------- !
    integer :: iTest, iAns
    ! -------------------------------------------------------------------- !

    allocate(dl_prod(2, maxPolyDegree+1))
    allocate(dl_prodDiff(2,maxPolyDegree+1))

    dl_prod = 0.0_rk
    do iTest=1, maxpolyDegree+1
      iAns = iTest-2
      if (iAns >= 1) then
        dl_prod(1,iTest) = ply_scalProdDualLeg(iAns, iTest)
      end if
      dl_prod(2,iTest) = ply_scalProdDualLeg(iTest, iTest)
    end do

    dl_prodDiff = 0.0_rk
    do iTest=2, maxpolyDegree+1
      iAns = iTest-5
      if (iAns >= 1) then
        dl_prodDiff(1,iTest) = ply_scalProdDualLegDiff(iAns, iTest-1)
      end if
      dl_prodDiff(2,iTest) = ply_scalProdDualLegDiff(iTest-1, iTest-1)
    end do

  end subroutine compute_scalProd_DualLeg
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Subroutine do define a specific stencil for a certain scheme.
  !!
  subroutine atl_define_SchemeStencil(nDims, me)
    ! -------------------------------------------------------------------- !
    !> Number of dimensions to consider in the equation.
    integer, intent(in) :: nDims
    !> the neighbor list you want to init.
    type(tem_stencilHeader_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    select case(nDims)
    case (1)
      call init( me     = me,        &
        &        QQN    = 2,         &
        &        QQ     = 2,         &
        &        useAll = .true.,    &
        &        nDims  = 1,         &
        &        label  = 'd1q2',    &
        &        cxDir  = d1q2_cxDir )
    case (2)
      call init( me     = me,        &
        &        QQN    = 4,         &
        &        QQ     = 4,         &
        &        useAll = .true.,    &
        &        nDims  = 2,         &
        &        label  = 'd2q4',    &
        &        cxDir  = d2q4_cxDir )
    case (3)
      call init( me     = me,        &
        &        QQN    = 6,         &
        &        QQ     = 6,         &
        &        useAll = .true.,    &
        &        nDims  = 3,         &
        &        label  = 'd3q6',    &
        &        cxDir  = d3q6_cxDir )
    end select

    call tem_stencil_map_toTreelmDef(me)

  end subroutine atl_define_SchemeStencil
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> summary: routine to init the timestepping scheme.
  subroutine init_local_time_integration(me)
    ! -------------------------------------------------------------------- !
    !> the scheme you want to initialize.
    type(atl_local_timestep_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    me%dt = 0.0_rk
  end subroutine init_local_time_integration
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  function atl_schemeID2ndim(schemeID) result(ndim)
    ! -------------------------------------------------------------------- !
    integer,intent(in) :: schemeID
    integer :: ndim
    ! -------------------------------------------------------------------- !

  if(schemeID .eq. atl_modg_scheme_prp) then
    ndim = 3
  else if (schemeID .eq. atl_modg_2d_scheme_prp) then
    ndim = 2
  else
     ndim = 1
   end if
  end function atl_schemeID2ndim
  ! ------------------------------------------------------------------------ !

end module atl_scheme_module
