! Copyright (c) 2017-2019, 2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017-2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Peter Vitt and Daniel
! Fleischer for University of Siegen.
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
!> Adaptive sampling of polynomial data.
!!
!! This module implements the sampling of polynomials with data dependent
!! refinement.
!! Elements, where the polynomials vary above a certain threshold, will be
!! split into their eight children and the polynomial data will be projected
!! onto those.
!! The polynomials in the children can be restricted in their order to limit
!! the memory consumption.
!! In the end only one degree of freedom will be returned for each (refined)
!! element, these are always the mean of the solution in those (refined)
!! elements.
!!
!! The module provides a data type to describe the configuration of the
!! adaptive sampling: [[ply_sampling_adaptive_type]],
!! one routine to load this configuration [[ply_sampling_adaptive_load]] and
!! one routine to actually perform the adaptive sampling
!! [[ply_sample_adaptive]].
module ply_sampling_adaptive_module
  use mpi

  use iso_c_binding, only: c_f_pointer, c_loc
  use env_module, only: rk, labelLen

  use aotus_module, only: flu_state,   &
    &                     aot_get_val

  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort
  use tem_bc_prop_module, only: tem_bc_prop_type, &
    &                           tem_bc_prop_pos,  &
    &                           tem_bc_prop_sublist
  use tem_logging_module, only: logunit
  use tem_refining_module, only: tem_refine_global_subtree
  use tem_subtree_type_module, only: tem_subtree_type, &
    &                                tem_destroy_subtree
  use tem_subTree_module, only: tem_subTree_from,      &
    &                           tem_create_subTree_of, &
    &                           tem_create_tree_from_sub
  use tem_time_module, only: tem_time_type
  use tem_topology_module, only: tem_levelOf
  use tem_tracking_module, only: tem_tracking_instance_type, &
    &                            tem_tracking_config_type
  use tem_tools_module, only: upper_to_lower
  use tem_varsys_module, only: tem_varSys_proc_element,       &
    &                          tem_varSys_proc_point,         &
    &                          tem_varSys_proc_getParams,     &
    &                          tem_varSys_proc_setParams,     &
    &                          tem_varSys_proc_setupIndices,  &
    &                          tem_varSys_proc_getValOfIndex, &
    &                          tem_varsys_append_statevar,    &
    &                          tem_varSys_type,               &
    &                          tem_varSys_op_type,            &
    &                          tem_varSys_getParams_dummy,    &
    &                          tem_varSys_setParams_dummy

  use ply_sampling_varsys_module, only: ply_sampling_var_type,         &
    &                                   ply_sampling_varsys_for_track, &
    &                                   ply_sampling_var_allocate,     &
    &                                   ply_sampling_var_move,         &
    &                                   ply_sampling_var_compute_elemdev
  use ply_split_element_module, only: ply_split_element_init, &
    &                                 ply_split_element,      &
    &                                 ply_split_element_1D,   &
    &                                 ply_split_element_2D,   &
    &                                 ply_split_element_3D
  use ply_filter_element_module, only: ply_filter_element,      &
    &                                  ply_filter_element_type, &
    &                                  ply_filter_element_load

  implicit none

  private

  public :: ply_sample_adaptive
  public :: ply_sampling_adaptive_type
  public :: ply_sampling_adaptive_load


  !> Constant to indicate the factor reduction mode.
  integer, parameter :: redux_factor = 1

  !> Constant to indicate the decrement reduction mode.
  integer, parameter :: redux_decrement = 2

  !> Configuration of the adaptive sampling.
  !!
  !! The main setting is max_nlevels, which states the maximum number of
  !! levels that elements will be refined.
  type ply_sampling_adaptive_type
    !> Maximal number of levels by which any mesh element should be refined.
    !!
    !! A setting of 0 results in no sampling, and the original mesh elements
    !! will be used with the integral mean value (first degree of freedom).
    !! Higher levels provide a limit for the refinement of the mesh.
    !! Note, that even for large settings here, the overall mesh depth is
    !! restricted by the global limit, due to the available space in the
    !! integers representing the treeIDs.
    integer :: max_nlevels = 0

    !> Maximum allowed oscillation of the solution.
    !! For adaptive subsampling only.
    real(kind=rk) :: eps_osci

    !> Method to use for the reduction.
    !!
    !! This may either be:
    !!
    !! - redux_factor: multiply the maximal polynomial degree by the
    !!               given factor in each refinement level. This allows to
    !!               maintain the same total number of degrees of freedom by
    !!               halfing the modes during each refinement.
    !!               This is the default reduction mode
    !! - redux_decrement: Cut off the given dof_decrement last modes in each
    !!                    refinement. This can be used to filter off the most
    !!                    oscillatory modes while affecting the solution
    !!                    minimally otherwise.
    integer :: reduction_mode

    !> Indication whether to filter modes during refinement by ignoring
    !! all modes in the parent, that exceed the target polynomial degree
    !! of the child elements.
    !!
    !! This provides a simple lowpass filtering method if activated.
    !! Defaults to false.
    logical :: ignore_highmodes = .false.

    !> Number of modes to cut off in each refinement.
    !!
    !! If the decrement mode for reduction is used, this setting will
    !! be used to cut off as many modes from the refined elements.
    integer :: dof_decrement = 1

    !> Factor to Reduce dofs for every sampling level.
    !! Can be used to avoid too drastic increase of memory consumption.
    !! For adaptive subsampling only.
    real(kind=rk) :: dofReducFactor

    !> Indicator for the limitation of memory consumption.
    logical :: adaptiveDofReduction

    !> Absolute upper bound level to refine to.
    integer :: AbsUpperBoundLevel

    !> Filtering the poylnomial modes during adaptive refinement.
    !!
    !! This filtering provides the possibility to change the applied
    !! filtering based on the polynomials and thereby attempting to
    !! capture discontinuities more sharply.
    type(ply_filter_element_type) :: filter_element
  end type ply_sampling_adaptive_type


  !> Small helping type to allow arrays of arrays for the variable data.
  type realarray_type
    real(kind=rk), pointer :: dat(:) => NULL()
  end type realarray_type

  !> A container for the method data to hold the data in a scalar pointer for
  !! the C-pointer conversion.
  type sampled_method_data_type
    type(realarray_type), allocatable :: component(:)
  end type sampled_method_data_type


contains


  ! ------------------------------------------------------------------------- !
  !> Load the configuration for adaptive subsampling.
  subroutine ply_sampling_adaptive_load(me, conf, parent)
    ! -------------------------------------------------------------------- !
    !> Sampling definition to load.
    type(ply_sampling_adaptive_type), intent(out) :: me

    !> Configuration to read the sampling settings from.
    type(flu_State), intent(in) :: conf

    !> Parent table in which to look for the adaptive sampling settings.
    integer, intent(in), optional :: parent
    ! -------------------------------------------------------------------- !
    integer :: iError
    character(len=labelLen) :: reduction_mode
    ! -------------------------------------------------------------------- !

    call aot_get_val( L       = conf,           &
      &               thandle = parent,         &
      &               key     = 'nlevels',      &
      &               val     = me%max_nlevels, &
      &               ErrCode = iError,         &
      &               default = 0               )

    call aot_get_val( L       = conf,        &
      &               thandle = parent,      &
      &               key     = 'tolerance', &
      &               val     = me%eps_osci, &
      &               ErrCode = iError,      &
      &               default = 0.0_rk       )
    write(logunit(1),*) '  Using a tolerance of ', me%eps_osci

    call aot_get_val( L       = conf,                  &
      &               thandle = parent,                &
      &               key     = 'AbsUpperBoundLevel',  &
      &               val     = me%AbsUpperBoundLevel, &
      &               ErrCode = iError,                &
      &               default = 0                      )

    if (me%AbsUpperBoundLevel > 0) then
      write(logunit(1),*) '  Level ', me%AbsUpperBoundLevel, &
        &                 ' will be used as absolute upper bound for the'
      write(logunit(1),*) '  refinement. No element will be refined beyond this'
      write(logunit(1),*) '  level.'
    end if

    call aot_get_val( L       = conf,                &
      &               thandle = parent,              &
      &               key     = 'ignore_highmodes',  &
      &               val     = me%ignore_highmodes, &
      &               ErrCode = iError,              &
      &               default = .false.              )

    if (me%ignore_highmodes) then
      write(logunit(1),*) 'The highest modes that exceed the target degree'
      write(logunit(1),*) 'will be ignored during each refinement!'
    end if

    call aot_get_val( L       = conf,              &
      &               thandle = parent,            &
      &               key     = 'reduction_mode',  &
      &               val     = reduction_mode,    &
      &               ErrCode = iError,            &
      &               default = 'factor'           )

    reduction_mode = upper_to_lower(reduction_mode)
    select case(trim(reduction_mode))
    case ('factor')
      me%reduction_mode = redux_factor
      call aot_get_val( L       = conf,              &
        &               thandle = parent,            &
        &               key     = 'dof_reduction',   &
        &               val     = me%dofReducFactor, &
        &               ErrCode = iError,            &
        &               default = 0.5_rk             )

      if (me%dofReducFactor > 1.0_rk .or. me%dofReducFactor <= 0.0_rk) then
        write(logunit(1),*) 'ERROR: dof_reduction has invalid setting:', &
          &                 me%dofReducFactor
        call tem_abort( 'dof_reduction needs to be in ' &
          & // '0.0 < dof_reduction <= 1.0'             )
      end if

      call aot_get_val( L       = conf,                    &
        &               thandle = parent,                  &
        &               key     = 'adaptiveDofReduction',  &
        &               val     = me%adaptiveDofReduction, &
        &               ErrCode = iError,                  &
        &               default = .FALSE.                  )

      if (me%adaptiveDofReduction) then

        write(logunit(1),*) '  Reducing the degrees of freedom adaptively.'
        write(logunit(1),*) '  This option tries to keep as many degrees of' &
          &                 // ' freedom'
        write(logunit(1),*) '  as possible while not increasing the required'
        write(logunit(1),*) '  memory.'
        write(logunit(1),*) '  However, the factor for the reduction will be'
        write(logunit(1),*) '  at least ', me%dofReducFactor
        write(logunit(1),*) '  Even if this results in an increased memory'
        write(logunit(1),*) '  consumption after the refinement.'

      else

        write(logunit(1),*) 'Reducing the degrees of freedom on each ' &
          &                 // 'refinement by a factor of', &
          &                 me%dofReducFactor

      end if

    case ('decrement')
      me%reduction_mode = redux_decrement

      call aot_get_val( L       = conf,             &
        &               thandle = parent,           &
        &               key     = 'dof_decrement',  &
        &               val     = me%dof_decrement, &
        &               ErrCode = iError,           &
        &               default = 1                 )

      write(logunit(1),*) '  Degrees of freedom will be decremented', &
        &                 ' on each level by ', me%dof_decrement

    case default
      write(logunit(1),*) 'ERROR: Unknown reduction mode: ', &
        &                 trim(reduction_mode)
      write(logunit(1),*) 'Available reduction modes are:'
      write(logunit(1),*) '* factor'
      write(logunit(1),*) '* decrement'
      call tem_abort( 'Unknown reduction mode!' )
    end select

    call ply_filter_element_load(      &
      &    me     = me%filter_element, &
      &    conf   = conf,              &
      &    parent = parent             )

  end subroutine ply_sampling_adaptive_load
  ! ------------------------------------------------------------------------- !
  ! ------------------------------------------------------------------------- !


  ! ------------------------------------------------------------------------- !
  !> Sample data described by varsys in orig_mesh according to the tracking
  !! object trackInst with adaptive refinements.
  !!
  !! Only works for Q-Polynomials.
  subroutine ply_sample_adaptive( me, ndims, orig_mesh, orig_bcs, varsys,   &
    &                             var_degree, lvl_degree, trackInst,        &
    &                             trackConfig, time, new_mesh, resvars      )
    ! -------------------------------------------------------------------- !
    !> A ply_sampling_type to describe the sampling method.
    type(ply_sampling_adaptive_type), intent(in) :: me

    !> The original mesh to be refined.
    type(treelmesh_type), intent(in) :: orig_mesh

    !> Boundary conditions for the original mesh.
    type(tem_BC_prop_type), intent(in) :: orig_bcs

    !> Variable system of the original data to do the sampling on.
    type(tem_varsys_type), intent(in) :: varsys

    !> Maximal polynomial degree for each variable.
    !!
    !! Needs to be matching the variable definition in the variable system.
    !! @todo Needs to be changed to be an information per element per variable!
    !!       Possibly by defining a variable in the varsys, providing the
    !!       degree.
    integer, intent(in) :: var_degree(:)

    !> Maximal polynomial degree for each level.
    integer, intent(in) :: lvl_degree(:)

    !> Number of dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    !> Tracking object describing what to sample.
    type(tem_tracking_instance_type), intent(in) :: trackInst

    !> Tracking configuration with the geometry to obtain from the overall mesh.
    type(tem_tracking_config_type), intent(in) :: trackConfig

    !> Point in time to get the data for.
    type(tem_time_type), intent(in) :: time

    !> The new mesh with the refined elements.
    type(treelmesh_type), intent(out) :: new_mesh

    !> Resulting system of variables describing the data in the arrays of
    !! subsampled elements.
    type(tem_varsys_type), intent(out) :: resvars
    ! -------------------------------------------------------------------- !
    type(ply_sampling_var_type), pointer :: var(:) => NULL()
    type(ply_sampling_var_type), pointer :: prev(:) => NULL()

    type(sampled_method_data_type), pointer :: vardat => NULL()

    real(kind=rk), allocatable :: reduction_factor(:)
    real(kind=rk), allocatable :: maxmean(:)
    real(kind=rk), allocatable :: minmean(:)
    real(kind=rk) :: memprefac

    real(kind=rk), pointer :: parent_data(:,:) => NULL()
    real(kind=rk), pointer :: child_data(:,:) => NULL()

    integer :: maxtarget
    integer :: targetdeg
    integer :: containersize

    integer :: nMaxModes
    integer :: nChildren
    integer :: nVars
    integer :: nScalars
    integer :: nComponents
    integer :: nDofs
    integer :: nOldDofs
    integer :: nElems
    integer :: newElems
    integer :: refinedElems

    integer :: varpos
    integer :: bcpos
    integer :: firstdof, lastdof
    integer :: oldfirst, oldlast

    integer :: iVar
    integer :: iScalar
    integer :: iComponent
    integer :: iRefLevel
    integer :: iElem
    integer :: iRefElem
    integer :: iNewElem
    integer :: iChild
    integer :: elemlevel

    integer :: iError

    integer, allocatable :: ReducableElems(:)
    integer, allocatable :: map2global(:)
    integer, allocatable :: maxdeg(:)
    integer, allocatable :: origsize(:)

    logical, allocatable :: reached_limit(:)
    logical, allocatable :: is_varying(:)
    logical :: lastRefine
    logical :: need2refine

    type(treelmesh_type), pointer :: curmesh => NULL()
    type(treelmesh_type), pointer :: oldmesh => NULL()
    type(tem_bc_prop_type), pointer :: curbcs => NULL()
    type(tem_bc_prop_type), pointer :: oldbcs => NULL()

    type(tem_subTree_type) :: refine_subtree
    type(tem_subTree_type) :: tracked_subtree

    procedure(ply_split_element), pointer :: split_element
    procedure(ply_filter_element), pointer :: filtering
    procedure(tem_varSys_proc_element), pointer :: get_element
    procedure(tem_varSys_proc_point), pointer :: get_point
    procedure(tem_varSys_proc_setParams), pointer :: set_params
    procedure(tem_varSys_proc_getParams), pointer :: get_params
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex
    ! -------------------------------------------------------------------- !

    nMaxModes = maxval(var_degree+1)
    nChildren = 2**nDims
    nullify(split_element)
    nullify(filtering)

    select case(nDims)
    case (1)
      split_element => ply_split_element_1D
      filtering     => me%filter_element%filter1D
    case (2)
      split_element => ply_split_element_2D
      filtering     => me%filter_element%filter2D
    case (3)
      split_element => ply_split_element_3D
      filtering     => me%filter_element%filter3D
    end select

    call ply_split_element_init(nMaxModes)

    call ply_sampling_varsys_for_track( varsys        = varsys,     &
      &                                 trackInst     = trackInst,  &
      &                                 mesh          = orig_mesh,  &
      &                                 nDims         = nDims,      &
      &                                 lvl_degree    = lvl_degree, &
      &                                 sample_varsys = resvars,    &
      &                                 var           = var,        &
      &                                 time          = time        )

    ! Create a mesh describing the original selection of elements to sample.
    ! We use two pointers to hold the current and the previous mesh information:
    ! curmesh, curbcs and oldmesh, oldbcs respectively

    allocate(curmesh)
    allocate(curbcs)

    call tem_create_tree_from_sub( intree     = orig_mesh,         &
      &                            subtree    = trackInst%subtree, &
      &                            newtree    = curmesh,           &
      &                            keep_props = .true.             )

    bcpos = tem_bc_prop_pos(orig_mesh)

    nScalars = size(var)
    nElems = curmesh%nElems
    allocate(maxdeg(nScalars))
    allocate(reduction_factor(nScalars))
    allocate(origsize(nScalars))
    allocate(reducableelems(nScalars))
    allocate(maxmean(nScalars))
    allocate(minmean(nScalars))

    do iScalar=1,nScalars
      origsize(iScalar) = size(var(iScalar)%dat)
      maxmean(iScalar) = maxval(var(iScalar)%dat(var(iScalar)%first(:nElems)))
    end do

    ! If the maximal mean of the variable across all elements is too close to
    ! 0, lift it to ensure some safeguard.
    maxmean = max( maxmean, 256*epsilon(maxmean(1)) )
    minmean = maxmean * me%eps_osci


    if (bcpos > 0) then
      ! Take care of the boundary properties
      if (trackInst%subtree%useGlobalMesh) then

        ! For global meshes we can just reuse the original boundary property.
        curbcs = orig_bcs
        curbcs%header => curmesh%global%property(bcpos)
        curbcs%property => curmesh%property(bcpos)

      else

        ! Otherwise we need to create a new boundary property with just the
        ! elements of the subtree.
        call tem_bc_prop_sublist( tree     = orig_mesh,                        &
          &                       bc       = orig_bcs,                         &
          &                       header   = orig_mesh%global%property(bcpos), &
          &                       property = orig_mesh%property(bcpos),        &
          &                       sublist  = trackInst%subtree%map2global,     &
          &                       sub_bc   = curbcs                            )

      end if
    end if

    ! All preparations done, we now have an array holding the original
    ! polynomial data for all requested variables and an accompanying mesh
    ! with the boundary properties.
    ! We now can go on and do the actual refinement:
    refining: do iRefLevel=1,me%max_nlevels

      lastrefine = ( iRefLevel == me%max_nLevels )
      newElems = 0

      do iScalar=1,nScalars
        call ply_sampling_var_compute_elemdev( var       = var(iScalar),    &
          &                                    threshold = me%eps_osci,     &
          &                                    min_mean  = minmean(iScalar) )
      end do

      if (allocated(reached_limit)) deallocate(reached_limit)
      allocate(reached_limit(curmesh%nElems))
      if (allocated(is_varying)) deallocate(is_varying)
      allocate(is_varying(curmesh%nElems))
      refinedElems = 0

      write(logunit(5),*) 'Adaptive sampling refinement ', iRefLevel
      write(logunit(5),*) '        Parent local mesh has ', curmesh%nelems, &
        &                 ' elements.'
      ! The decision whether an element needs to be refined or not depends
      ! on all variable components. Thus, this needs to be checked once
      ! beforehand.
      ! For each element we also determine whether this will be the last
      ! refinement, as in this case the data will be reduced to a single dof
      ! later on independent of the variation of the data.
      do iElem=1,curmesh%nElems

        elemlevel = tem_levelOf(curmesh%treeID(iElem))
        reached_limit(iElem) = lastrefine                                     &
          &                    .or. ( (me%AbsUpperBoundLevel > 0)             &
          &                    .and. (elemlevel >= me%AbsUpperBoundLevel - 1) )

        ! Now get the spectral variation (sum of absolute values of all
        ! higher modes).
        is_varying(iElem) = .false.
        ! Only consider refining, if absupperboundlevel is not reached.
        if ( (elemlevel < me%AbsUpperBoundLevel) &
          & .or. (me%AbsUpperBoundLevel == 0)    ) then
          do iScalar=1,nScalars
            is_varying(iElem) = is_varying(iElem) &
              &                 .or. var(iScalar)%deviates(iElem)
          end do
        end if

        if ( is_varying(iElem) ) then
          ! At least one of the variables in the data has a variation above
          ! the threshold and we need to refine this element.
          newElems = newElems + nChildren
          refinedElems = refinedElems+1
        else
          ! No or only small variations in all variables, no need to
          ! refine, just keep the element with the first degree of freedom
          ! only.
          newElems = newElems + 1
        end if

      end do

      write(logunit(5),*) '        The new mesh will have ', newelems, &
        &                 ' elements.'
      flush(logunit(5))

      need2refine = (newElems > curmesh%nElems)

      call MPI_Allreduce( MPI_IN_PLACE,        & !sendbuf
        &                 need2refine,         & !recvbuf
        &                 1,                   & !count
        &                 MPI_LOGICAL,         & !datatype
        &                 MPI_LOR,             & !op
        &                 curmesh%global%comm, & !comm
        &                 iError               ) !ierror

      if (need2refine) then
        write(logunit(5),*) '        Need to do a refinement!'
        ! Mesh needs to be refined, and we need to create a new one.
        if (allocated(map2global)) deallocate(map2global)
        allocate(map2global(refinedElems))
        iRefElem=0
        do iElem=1,curmesh%nElems
          if (is_varying(iElem)) then
            iRefElem = iRefElem+1
            map2global(iRefElem) = iElem
          end if
        end do
        call tem_subtree_from( me         = refine_subtree,     &
          &                    map2global = map2global,         &
          &                    comm       = curmesh%global%comm )

        if (associated(oldmesh)) then
          deallocate(oldmesh%global%property)
          deallocate(oldmesh%property)
          deallocate(oldmesh)
        end if
        oldmesh => curmesh
        nullify(curmesh)

        if (associated(oldbcs)) deallocate(oldbcs)
        oldbcs => curbcs
        nullify(curbcs)

        allocate(curmesh)
        allocate(curbcs)

        call tem_refine_global_subtree( orig_mesh       = oldmesh,        &
          &                             orig_bcs        = oldbcs,         &
          &                             subtree         = refine_subtree, &
          &                             ndims           = ndims,          &
          &                             new_mesh        = curmesh,        &
          &                             new_bcs         = curbcs,         &
          &                             restrict_to_sub = .false.         )

        call tem_destroy_subtree(refine_subtree)
      end if

      ! Apply the same reduction factor for all variables.
      ! With adaptive dof reduction this might be changed below.
      reduction_factor = me%dofReducFactor

      do iScalar=1,nScalars
        maxdeg(iScalar) = maxval(var(iScalar)%degree)
        ReducableElems(iScalar) = var(iScalar)%nDeviating*nChildren
      end do

      if (me%reduction_mode == redux_factor &
        & .and. me%adaptiveDofReduction     &
        & .and. need2refine                 ) then
        ! If adaptive reduction is active, we may increase the factor and
        ! keep more dofs after the refinement for improved accuracy without
        ! increased memory. This is computed individually for each variable
        ! as they all live in separate arrays.
        ! We only need to compute this factor, if there are actually elements
        ! to refine, otherwise the target polynomial will always have a degree
        ! of 0.
        do iScalar=1,nScalars
          ! The maximal amount of memory required will be:
          ! `totaldofs = newelems + reducableElems*((r*(maxdeg+1))**nDims - 1)`
          ! Where `r` is the dof reduction factor. Solving for this factor we
          ! get:
          if ( (newelems <= origsize(iScalar))     &
            &  .and. (reducableElems(iScalar) > 0) ) then
            memprefac = ( real(1 + origsize(iScalar) - newelems, kind=rk) &
              &               / ReducableElems(iScalar) )**(1.0_rk/nDims) &
              &         / real(maxdeg(iScalar)+1, kind=rk)
            ! Do not increase the number of dofs, even if memory would allow it.
            ! (limit factor to 1).
            memprefac = min(memprefac, 1.0_rk)
          else
            ! Even a single degree of freedom per element exceeds the original
            ! memory. We can not increase the factor and maintain more degrees
            ! of freedom without increasing the memory demand, thus we stick to
            ! the configured reduction factor.
            memprefac = reduction_factor(iScalar)
          end if

          ! Reduction factor has to be at least as large as given by the user,
          ! but if possible we'll use a larger factor and preserve more modes
          ! after the refinement.
          reduction_factor(iScalar) = max( memprefac,                &
            &                              reduction_factor(iScalar) )
        end do
      end if

      ! Move data of previous iteration to maintain variable names
      ! (var is moved to prev and the old content of prev is discarded)
      call ply_sampling_var_move( source      = var, &
        &                         destination = prev )
      allocate(var(nScalars))

      variables: do iScalar=1,nScalars

        ! The following code was moved in front of the condition to silence a
        ! compiler warning about a potentially uninitialized maxtarget variable.
        !
        ! No refinement to be done, just a single degree of freedom per
        ! element needed.
        maxtarget = 1
        containersize = newelems

        ! Allocate an array of sufficient size for the refined data.
        if (need2refine) then
          select case(me%reduction_mode)
          case(redux_factor)
            maxtarget = ceiling( reduction_factor(iScalar) &
              &                  * (maxdeg(iScalar)+1) )
          case(redux_decrement)
            maxtarget = max(maxdeg(iScalar) - me%dof_decrement + 1, 1)
          end select
          containersize = newelems &
            &           + reducableElems(iScalar) * (maxtarget**nDims - 1)
        end if

        call ply_sampling_var_allocate( var     = var(iScalar), &
          &                             nElems  = newElems,     &
          &                             datalen = containersize )

        iNewElem = 1
        var(iScalar)%first(1) = 1

        if (.not.need2refine) then
          ! No need to refine the mesh, just copy the first degree of freedom
          ! in each element and skip to next variable.
          do iElem=1,curmesh%nElems
            firstdof = var(iScalar)%first(iNewElem)
            oldfirst = prev(iScalar)%first(iElem)
            var(iScalar)%first(iNewElem+1) = firstdof + 1
            var(iScalar)%dat(firstdof) = prev(iScalar)%dat(oldfirst)
            iNewElem = iNewElem+1
          end do
          CYCLE variables
        end if

        ! Iterate over the elements of the old mesh to refine each one as
        ! needed.
        do iElem=1,oldmesh%nElems

          oldfirst = prev(iScalar)%first(iElem)
          firstdof = var(iScalar)%first(iNewElem)

          varelem: if ( is_varying(iElem) ) then
            ! There is at least one variable that varies in this element,
            ! need to create all child elements.

            if (prev(iScalar)%deviates(iElem)) then
              ! The data of this scalar varies, we need to project it to the
              ! child elements.

              if (reached_limit(iElem)) then
                ! If we reached the limit for refinements, we only keep one
                ! degree of freedom.
                targetdeg = 0
              else
                select case(me%reduction_mode)
                case(redux_factor)
                  targetdeg = ceiling( reduction_factor(iScalar)             &
                    &                  * (prev(iScalar)%degree(iElem)+1) ) - 1
                case(redux_decrement)
                  targetdeg = prev(iScalar)%degree(iElem) - me%dof_decrement
                end select
                targetdeg = max(targetdeg, 0)
              end if

              oldlast = prev(iScalar)%first(iElem+1) - 1
              lastdof = var(iScalar)%first(iNewElem) &
                &       - 1 + nChildren              &
                &             * (targetdeg+1)**nDims

              ndofs = (targetdeg+1)**ndims
              nOlddofs = oldlast - oldfirst + 1

              parent_data(1:nOlddofs,1:1)                 &
                & => prev(iScalar)%dat(oldfirst:oldlast)
              child_data(1:ndofs,1:nChildren) &
                & => var(iScalar)%dat(firstdof:lastdof)

              if (associated(filtering)) then
                call filtering(                                      &
                  &    me             = me%filter_element,           &
                  &    element_degree = prev(iScalar)%degree(iElem), &
                  &    element_data   = parent_data                  )
              end if

              call split_element( parent_degree    = prev(iScalar)        &
                &                                    %degree(iElem),      &
                &                 child_degree     = targetdeg,           &
                &                 ignore_highmodes = me%ignore_highmodes, &
                &                 parent_data      = parent_data,         &
                &                 child_data       = child_data           )

            else
              ! This scalar does not vary more than the threshold, just keep
              ! the first degree of freedom.

              ndofs = 1
              targetdeg = 0
              lastdof = firstdof+nChildren-1
              var(iScalar)%dat(firstdof:lastdof) = prev(iScalar)%dat(oldfirst)

            end if

            do iChild=1,nChildren
              ! todo: Filter modes
              ! if (filter_tolerance > 0 .and. targetdeg > 1) then
              !   ! only filter highest modes if there is a tolerance and there
              !   ! is more than one mode to filter.
              !   sum abs of modes from last backwards until filter tolerance
              !   is reached and cut off modes above that.
              ! end if
              var(iScalar)%degree(iNewElem) = targetdeg
              var(iScalar)%first(iNewElem+1) = firstdof + iChild*ndofs
              iNewElem = iNewElem+1
            end do

          else varelem
            ! No variation in this element, keep it with a single degree of
            ! freedom.

            var(iScalar)%degree(iNewElem) = 0
            var(iScalar)%dat(firstdof) = prev(iScalar)%dat(oldfirst)
            var(iScalar)%first(iNewElem+1) = firstdof + 1
            iNewElem = iNewElem+1

          end if varelem

        end do

      end do variables

      if (need2refine) then

        if (.not. trackInst%subtree%useGlobalMesh) then
          ! Mesh was refined, but we do not track the complete mesh, restrict
          ! data to the tracked elements again.
          call tem_create_subTree_of( inTree  = curmesh,             &
            &                         bc_prop = curbcs,              &
            &                         subtree = tracked_subtree,     &
            &                         inShape = trackConfig%geometry )

          ! Move data of unrestricted mesh to store the restricted data in var.
          call ply_sampling_var_move( source      = var, &
            &                         destination = prev )
          allocate(var(nScalars))

          newelems = tracked_subtree%nElems
          do iScalar=1,nScalars
            containersize = sum( (prev(iScalar)                     &
              &                   %degree(tracked_subtree           &
              &                           %map2global) + 1 )**nDims )
            call ply_sampling_var_allocate( var     = var(iScalar), &
              &                             nElems  = newElems,     &
              &                             datalen = containersize )
            var(iScalar)%first(1) = 1
            do iNewElem=1,newelems
              iElem = tracked_subtree%map2global(iNewElem)
              var(iScalar)%degree(iNewElem) = prev(iScalar)%degree(iElem)
              var(iScalar)%first(iNewElem+1)             &
                & = var(iScalar)%first(iNewElem)         &
                & + (prev(iScalar)%degree(iElem)+1)**nDims

              firstdof = var(iScalar)%first(iNewElem)
              lastdof = var(iScalar)%first(iNewElem+1) - 1
              oldfirst = prev(iScalar)%first(iElem)
              oldlast = prev(iScalar)%first(iElem+1) - 1

              var(iScalar)%dat(firstdof:lastdof) &
                & = prev(iScalar)%dat(oldfirst:oldlast)
            end do
          end do

          if (associated(oldmesh)) then
            deallocate(oldmesh%global%property)
            deallocate(oldmesh%property)
            deallocate(oldmesh)
          end if
          oldmesh => curmesh
          nullify(curmesh)

          if (associated(oldbcs)) deallocate(oldbcs)
          oldbcs => curbcs
          nullify(curbcs)

          allocate(curmesh)
          allocate(curbcs)

          call tem_create_tree_from_sub( intree  = oldmesh,         &
            &                            subtree = tracked_subtree, &
            &                            newtree = curmesh          )

        end if

      else

        ! Nothing to refine anymore, leave the loop.
        EXIT refining

      end if

    end do refining

    ! Now the refined data is stored in var%dat, and there is only one degree
    ! of freedom left for each element.
    ! The final mesh is stored in curmesh.

    ! Discard old data to free memory.
    if (associated(oldmesh)) then
      deallocate(oldmesh%global%property)
      deallocate(oldmesh%property)
      deallocate(oldmesh)
    end if
    if (associated(oldbcs)) deallocate(oldbcs)

    call ply_sampling_var_move( source      = var, &
      &                         destination = prev )

    get_element => get_sampled_element
    get_params => tem_varSys_getparams_dummy
    set_params => tem_varSys_setparams_dummy
    nullify(get_point)
    nullify(setup_indices)
    nullify(get_valOfIndex)

    ! Finally link the data into a structure that is suitable for the method
    ! data of each variable
    nVars = trackInst%varmap%varPos%nVals
    iScalar = 0
    do iVar=1,nVars

      varpos = trackInst%varmap%varpos%val(iVar)
      nComponents = varsys%method%val(varpos)%nComponents
      allocate(vardat)
      allocate(vardat%component(nComponents))
      do iComponent=1,nComponents
        iScalar = iScalar+1
        vardat%component(iComponent)%dat => prev(iScalar)%dat
      end do

      call tem_varSys_append_stateVar(                 &
        & me             = resvars,                    &
        & varname        = varsys%varname%val(varpos), &
        & nComponents    = nComponents,                &
        & method_data    = c_loc(vardat),              &
        & set_params     = set_params,                 &
        & get_point      = get_point,                  &
        & get_element    = get_element,                &
        & get_params     = get_params,                 &
        & setup_indices  = setup_indices,              &
        & get_valofindex = get_valofindex              )

      nullify(vardat)

    end do

    do iScalar=1,nScalars
      nullify(prev(iScalar)%dat)
      deallocate(prev(iScalar)%degree)
      deallocate(prev(iScalar)%first)
    end do
    deallocate(prev)

    new_mesh = curmesh

  end subroutine ply_sample_adaptive
  ! ------------------------------------------------------------------------- !
  ! ------------------------------------------------------------------------- !


  ! ------------------------------------------------------------------------- !
  !> Get sampled data.
  !!
  !! This routine provides the get_element function of the variable definition
  !! to access the sampled data array obtained by ply_sample_data.
  subroutine get_sampled_element( fun, varsys, elempos, time, tree, n, &
    &                             nDofs, res                           )
    ! -------------------------------------------------------------------- !
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of elements to obtain for this variable (vectorized access).
    integer, intent(in) :: n

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! -------------------------------------------------------------------- !
    type(sampled_method_data_type), pointer :: p
    integer :: iElem
    integer :: iComp
    integer :: nComps
    ! -------------------------------------------------------------------- !
    nComps = fun%nComponents

    call c_f_pointer(fun%method_data, p)

    do iComp=1,nComps
      do iElem=1,n
        res(iComp+(iElem-1)*nComps) &
          & = p%component(iComp)%dat(elempos(iElem))
      end do
    end do

  end subroutine get_sampled_element
  ! ------------------------------------------------------------------------- !
  ! ------------------------------------------------------------------------- !

end module ply_sampling_adaptive_module
