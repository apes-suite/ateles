! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

!> author: Jens Zudrop
!! Routines for stabilization of a spectral scheme.
module atl_stabilize_module
  use env_module,                    only: rk
  ! Treelm modules
  use tem_aux_module,                only: tem_abort
  use tem_logging_module,            only: logUnit
  use tem_element_module,            only: eT_fluid
  use tem_general_module,            only: tem_general_type
  use tem_timer_module,              only: tem_startTimer, tem_stopTimer
  use treelmesh_module,              only: treelmesh_type

  ! Ateles modules
  use atl_cube_elem_module,          only: atl_cube_elem_type
  use atl_kerneldata_module,         only: atl_statedata_type
  use atl_scheme_module,             only: atl_scheme_type,        &
    &                                      atl_modg_scheme_prp,    &
    &                                      atl_modg_2d_scheme_prp, &
    &                                      atl_modg_1d_scheme_prp
  use atl_stabilization_module,      only: atl_no_stab_prp,                 &
    &                                      atl_spectral_visc_prp,           &
    &                                      atl_positivity_preserv_prp,      &
    &                                      atl_cons_positivity_preserv_prp, &
    &                                      atl_cheb_spectral_visc_prp,      &
    &                                      atl_covolume_prp
  use atl_spectral_viscosity_module, only: atl_spectral_visc_type,    &
    &                                      atl_exp_spectral_visc_prp, &
    &                                      atl_poly_spectral_visc_prp
  use atl_positivity_preserv_module, only: atl_positivity_preserv_type
  use atl_cons_positivity_preserv_module, &
    &                                only: atl_cons_positivity_preserv_type
  use atl_equation_module,           only: atl_equations_type
  use atl_eqn_euler_module,          only: atl_euler_type
  use atl_boundary_module,           only: atl_level_boundary_type, &
    &                                      atl_get_numBndElems
  use atl_bc_header_module,          only: atl_boundary_type
  use atl_covolume_module,           only: atl_covolume_type
  use atl_covolume_boundary_module,  only: atl_set_covolume_bnd
  use atl_covolume_projection_module, only: atl_primal_to_covolume_projection, &
                                          & atl_covolume_to_primal_projection, &
                                          & atl_primal_to_covolume_projection_2d, &
                                          & atl_covolume_to_primal_projection_2d, &
                                          & atl_primal_to_covolume_projection_1d, &
                                          & atl_covolume_to_primal_projection_1d
  use atl_modg_multiLevel_module,    only: atl_modg_fineToCoarseElem, &
                                         & atl_modg_coarseToFineElem
  use atl_modg_2d_multiLevel_module, only: atl_modg_fineToCoarseElem_2d, &
                                         & atl_modg_coarseToFineElem_2d
  use atl_modg_1d_multiLevel_module, only: atl_modg_fineToCoarseElem_1d, &
                                         & atl_modg_coarseToFineElem_1d
  use atl_timer_module,              only: atl_timerHandles
  use atl_materialPrp_module,        only: atl_material_type

  ! Polynomials modules
  use ply_poly_project_module,       only: ply_poly_project_type, &
    &                                      ply_poly_project_n2m,  &
    &                                      ply_poly_project_m2n
  use ply_oversample_module,         only: ply_convert2oversample,  &
    &                                      ply_convertFromoversample
  use ply_polyBaseExc_module,        only: ply_fpt_exec_striped, &
    &                                      assignment(=)
  use ply_dof_module,                only: ply_change_poly_space, &
    &                                      Q_space,               &
    &                                      P_space

  implicit none
  private

  type atl_adaptive_orders_type
    real(kind=rk), allocatable :: orders(:)
  end type atl_adaptive_orders_type

  public :: atl_stabilize, atl_positivity_preserv, &
    &       atl_positivity_preserv_2d,             &
    &       atl_cons_positivity_preserv,           &
    &       atl_cons_positivity_preserv_2d,        &
    &       atl_primal_to_covolume_projection,     &
    &       atl_covolume_to_primal_projection,     &
    &       atl_primal_to_covolume_projection_2d,  &
    &       atl_covolume_to_primal_projection_2d,  &
    &       atl_primal_to_covolume_projection_1d,  &
    &       atl_covolume_to_primal_projection_1d

contains


  !> Subroutine to apply the stabilization procedure to the state vector.
  subroutine atl_stabilize( minlevel, maxlevel, statedata_list,              &
    & statedata_stab_list, mesh_list, scheme_list, poly_proj_pos,            &
    & poly_proj_list, equation, tree, bc, boundary, general, commStateTimer, &
    & material_list)
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: minlevel
    integer, intent(in) :: maxlevel
    type(atl_statedata_type), intent(inout) :: statedata_list(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: statedata_stab_list(minlevel:maxlevel)
    type(atl_cube_elem_type), intent(inout) :: mesh_list(minlevel:maxlevel)
    type(atl_scheme_type), intent(inout) :: scheme_list(minlevel:maxlevel)
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    type(atl_material_type), intent (in) :: material_list(minlevel:maxlevel)
    ! --------------------------------------------------------------------------
    integer :: iLevel, iStab, iDir
    integer :: nStabs, nElems, nScalars, nDofs, nTotal
    real(kind=rk) :: ref_order, recovery_order
    real(kind=rk) :: pressure, mach
    real(kind=rk), allocatable :: velocity(:)
    integer :: oversamp_dofs, iElem, nquadpoints, iPoint
    type(atl_adaptive_orders_type) :: adaptive_orders(minlevel:maxlevel)
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    real(kind=rk), allocatable :: pointVal(:,:)
    type(atl_statedata_type), allocatable :: Q_statedata_list(:)
    type(atl_statedata_type), allocatable :: Q_statedata_stab_list(:,:)
    integer :: nBndStabElems(minLevel:maxLevel,1:3)
    integer :: maxPolyDeg
    ! --------------------------------------------------------------------------

    call tem_startTimer( timerHandle = atl_timerHandles%stabalize )
    do iLevel = minlevel,maxlevel
      nElems = mesh_list(ilevel)%descriptor%elem%nElems(eT_fluid)

      do iElem = 1, nElems
        if (material_list(iLevel)%material_dat%mode_reducable(iElem)) then
          statedata_list(iLevel)%state(iElem,2:,:) = 0.0_rk
        end if
      end do
      if(any(scheme_list(ilevel)%stabilization(:)%spectral_visc%isAdaptive) .or. &
        & any(scheme_list(ilevel)%stabilization(:)%covolume%isAdaptive)) then

        allocate( adaptive_orders(iLevel)%orders( nElems ) )

        if( scheme_list(iLevel)%stabilization(1)%stab_kind .ne. atl_spectral_visc_prp ) then
          write(logUnit(1),*) 'ERROR: adaptive stabilization: first filter has to be '
          write(logUnit(1),*) '       of spectral viscosity type, stopping ...'
          call tem_abort()
        end if
        ref_order = scheme_list(iLevel)%stabilization(1)%spectral_visc%order
        recovery_order = scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_order

        ! Initialize all elements with the reference order
        adaptive_orders(iLevel)%orders(:) = ref_order

        select case(scheme_list(ilevel)%scheme)
        case(atl_modg_1d_scheme_prp)

          oversamp_dofs = poly_proj_list(poly_proj_pos(iLevel))%body_1D%oversamp_dofs
          nquadpoints = poly_proj_list(poly_proj_pos(iLevel))%body_1D%nquadpoints
          allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars) )
          allocate( pointVal(nQuadPoints,equation%varSys%nScalars) )
          allocate( velocity(1) )

          do iElem = 1, nElems

            call ply_convert2oversample(state= statedata_list(ilevel)%state(iElem,:,:),  &
              &                         ndim = 1,                                         &
              &                         poly_proj= poly_proj_list(poly_proj_pos(iLevel)), &
              &                         modalCoeffs = modalCoeffs           )
            call ply_poly_project_m2n(me = poly_proj_list(poly_proj_pos(iLevel)), &
             &                       dim = 1 ,                         &
             &                       nVars = equation%varSys%nScalars, &
             &                       nodal_data= pointVal,             &
             &                       modal_data= modalCoeffs           )

            select case(equation%eq_kind)
            case('euler_1d', 'lineareuler_1d')

               do iPoint = 1, nquadpoints
                 pressure = (equation%euler%isen_coef-1.0_rk)*( pointVal(iPoint,4) &
                          & -0.5_rk * pointVal(iPoint,2)**2 /pointVal(iPoint,1) )
                 velocity(1) = pointVal(iPoint,2)/pointVal(iPoint,1)
                 mach = sqrt( velocity(1)**2 ) &
                      & / ( sqrt( equation%euler%isen_coef * pressure / pointVal(iPoint,1)  ) )

                 ! Check how critical the state is.
                 if( pressure .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_pressure &
                   & .or. mach .gt. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_mach &
                   & .or. pointVal(iPoint,1) .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_density &
                   & ) then
                   adaptive_orders(iLevel)%orders(iElem) = recovery_order
                 end if
               end do
            case default
              write(logUnit(1),*) 'ERROR in adaptive strategy of 1D atl_stabilize_module: ' // &
                & 'Unknown equation, stopping ...'
              call tem_abort()
            end select
          end do

          deallocate( modalCoeffs )
          deallocate( pointVal )
          deallocate( velocity )

        case(atl_modg_2d_scheme_prp)

          oversamp_dofs = poly_proj_list(poly_proj_pos(iLevel))%body_2D%oversamp_dofs
          nquadpoints = poly_proj_list(poly_proj_pos(iLevel))%body_2D%nquadpoints
          allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars) )
          allocate( pointVal(nQuadPoints,equation%varSys%nScalars) )
          allocate( velocity(2) )

          do iElem = 1, nElems

            call ply_convert2oversample(state= statedata_list(ilevel)%state(iElem,:,:),  &
              &                         ndim = 2,                                         &
              &                         poly_proj= poly_proj_list(poly_proj_pos(iLevel)), &
              &                         modalCoeffs = modalCoeffs           )
            call ply_poly_project_m2n(me = poly_proj_list(poly_proj_pos(iLevel)), &
             &                       dim = 2 ,                         &
             &                       nVars = equation%varSys%nScalars, &
             &                       nodal_data= pointVal,             &
             &                       modal_data= modalCoeffs           )

            select case(equation%eq_kind)
            case('euler_2d', 'navier_stokes_2d', 'lineareuler_2d')

               do iPoint = 1, nquadpoints
                 pressure = (equation%euler%isen_coef-1.0_rk)*( pointVal(iPoint,4) &
                          & -0.5_rk * sum(pointVal(iPoint,2:3)**2)/pointVal(iPoint,1) )
                 velocity(1:2) = pointVal(iPoint,2:3)/pointVal(iPoint,1)
                 mach = sqrt( sum(velocity(1:2)**2) ) &
                      & / ( sqrt( equation%euler%isen_coef * pressure / pointVal(iPoint,1)  ) )

                 ! Check how critical the state is.
                 if( pressure .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_pressure &
                   & .or. mach .gt. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_mach &
                   & .or. pointVal(iPoint,1) .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_density &
                   & ) then
                   adaptive_orders(iLevel)%orders(iElem) = recovery_order
                 end if
               end do
            case('filtered_navier_stokes_2d')

               do iPoint = 1, nquadpoints
                 pressure = (equation%euler%isen_coef-1.0_rk)*( pointVal(iPoint,5) &
                         & -0.5_rk * sum(pointVal(iPoint,2:3)**2)/pointVal(iPoint,1) &
                         & - pointVal(iPoint,5) )
                 velocity(1:2) = pointVal(iPoint,2:3)/pointVal(iPoint,1)
                 mach = sqrt( sum(velocity(1:2)**2) ) &
                    & / ( sqrt( equation%euler%isen_coef * pressure / pointVal(iPoint,1)  ) )

                 ! Check how critical the state is.
                 if( pressure .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_pressure &
                   & .or. mach .gt. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_mach &
                   & .or. pointVal(iPoint,1) .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_density &
                   & ) then
                   adaptive_orders(iLevel)%orders(iElem) = recovery_order
                 end if
              end do

            case default
              write(logUnit(1),*) 'ERROR in adaptive strategy of 2D atl_stabilize_module: ' // &
                & 'Unknown equation, stopping ...'
              call tem_abort()
            end select
          end do

          deallocate( modalCoeffs )
          deallocate( pointVal )
          deallocate( velocity )

        case(atl_modg_scheme_prp)

          oversamp_dofs = poly_proj_list(poly_proj_pos(iLevel))%body_3D%oversamp_dofs
          nquadpoints = poly_proj_list(poly_proj_pos(iLevel))%body_3D%nquadpoints
          allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars) )
          allocate( pointVal(nQuadPoints,equation%varSys%nScalars) )
          allocate( velocity(3) )

          do iElem = 1, nElems

            call ply_convert2oversample(state= statedata_list(ilevel)%state(iElem,:,:),  &
              &                         ndim = 3,                                         &
              &                         poly_proj= poly_proj_list(poly_proj_pos(iLevel)), &
              &                         modalCoeffs = modalCoeffs           )
            call ply_poly_project_m2n(me = poly_proj_list(poly_proj_pos(iLevel)), &
             &                       dim = 3 ,                         &
             &                       nVars = equation%varSys%nScalars, &
             &                       nodal_data= pointVal,             &
             &                       modal_data= modalCoeffs           )

            select case(equation%eq_kind)
            case('euler', 'navier_stokes')

               do iPoint = 1, nquadpoints
                 pressure = (equation%euler%isen_coef-1.0_rk)*( pointVal(iPoint,5) &
                         & -0.5_rk * sum(pointVal(iPoint,2:4)**2)/pointVal(iPoint,1) )
                 velocity(1:3) = pointVal(iPoint,2:4)/pointVal(iPoint,1)
                 mach = sqrt( sum(velocity(1:3)**2) ) &
                    & / ( sqrt( equation%euler%isen_coef * pressure / pointVal(iPoint,1)  ) )

                 ! Check how critical the state is.
                 if( pressure .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_pressure &
                   & .or. mach .gt. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_mach &
                   & .or. pointVal(iPoint,1) .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_density &
                   & ) then
                   adaptive_orders(iLevel)%orders(iElem) = recovery_order
                 end if
              end do

            case('filtered_navier_stokes')

               do iPoint = 1, nquadpoints
                 pressure = (equation%euler%isen_coef-1.0_rk)*( pointVal(iPoint,5) &
                         & -0.5_rk * sum(pointVal(iPoint,2:4)**2)/pointVal(iPoint,1) &
                         & - pointVal(iPoint,6) )
                 velocity(1:3) = pointVal(iPoint,2:4)/pointVal(iPoint,1)
                 mach = sqrt( sum(velocity(1:3)**2) ) &
                    & / ( sqrt( equation%euler%isen_coef * pressure / pointVal(iPoint,1)  ) )

                 ! Check how critical the state is.
                 if( pressure .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_pressure &
                   & .or. mach .gt. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_mach &
                   & .or. pointVal(iPoint,1) .le. &
                   & scheme_list(iLevel)%stabilization(1)%spectral_visc%recovery_density &
                   & ) then
                   adaptive_orders(iLevel)%orders(iElem) = recovery_order
                 end if
              end do
            case default
              write(logUnit(1),*) 'ERROR in adaptive strategy of 3D atl_stabilize_module: ' // &
                & 'Unknown equation, stopping ...'
              call tem_abort()
            end select
          end do

          deallocate( modalCoeffs )
          deallocate( pointVal )
          deallocate( velocity )

        case default
          write(logUnit(1),*) 'ERROR in adaptive strategy of atl_stabilize_module: ' // &
            & 'Unknown scheme, stopping ...'
          call tem_abort()
        end select


      end if
    end do

    nStabs = size(scheme_list(minlevel)%stabilization)
    do iStab = 1, nStabs


      select case(scheme_list(minlevel)%scheme)
      case(atl_modg_1d_scheme_prp)
        select case(scheme_list(minlevel)%stabilization(iStab)%stab_kind)
        case(atl_no_stab_prp)
          ! nothing to be done
        case(atl_spectral_visc_prp)
          do iLevel = minlevel, maxlevel
            call atl_spectral_visc_1d(                               &
            & state      = statedata_list(iLevel),                   &
            & mesh       = mesh_list(iLevel),                        &
            & filter     = scheme_list(iLevel)%stabilization(iStab)  &
            &                             %spectral_visc,            &
            & maxPolyDeg = scheme_list(iLevel)%modg_1d%maxPolyDegree )
          end do
        case(atl_cheb_spectral_visc_prp)
          do iLevel = minlevel, maxlevel
            call atl_cheb_spectral_visc_1d(                          &
            & state      = statedata_list(iLevel),                   &
            & mesh       = mesh_list(iLevel),                        &
            & filter     = scheme_list(iLevel)%stabilization(iStab)  &
            &                                 %spectral_visc,        &
            & poly_proj  = poly_proj_list(poly_proj_pos(iLevel)),    &
            & maxPolyDeg = scheme_list(iLevel)%modg_1d%maxPolyDegree )
          end do
        case(atl_covolume_prp)
          call atl_covolume_1d(                                    &
          & minlevel   = minlevel,                                 &
          & maxlevel   = maxlevel,                                 &
          & state      = statedata_list,                           &
          & state_stab = statedata_stab_list,                      &
          & mesh       = mesh_list,                                &
          & filter     = scheme_list(minlevel)%stabilization(iStab)&
          &                             %covolume,                 &
          & scheme = scheme_list,                                  &
          & equation = equation,                                   &
          & tree   = tree,                                         &
          & poly_proj  = poly_proj_list,                           &
          & poly_proj_pos = poly_proj_pos,                         &
          & bc = bc,                                               &
          & boundary = boundary,                                   &
          & general    = general,                                  &
          & commStateTimer = commStateTimer                        )
        case default
          write(logUnit(1),*) 'ERROR in atl_stabilize_module: Unknown' // &
            & ' stabilization for 1D MODG scheme, stopping ...'
          call tem_abort()
        end select
      case(atl_modg_2d_scheme_prp)
        select case(scheme_list(minlevel)%stabilization(iStab)%stab_kind)

        case(atl_no_stab_prp)
          ! nothing to be done

        case(atl_spectral_visc_prp)
          do iLevel = minlevel, maxlevel
            if(scheme_list(minlevel)%stabilization(iStab)%spectral_visc%isAdaptive) then
              call atl_spectral_visc_2d(                               &
              & state      = statedata_list(iLevel),                   &
              & mesh       = mesh_list(iLevel),                        &
              & filter     = scheme_list(iLevel)%stabilization(iStab)  &
              &                                 %spectral_visc,        &
              & orders     = adaptive_orders(iLevel)%orders,           &
              & maxPolyDeg = scheme_list(iLevel)%modg_2d%maxPolyDegree )
            else
              call atl_spectral_visc_2d(                               &
              & state      = statedata_list(iLevel),                   &
              & mesh       = mesh_list(iLevel),                        &
              & filter     = scheme_list(iLevel)%stabilization(iStab)  &
              &                                 %spectral_visc,        &
              & maxPolyDeg = scheme_list(iLevel)%modg_2d%maxPolyDegree )
            end if
          end do

        case(atl_positivity_preserv_prp)
          if(.not.(trim(equation%eq_kind).eq.'euler_2d') .and. &
            & .not.(trim(equation%eq_kind).eq.'navier_stokes_2d')) then
            write(logUnit(1),*) 'ERROR in atl_stabilize_module: '     // &
              & 'Positivity preserving stablization works only with ' // &
              & 'Euler_2d and Navier_Stokes_2d, stopping ...'
            call tem_abort()
          end if
          do iLevel = minlevel, maxlevel
            if(.not. poly_proj_list(poly_proj_pos(iLevel))%lobattoPoints) then
              write(logUnit(1),*) 'ERROR in atl_stabilize_module: '     // &
                & 'Positivity preserving stablization works only with ' // &
                & 'Lobatto points, stopping ...'
              call tem_abort()
            end if
            call atl_positivity_preserv_2d(                          &
              & state     = statedata_list(iLevel),                  &
              & mesh      = mesh_list(iLevel),                       &
              & filter    = scheme_list(iLevel)%stabilization(iStab) &
              &                                %positivity_preserv,  &
              & euler     = equation%euler,                          &
              & poly_proj = poly_proj_list(poly_proj_pos(iLevel) )   )
          end do

        case(atl_cons_positivity_preserv_prp)
          select case(trim(equation%eq_kind))
          case('euler_2d','navier_stokes_2d')
            do iLevel = minlevel, maxlevel
              if(.not.poly_proj_list(poly_proj_pos(iLevel))%lobattoPoints) then
                write(logUnit(1),*) 'ERROR in atl_stabilize_module: ' // &
                  & 'Conservative positivity preserving stablization' // &
                  & ' works only with Lobatto points, stopping ...'
                call tem_abort()
              end if
              call atl_cons_positivity_preserv_2d(                         &
                & state     = statedata_list(iLevel),                      &
                & mesh      = mesh_list(iLevel),                           &
                & filter    = scheme_list(iLevel)%stabilization(iStab)     &
                &                                %cons_positivity_preserv, &
                & euler     = equation%euler,                              &
                & poly_proj = poly_proj_list(poly_proj_pos(iLevel))        )
            end do
          case default
              write(logUnit(1),*) 'ERROR in atl_stabilize_module: ' // &
                & 'Conservative positivity preserving stablization' // &
                & ' does not work with '//trim(equation%eq_kind)       &
                & //', stopping ...'
              call tem_abort()
          end select

        case(atl_covolume_prp)
          if(scheme_list(minlevel)%stabilization(iStab)%covolume%isAdaptive) then
            call atl_covolume_2d(                                    &
            & minlevel   = minlevel,                                 &
            & maxlevel   = maxlevel,                                 &
            & state      = statedata_list,                           &
            & state_stab = statedata_stab_list,                      &
            & mesh       = mesh_list,                                &
            & filter     = scheme_list(minlevel)%stabilization(iStab)&
            &                             %covolume,                 &
            & adaptive_orders = adaptive_orders,                     &
            & scheme = scheme_list,                                  &
            & equation = equation,                                   &
            & tree           = tree,                                 &
            & poly_proj  = poly_proj_list,                           &
            & poly_proj_pos = poly_proj_pos,                         &
            & bc = bc,                                               &
            & boundary = boundary,                                   &
            & general    = general,                                  &
            & commStateTimer = commStateTimer                        )
          else
            call atl_covolume_2d(                                    &
            & minlevel   = minlevel,                                 &
            & maxlevel   = maxlevel,                                 &
            & state      = statedata_list,                           &
            & state_stab = statedata_stab_list,                      &
            & mesh       = mesh_list,                                &
            & filter     = scheme_list(minlevel)%stabilization(iStab)&
            &                             %covolume,                 &
            & scheme = scheme_list,                                  &
            & equation = equation,                                   &
            & tree           = tree,                                 &
            & poly_proj  = poly_proj_list,                           &
            & poly_proj_pos = poly_proj_pos,                         &
            & bc = bc,                                               &
            & boundary = boundary,                                   &
            & general    = general,                                  &
            & commStateTimer = commStateTimer                        )
          end if

        case default
          write(logUnit(1),*) 'ERROR in atl_stabilize_module: Unknown ' // &
            & 'stabilization for 2D MODG scheme, stopping ...'
          call tem_abort()
        end select

      case(atl_modg_scheme_prp)
        select case(scheme_list(minlevel)%stabilization(iStab)%stab_kind)
        case(atl_no_stab_prp)
          ! nothing to be done

        case(atl_spectral_visc_prp)
          do iLevel = minlevel, maxlevel
            if(scheme_list(minlevel)%stabilization(iStab)%spectral_visc%isAdaptive) then
              call atl_spectral_visc_3d(                                &
                & state      = statedata_list(iLevel),                  &
                & mesh       = mesh_list(iLevel),                       &
                & filter     = scheme_list(iLevel)%stabilization(iStab) &
                &                                 %spectral_visc,       &
                & orders     = adaptive_orders(iLevel)%orders,          &
                & maxPolyDeg = scheme_list(iLevel)%modg%maxPolyDegree   )
            else
              call atl_spectral_visc_3d(                                &
                & state      = statedata_list(iLevel),                  &
                & mesh       = mesh_list(iLevel),                       &
                & filter     = scheme_list(iLevel)%stabilization(iStab) &
                &                                 %spectral_visc,       &
                & maxPolyDeg = scheme_list(iLevel)%modg%maxPolyDegree   )
            end if
          end do

        case(atl_positivity_preserv_prp)
          if(.not.(trim(equation%eq_kind).eq.'euler') .and. &
            & .not.(trim(equation%eq_kind).eq.'navier_stokes')) then
            write(logUnit(1),*) 'ERROR in atl_stabilize_module: ' // &
              & 'Positivity preserving stablization '             // &
              & 'works only with Euler and Navier_Stokes, stopping ...'
            call tem_abort()
          end if
          do iLevel = minlevel, maxlevel
            if(.not.poly_proj_list(poly_proj_pos(iLevel))%lobattoPoints) then
              write(logUnit(1),*) 'ERROR in atl_stabilize_module: ' // &
                & 'Positivity preserving stablization '             // &
                & 'works only with Lobatto points, stopping ...'
              call tem_abort()
            end if
            call atl_positivity_preserv(                          &
              & state = statedata_list(iLevel),                   &
              & mesh = mesh_list(iLevel),                         &
              & filter = scheme_list(iLevel)%stabilization(iStab) &
              &                             %positivity_preserv,  &
              & euler = equation%euler,                           &
              & poly_proj = poly_proj_list(poly_proj_pos(iLevel)) )
          end do

        case(atl_cons_positivity_preserv_prp)
          if(.not.(trim(equation%eq_kind).eq.'euler') .and. &
            & .not.(trim(equation%eq_kind).eq.'navier_stokes')) then
            write(logUnit(1),*) 'ERROR in atl_stabilize_module: '  // &
              & 'Conservative positivity preserving stablization ' // &
              & 'works only with Euler and Navier_Stokes, stopping ...'
            call tem_abort()
          end if
          do iLevel = minlevel, maxlevel
            if(.not. poly_proj_list(poly_proj_pos(iLevel))%lobattoPoints) then
              write(logUnit(1),*) 'ERROR in atl_stabilize_module: '  // &
                & 'Conservative positivity preserving stablization ' // &
                & 'works only with Lobatto points, stopping ...'
              call tem_abort()
            end if
            call atl_cons_positivity_preserv(                         &
              & state = statedata_list(iLevel),                       &
              & mesh = mesh_list(iLevel),                             &
              & filter = scheme_list(iLevel)%stabilization(iStab)     &
              &                             %cons_positivity_preserv, &
              & euler = equation%euler,                               &
              & poly_proj = poly_proj_list(poly_proj_pos(iLevel))     )
          end do

        case(atl_covolume_prp)
          if(scheme_list(minlevel)%modg%basisType .eq. Q_space) then !Q_space
            if(scheme_list(minlevel)%stabilization(iStab)%covolume%isAdaptive) then
              call atl_covolume(                                         &
                & minlevel   = minlevel,                                 &
                & maxlevel   = maxlevel,                                 &
                & state      = statedata_list,                           &
                & state_stab = statedata_stab_list,                      &
                & mesh       = mesh_list,                                &
                & filter     = scheme_list(minlevel)%stabilization(iStab)&
                &                             %covolume,                 &
                & adaptive_orders = adaptive_orders,                     &
                & scheme = scheme_list,                                  &
                & equation = equation,                                   &
                & tree           = tree,                                 &
                & poly_proj  = poly_proj_list,                           &
                & poly_proj_pos = poly_proj_pos,                         &
                & bc = bc,                                               &
                & boundary = boundary,                                   &
                & general    = general,                                  &
                & commStateTimer = commStateTimer                        )
            else
              call atl_covolume(                                             &
                & minlevel       = minlevel,                                 &
                & maxlevel       = maxlevel,                                 &
                & state          = statedata_list,                           &
                & state_stab     = statedata_stab_list,                      &
                & mesh           = mesh_list,                                &
                & filter         = scheme_list(minlevel)%stabilization(iStab)&
                &                                 %covolume,                 &
                & scheme         = scheme_list,                              &
                & equation       = equation,                                 &
                & tree           = tree,                                     &
                & poly_proj      = poly_proj_list,                           &
                & poly_proj_pos  = poly_proj_pos,                            &
                & bc             = bc,                                       &
                & boundary       = boundary,                                 &
                & general        = general,                                  &
                & commStateTimer = commStateTimer                            )
            end if
          else !P_space

            ! If co-volume filter is used with P-polynomials we just
            ! project the state, which is given in P-space, to
            ! a temporary state in Q-space format and perform the
            ! filtering on this temporary state. After the filtering is done
            ! we can simply project the filtered state back to P-space.
            ! In this way we can not benefit from the performance advantages
            ! that come with P-polynomials but we can use the same subroutines
            ! for both cases.

            ! Allocate and initialize the temporary state and the
            ! stab-state in Q-space
            allocate(Q_statedata_list(minlevel:maxlevel))
            allocate(Q_statedata_stab_list(minlevel:maxlevel,3))

            nScalars = equation%varSys%nScalars

            do iLevel = minlevel, maxlevel
              maxPolyDeg = scheme_list(iLevel)%modg%maxPolyDegree
              nDofs = (maxPolyDeg+1)**3
              nElems = mesh_list(ilevel)%descriptor%elem%nElems(eT_fluid)

              allocate(Q_statedata_list(iLevel)%state( nElems,   &
                &                                      nDoFs,    &
                &                                      nScalars) )
              Q_statedata_list(iLevel)%local_time = &
                & statedata_list(maxlevel)%local_time
              Q_statedata_list(iLevel)%state = 0.0_rk

              ! Projection of the state from P-space to Q-space
              ! Higher modes in the Q-space will be filled with zeros.
              call ply_change_poly_space(                      &
                & inspace    = P_space,                        &
                & instate    = statedata_list(iLevel)%state,   &
                & outstate   = Q_statedata_list(iLevel)%state, &
                & maxPolyDeg = maxPolydeg,                     &
                & nElems     = nElems,                         &
                & nVars      = nScalars,                       &
                & nDims      = 3                               )

              nBndStabElems = atl_get_numBndElems( &
                & minLevel      = minLevel,        &
                & maxLevel      = maxLevel,        &
                & boundary_list = boundary         )

              ! Allocate the stab-state in Q-space and prepare it for the
              ! dim-by-dim procedure performed in the co-volume filtering.
              do iDir=1,3
                nTotal = mesh_list(iLevel)%faces%dimByDimDesc(iDir)%nElems &
                  &    + nBndStabElems(iLevel,iDir)

                allocate(Q_statedata_stab_list(iLevel,iDir)%state( nTotal,   &
                  &                                                nDoFs,    &
                  &                                                nScalars) )
                Q_statedata_stab_list(iLevel,iDir)%local_time = &
                  & statedata_stab_list(maxlevel)%local_time
                Q_statedata_stab_list(iLevel,iDir)%state = 0.0_rk
              end do

            end do

            ! Apply the covolume filtering to the temporary states in Q-space
            ! just like normal.
            if(scheme_list(minlevel)%stabilization(iStab)%covolume%isAdaptive) then
              call atl_covolume(                                         &
                & minlevel   = minlevel,                                 &
                & maxlevel   = maxlevel,                                 &
                & state      = Q_statedata_list,                         &
                & state_stab = Q_statedata_stab_list,                    &
                & mesh       = mesh_list,                                &
                & filter     = scheme_list(minlevel)%stabilization(iStab)&
                &                             %covolume,                 &
                & adaptive_orders = adaptive_orders,                     &
                & scheme = scheme_list,                                  &
                & equation = equation,                                   &
                & tree           = tree,                                 &
                & poly_proj  = poly_proj_list,                           &
                & poly_proj_pos = poly_proj_pos,                         &
                & bc = bc,                                               &
                & boundary = boundary,                                   &
                & general    = general,                                  &
                & commStateTimer = commStateTimer                        )
            else
              call atl_covolume(                                             &
                & minlevel       = minlevel,                                 &
                & maxlevel       = maxlevel,                                 &
                & state          = Q_statedata_list,                         &
                & state_stab     = Q_statedata_stab_list,                    &
                & mesh           = mesh_list,                                &
                & filter         = scheme_list(minlevel)%stabilization(iStab)&
                &                                 %covolume,                 &
                & scheme         = scheme_list,                              &
                & equation       = equation,                                 &
                & tree           = tree,                                     &
                & poly_proj      = poly_proj_list,                           &
                & poly_proj_pos  = poly_proj_pos,                            &
                & bc             = bc,                                       &
                & boundary       = boundary,                                 &
                & general        = general,                                  &
                & commStateTimer = commStateTimer                            )
            end if

            ! Projection of the stab-state back from Q-space to P-space.
            ! This is done by cutting off the higher modes.
            do iLevel = minlevel, maxlevel
              call ply_change_poly_space(                             &
                & inspace    = Q_space,                               &
                & instate    = Q_statedata_stab_list(iLevel,3)%state, &
                & outstate   = statedata_stab_list(iLevel)%state,     &
                & maxPolyDeg = maxPolydeg,                            &
                & nElems     = nElems,                                &
                & nVars      = nScalars,                              &
                & nDims      = 3                                      )
            end do

            deallocate(Q_statedata_list)
            deallocate(Q_statedata_stab_list)

          end if

        case default
          write(logUnit(1),*) 'ERROR in atl_stabilize_module: ' // &
            & 'Unknown stabilization for MODG scheme, stopping ...'
          call tem_abort()
        end select

      case default
        write(logUnit(1),*) 'ERROR in atl_stabilize_module: ' // &
          & 'Unknown scheme, stopping ...'
        call tem_abort()
      end select

    end do
  call tem_stopTimer( timerHandle = atl_timerHandles%stabalize )

  end subroutine atl_stabilize

  !> Apply conservative limitation of denisty and energy to ensure positivity
  !! for density and pressure.
  subroutine atl_cons_positivity_preserv( state, mesh, filter, euler, &
    &                                     poly_proj )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_cons_positivity_preserv_type), intent(in) :: filter
    type(atl_Euler_type), intent(in) :: euler
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:), pointVal(:,:), hatDens(:), &
      &                           t_vec(:), limitedPntVal(:,:)
    real(kind=rk) :: mean(5)
    ! Loop vars
    integer :: iElem, iPoint, iVar, i
    real(kind=rk) :: minDens, theta_1, theta_2, momSq, tmpPressure, &
      &              pressMean
    integer :: nquadpoints
    integer :: oversamp_dofs
    ! --------------------------------------------------------------------------
    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_3D%nquadpoints
    oversamp_dofs = poly_proj%body_3D%oversamp_dofs

    allocate( modalCoeffs(oversamp_dofs,5) )
    allocate( pointVal(nQuadPoints,5) )
    allocate( limitedPntVal(nQuadPoints,5) )
    allocate( hatDens(nQuadPoints))
    allocate( t_vec(nQuadPoints ))


    elemLoop: do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! get the modal coefficients of the current cell (for all variables
      ! of the Euler equation, therefore we use ":" for the third index).
      ! ATTENTION: have to be duplicated as the FPT is modifying the input
      ! vector.

      ! Read out the mean values (before applying the limiter)
      mean(:) = state%state(iElem,1,:)
      pressMean = (euler%isen_coef-1.0_rk) * ( mean(5) &
        & - (0.5_rk/mean(1)) * ( mean(2)*mean(2)       &
        & + mean(3)*mean(3) + mean(4)*mean(4) ) )

      highmean: if (mean(1) > filter%eps .and. pressMean > filter%eps) then
        ! Only if the mean values are above the threshold the filtering
        ! is reasonable to apply.

        ! --> modal space
        call ply_convert2oversample(state       = state%state(iElem,:,:), &
          &                         poly_proj   = poly_proj,              &
          &                         nDim        = 3,                      &
          &                         modalCoeffs = modalCoeffs             )
        ! --> oversamp modal space

        ! Now, we transform the modal representation of this element to nodal
        ! space by making use of fast polynomial transformations (FPT)
        call ply_poly_project_m2n(me = poly_proj ,       &
           &                      dim = 3 ,              &
           &                      nVars = 5,             &
           &                      nodal_data=pointVal,   &
           &                      modal_data=modalCoeffs )
        ! --> oversamp nodal space

        ! Step 1: Limit the density
        minDens = minval( pointVal(:,1) )
        if ( (mean(1) - minDens) > epsilon(minDens)*mean(1)) then
          theta_1 = min(abs((mean(1) - filter%eps)/(mean(1) - minDens)),1.0_rk)
          do iPoint = 1, nQuadPoints
            hatDens(iPoint) = theta_1*(pointVal(iPoint,1) - mean(1)) + mean(1)
          end do
        else
          hatDens = pointVal(:,1)
        end if

        ! Step 2: limit the pressure
        do iPoint = 1,nQuadPoints
          momSq = sum(pointVal(iPoint,2:4)*pointVal(iPoint,2:4))
          tmpPressure = (euler%isen_coef - 1.0_rk) * (pointVal(iPoint,5) &
            & - (1.0_rk/2.0_rk)*((momSq)/hatDens(iPoint)))
          if (tmpPressure < filter%eps) then
            ! Intersection point of mean value state and admissible set G
            ! (see paper Shu, et al)
            t_vec(iPoint) = solve_admissible_state( euler%isen_coef,     &
              &                                     filter%eps,          &
              &                                     mean,                &
              &                                     hatDens(iPoint),     &
              &                                     pointVal(iPoint,2:5) )
          else
            t_vec(iPoint) = 1.0_rk
          end if
        end do
        theta_2 = minval(t_vec)

        do i=lbound(limitedPntVal, 1),ubound(limitedPntVal,1)
          limitedPntVal(i,1) = theta_2 * ( hatDens(i) - mean(1) ) + mean(1)
        end do

        do iVar = 2,5
          limitedPntVal(:,iVar) = theta_2 * ( pointVal(:,iVar) &
            & - mean(iVar) ) + mean(iVar)
        end do

        ! transform back and write the results to the state
        call ply_poly_project_n2m(me        = poly_proj,     &
          &                      dim        = 3,             &
          &                      nVars      = 5,             &
          &                      nodal_data = limitedPntVal, &
          &                      modal_data = modalCoeffs    )
        ! --> oversamp modal space
        call ply_convertFromOversample(modalCoeffs = modalCoeffs,           &
          &                            poly_proj   = poly_proj,             &
          &                            nDim        = 3,                     &
          &                            state       = state%state(iElem,:,:) )
        ! --> modal space

      else highmean
        ! Density or pressure are so small, that the mean is below the
        ! threshold. This should basically never happen, but if we run into
        ! this we fall back to first order approximation and push the state
        ! above the threshold (not conservative).
        state%state(iElem,:,:) = 0.0_rk
        mean(1) = max(mean(1), filter%eps)
        pressMean = max(pressMean, filter%eps)
        mean(5) = pressMean/(euler%isen_coef-1.0_rk) &
          &     + 0.5_rk/mean(1) * (  mean(2)**2     &
          &                         + mean(3)**2     &
          &                         + mean(4)**2     )
        state%state(iElem,1,:) = mean
      end if highmean

    end do elemLoop


  end subroutine atl_cons_positivity_preserv

  !> Solve for admissible state of the conservative positivity preserving
  !! limiter.
  function solve_admissible_state( isen_coef, my_eps, mean, limitedDens, &
      &                            momEnerg ) result(t)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: isen_coef
    real(kind=rk), intent(in) :: my_eps
    real(kind=rk), intent(in) :: mean(5)
    real(kind=rk), intent(in) :: limitedDens
    real(kind=rk), intent(in) :: momEnerg(2:5)
    real(kind=rk) :: t
    ! --------------------------------------------------------------------------
    real(kind=rk) :: p_dash, p, q, t_1, t_2, tmp
    ! --------------------------------------------------------------------------

    ! we solve p_dash*t^2 + p*t + q = 0.
    ! definte the constants of the polynomial

    p_dash = (isen_coef-1.0_rk) * (mean(5) - momEnerg(5)) * (mean(1) &
      & - limitedDens) - ((isen_coef-1.0_rk)/2.0_rk) * ((mean(2)     &
      & - momEnerg(2))**2 + (mean(3) - momEnerg(3))**2 + (mean(4)    &
      & - momEnerg(4))**2)

    p = mean(1)*my_eps - limitedDens * my_eps                            &
      & + (isen_coef-1.0_rk) * (mean(5) * limitedDens - 2.0_rk * mean(5) &
      & * mean(1) + momEnerg(5) * mean(1))                               &
      & - ((isen_coef-1.0_rk)/2.0_rk)                                    &
      &   * (   2.0_rk*mean(2)*(momEnerg(2)-mean(2))                     &
      &       + 2.0_rk*mean(3)*(momEnerg(3)-mean(3))                     &
      &       + 2.0_rk*mean(4)*(momEnerg(4)-mean(4)) )

    q = (-my_eps) * mean(1) + (isen_coef-1.0_rk) * mean(1) * mean(5)  &
      & - ((isen_coef-1.0_rk)/2.0_rk) * ( mean(2) * mean(2) + mean(3) &
      & * mean(3) + mean(4) * mean(4) )

    ! Now, solve the quadratic equation by using the p-q-formula
    p = p / p_dash
    q = q / p_dash

    if( ((p/2.0_rk)**2 - q) < tiny(p) )then
      write(*,*) 'p/2.0_rk)**2-q =',  (p/2.0_rk)**2 - q
    end if
    tmp = max( (p/2.0_rk)**2 - q, tiny(p) )

    t_1 = - (p/2.0_rk) + sqrt( tmp )
    t_2 = - (p/2.0_rk) - sqrt( tmp )

    ! Identify the valid t
    if( t_1 <= 1.0_rk .and. t_1 >= 0.0_rk ) then
      t = t_1
    elseif(t_2 <= 1.0_rk .and. t_2 >= 0.0_rk ) then
      t = t_2;
    else
      call tem_abort( 'Error in solve_admissible_state:' &
        & // ' t is not in interval [0;1]'               )
    end if

  end function solve_admissible_state

  !> Apply conservative limitation of denisty and energy to ensure positivity
  !! for density and pressure.
  subroutine atl_cons_positivity_preserv_2d( state, mesh, filter, euler, &
    &                                        poly_proj )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_cons_positivity_preserv_type), intent(in) :: filter
    type(atl_Euler_type), intent(in) :: euler
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:), &
      &                           pointVal(:,:),    &
      &                           hatDens(:),       &
      &                           t_vec(:),         &
      &                           limitedPntVal(:,:)
    real(kind=rk) :: mean(4)
    ! Loop vars
    integer :: iElem, iPoint, iVar, j
    real(kind=rk) :: minDens, theta_1, theta_2, momSq, tmpPressure
    real(kind=rk) :: my_eps, pressMean
    integer :: nquadpoints
    integer :: oversamp_dofs
    ! --------------------------------------------------------------------------

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs= poly_proj%body_2d%oversamp_dofs

    allocate( modalCoeffs(oversamp_dofs,4) )
    allocate( pointVal(nQuadPoints,4) )
    allocate( limitedPntVal(nQuadPoints,4) )
    allocate( hatDens(nQuadPoints) )
    allocate( t_vec(nQuadPoints) )


    elemLoop: do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! get the modal coefficients of the current cell (for all variables
      ! of the Euler equation, therefore we use ":" for the third index).
      ! ATTENTION: have to be duplicated as the FPT is modifying the input
      !            vector.

      ! Read out the mean values (before applying the limiter)
      mean(:) = state%state(iElem,1,:)
      pressMean = (euler%isen_coef-1.0_rk) * ( mean(4) &
        & - (0.5_rk/mean(1)) * (mean(2) * mean(2) + mean(3) * mean(3)) )

      highmean: if (mean(1) > filter%eps .and. pressMean > filter%eps) then
        ! Only if the mean values are above the threshold the filtering
        ! is reasonable to apply.

        ! --> modal oversample
        call ply_convert2oversample(state= state%state(iElem,:,:), &
          &                         poly_proj= poly_proj,          &
          &                         nDim = 2,                      &
          &                         modalCoeffs = modalCoeffs      )
        ! --> oversample modal oversample

        ! Smallness parameter identification is missing ...
        my_eps = min( filter%eps, mean(1), pressMean )

        ! Now, we transform the modal representation of this element to nodal
        ! space by making use of fast polynomial transformations (FPT)
        call ply_poly_project_m2n(me = poly_proj ,      &
          &                      dim = 2 ,              &
          &                      nVars = 4,             &
          &                      nodal_data=pointVal,   &
          &                      modal_data=modalCoeffs )
        ! --> oversamp nodal space

        ! Step 1: Limit the density
        minDens = minval( pointVal(:,1) )
        if ( (mean(1) - minDens) > epsilon(minDens)*mean(1)) then
          theta_1 = min(abs((mean(1) - my_eps)/(mean(1) - minDens)),1.0_rk)
          do iPoint = 1, nQuadPoints
            hatDens(iPoint) = theta_1*(pointVal(iPoint,1) - mean(1)) + mean(1)
          end do
        else
          hatDens = pointVal(:,1)
        end if

        ! Step 2: limit the pressure
        do iPoint = 1, nQuadPoints
          momSq = sum(pointVal(iPoint,2:3)*pointVal(iPoint,2:3))
          tmpPressure = (euler%isen_coef - 1.0_rk) * (pointVal(iPoint,4) &
            & - (1.0_rk / 2.0_rk) * ((momSq) / hatDens(iPoint)))
          if(tmpPressure < my_eps) then
            ! Intersection point of mean value state and admissible set G (see
            ! paper Shu, et al)
            t_vec(iPoint) = solve_admissible_state_2d( euler%isen_coef,     &
              &                                        my_eps,              &
              &                                        mean,                &
              &                                        hatDens(iPoint),     &
              &                                        pointVal(iPoint,2:4) )
          else
            t_vec(iPoint) = 1.0_rk
          end if
        end do
        theta_2 = minval(t_vec)

        do j=lbound(limitedPntVal,1),ubound(limitedPntVal,1)
          limitedPntVal(j,1) = theta_2 * ( hatDens(j) - mean(1) ) + mean(1)
        end do



        do iVar = 2,4
          limitedPntVal(:,iVar) = theta_2 * ( pointVal(:,iVar) - mean(iVar) ) &
            & + mean(iVar)
        end do

        ! transform back and write the results to the state
        call ply_poly_project_n2m(me        = poly_proj,     &
          &                      dim        = 2,             &
          &                      nVars      = 4,             &
          &                      nodal_data = limitedpntval, &
          &                      modal_data = modalcoeffs    )
        ! --> oversamp modal space
        ! --> oversamp modal space
        call ply_convertFromOversample(modalCoeffs = modalCoeffs,           &
          &                            poly_proj   = poly_proj,             &
          &                            nDim        = 2,                     &
          &                            state       = state%state(iElem,:,:) )
        ! --> modal space

      else highmean
        ! Density or pressure are so small, that the mean is below the
        ! threshold. This should basically never happen, but if we run into
        ! this we fall back to first order approximation and push the state
        ! above the threshold (not conservative).
        state%state(iElem,:,:) = 0.0_rk
        mean(1) = max(mean(1), filter%eps)
        pressMean = max(pressMean, filter%eps)
        mean(4) = pressMean/(euler%isen_coef-1.0_rk) &
          &     + 0.5_rk/mean(1) * (  mean(2)**2     &
          &                         + mean(3)**2     )
        state%state(iElem,1,:) = mean
      end if highmean

    end do elemLoop


  end subroutine atl_cons_positivity_preserv_2d

  !> Solve for admissible state of the conservative positivity preserving limiter.
  function solve_admissible_state_2d( isen_coef, my_eps, mean, limitedDens, &
      &                               momEnerg ) result(t)
    ! --------------------------------------------------------------------------
    real(kind=rk), intent(in) :: isen_coef
    real(kind=rk), intent(in) :: my_eps
    real(kind=rk), intent(in) :: mean(4)
    real(kind=rk), intent(in) :: limitedDens
    real(kind=rk), intent(in) :: momEnerg(2:4)
    real(kind=rk) :: t
    ! --------------------------------------------------------------------------
    real(kind=rk) :: p_dash, p, q, t_1, t_2, tmp
    ! --------------------------------------------------------------------------

    ! we solve p_dash*t^2 + p*t + q = 0.
    ! definte the constants of the polynomial

    p_dash = (isen_coef-1.0_rk) * (mean(4) - momEnerg(4)) * (mean(1) &
      & - limitedDens) - ((isen_coef - 1.0_rk) / 2.0_rk) &
      & * ((mean(2) - momEnerg(2))**2 + (mean(3) - momEnerg(3))**2)

    p = mean(1)*my_eps - limitedDens * my_eps                               &
      & + (isen_coef-1.0_rk) * (mean(4) * limitedDens - 2.0_rk * mean(4)    &
      & * mean(1) + momEnerg(4) * mean(1))                                  &
      & - ((isen_coef - 1.0_rk) / 2.0_rk)                                   &
      &   * ( 2.0_rk * mean(2) * (momEnerg(2) - mean(2)) + 2.0_rk * mean(3) &
      & * (momEnerg(3) - mean(3)) )

    q = (-my_eps) * mean(1) + (isen_coef - 1.0_rk) * mean(1) * mean(4)    &
      & - ((isen_coef - 1.0_rk) / 2.0_rk) * ( mean(2) * mean(2) + mean(3) &
      & * mean(3) )

    ! Now, solve the quadratic equation by using the p-q-formula
    p = p / p_dash
    q = q / p_dash
    tmp = max( (p/2.0_rk)**2 - q , tiny(p))

    t_1 = - (p/2.0_rk) + sqrt( tmp )
    t_2 = - (p/2.0_rk) - sqrt( tmp )

    ! Identify the valid t
    if( t_1 <= 1.0_rk .and. t_1 >= 0.0_rk ) then
      t = t_1
    elseif(t_2 <= 1.0_rk .and. t_2 >= 0.0_rk ) then
      t = t_2
    else
      call tem_abort( 'Error in solve_admissible_state_2d:' &
        & // ' t is not in interval [0;1]'                  )
    end if

  end function solve_admissible_state_2d


  !> Apply pointwise limitation of denisty and energy to ensure positivity
  !! for density and pressure.
  subroutine atl_positivity_preserv(state, mesh, filter, euler, poly_proj)
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_positivity_preserv_type), intent(in) :: filter
    type(atl_Euler_type), intent(in) :: euler
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:), pointVal(:,:)
    real(kind=rk) :: limitedPntVal(5)
    real(kind=rk) :: kinEnergy
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nquadpoints, ndofs
    integer :: oversamp_dofs
    ! --------------------------------------------------------------------------

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_3D%nquadpoints
    ndofs = poly_proj%body_3D%ndofs
    oversamp_dofs= poly_proj%body_3d%oversamp_dofs

    allocate( modalCoeffs(oversamp_dofs,5) )
    allocate( pointVal(nQuadPoints,5) )


    elemLoop: do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! get the modal coefficients of the current cell (for all variables
      ! of the Euler equation, therefore we use ":" for the third index).
      ! ATTENTION: have to be duplicated as the FPT is modifying the input
      ! vector.

      ! --> modal oversample
      call ply_convert2oversample(state       = state%state(iElem,:,:), &
        &                         poly_proj   = poly_proj,              &
        &                         nDim        = 3,                      &
        &                         modalCoeffs = modalCoeffs             )
      ! --> oversample modal oversample

      ! Now, we transform the modal representation of this element to nodal
      ! space by making use of fast polynomial transformations (FPT)
      call ply_poly_project_m2n(me         = poly_proj,  &
         &                      dim        = 3,          &
         &                      nVars      = 5,          &
         &                      nodal_data = pointVal,   &
         &                      modal_data = modalCoeffs )
      ! --> oversampled nodal space

      do iPoint = 1, nQuadPoints

        limitedPntVal(:) = pointVal(iPoint,:)

        ! Check if denisty is below threshold at this point
        if (limitedPntVal(1) < filter%eps) then
          limitedPntVal(1) = filter%eps*(1.0_rk+epsilon(filter%eps))
        end if

        ! Check if pressure is below threshold at this point
        kinEnergy = 0.5_rk * (                                     &
          &                      limitedPntVal(2)*limitedPntVal(2) &
          &                    + limitedPntVal(3)*limitedPntVal(3) &
          &                    + limitedPntVal(4)*limitedPntVal(4) &
          &                    ) / limitedPntVal(1)
        if ( limitedPntVal(5)                                     &
          & < kinEnergy + filter%eps / (euler%isen_coef - 1.0_rk) ) then
          limitedPntVal(5) = kinEnergy + filter%eps / (euler%isen_coef - 1.0_rk)
        end if

        ! Write back the limited values
        pointVal(iPoint,:) = limitedPntVal(:)

      end do


      ! transform back and write the results to the state
      call ply_poly_project_n2m(me         = poly_proj,  &
         &                      dim        = 3,          &
         &                      nVars      = 5,          &
         &                      nodal_data = pointVal,   &
         &                      modal_data = modalCoeffs )
      ! --> oversamp modal space
      call ply_convertFromOversample(modalCoeffs = modalCoeffs,           &
        &                            poly_proj   = poly_proj,             &
        &                            nDim        = 3,                     &
        &                            state       = state%state(iElem,:,:) )
      ! --> modal space

    end do elemLoop


  end subroutine atl_positivity_preserv

  !> Apply pointwise limitation of denisty and energy to ensure positivity
  !! for density and pressure.
  subroutine atl_positivity_preserv_2d( state, mesh, filter, euler, poly_proj )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_positivity_preserv_type), intent(in) :: filter
    type(atl_Euler_type), intent(in) :: euler
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    ! The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:), pointVal(:,:)
    real(kind=rk) :: limitedPntVal(4)
    real(kind=rk) :: kinEnergy
    ! Loop vars
    integer :: iElem, iPoint
    integer :: nQuadPoints
    integer :: oversamp_dofs
    ! --------------------------------------------------------------------------

    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nQuadPoints = poly_proj%body_2D%nquadpoints
    oversamp_dofs = poly_proj%body_2D%oversamp_dofs

    allocate( modalCoeffs(oversamp_dofs,4) )
    allocate( pointVal(nQuadPoints,4) )


    elemLoop: do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

      ! get the modal coefficients of the current cell (for all variables
      ! of the Euler equation, therefore we use ":" for the third index).
      ! ATTENTION: have to be duplicated as the FPT is modifying the input
      ! vector.
      ! --> oversamp modal space
      call ply_convert2oversample(state       = state%state(iElem,:,:), &
        &                         poly_proj   = poly_proj,              &
        &                         nDim        = 2,                      &
        &                         modalCoeffs = modalCoeffs             )
      ! --> oversamp modal space

      ! Now, we transform the modal representation of this element to nodal
      ! space by making use of fast polynomial transformations (FPT)
      call ply_poly_project_m2n(me = poly_proj ,       &
         &                      dim = 2 ,              &
         &                      nVars = 4,             &
         &                      nodal_data=pointVal,   &
         &                      modal_data=modalCoeffs )
      ! --> oversamp nodal space

      do iPoint = 1, nQuadPoints

        limitedPntVal(:) = pointVal(iPoint,:)

        ! Check if denisty is below threshold at this point
        if(limitedPntVal(1) < filter%eps) then
          limitedPntVal(1) = filter%eps
        end if

        ! Check if pressure is below threshold at this point
        kinEnergy = 0.5_rk * (limitedPntVal(2) * limitedPntVal(2) &
          & + limitedPntVal(3) * limitedPntVal(3)) / limitedPntVal(1)
        if( limitedPntVal(4) &
          & < kinEnergy + filter%eps / (euler%isen_coef - 1.0_rk) ) then
          limitedPntVal(4) =kinEnergy + filter%eps / (euler%isen_coef - 1.0_rk)
        end if

        ! Write back the limited values
        pointVal(iPoint,:) = limitedPntVal(:)

      end do


      ! transform back and write the results to the state
      call ply_poly_project_n2m(me = poly_proj ,       &
         &                      dim = 2 ,              &
         &                      nVars = 4,             &
         &                      nodal_data=pointVal,   &
         &                      modal_data=modalCoeffs )
      ! --> oversamp modal space
      call ply_convertFromOversample(modalCoeffs = modalCoeffs,           &
        &                            poly_proj   = poly_proj,             &
        &                            nDim        = 2,                     &
        &                            state       = state%state(iElem,:,:) )
      ! --> modal space

    end do elemLoop


  end subroutine atl_positivity_preserv_2d


  !> Covolume filtering for 3D equation.
  subroutine atl_covolume( minlevel, maxlevel, state, state_stab, mesh,        &
    & filter, scheme, equation, tree, poly_proj, poly_proj_pos, bc, boundary,  &
    & general, commStateTimer, adaptive_orders                                 )
    ! --------------------------------------------------------------------------
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> State to be filtered (input and output)
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> The actual co-volume filter to be applied.
    type(atl_covolume_type), intent(in) :: filter
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    !> The filters orders, if adaptive filter is applied
    type(atl_adaptive_orders_type), intent(in), optional :: adaptive_orders(minlevel:maxlevel)
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! We compute the co-volume filter in a dim-by-dim procedure. Hence
    ! we loop over all the directions and use the results of the previous
    ! direction as an input for the next direction.
    do iDir = 1,3

       ! Project from primal mesh to covolume mesh
       call atl_covolume_tocovolume( minlevel        = minlevel,        &
         &                           maxlevel        = maxLevel,        &
         &                           iDir            = iDir,            &
         &                           state           = state,           &
         &                           state_stab      = state_stab,      &
         &                           mesh            = mesh,            &
         &                           filter          = filter,          &
         &                           scheme          = scheme,          &
         &                           equation        = equation,        &
         &                           tree            = tree,            &
         &                           poly_proj       = poly_proj,       &
         &                           poly_proj_pos   = poly_proj_pos,   &
         &                           bc              = bc,              &
         &                           boundary        = boundary,        &
         &                           general         = general,         &
         &                           commStateTimer  = commStateTimer,  &
         &                           adaptive_orders = adaptive_orders  )

     end do

  end subroutine atl_covolume


  !> Recursive routine to project the state from primal mesh to covolume mesh.
  !!
  !! This routine is a recursive subroutine to project from the primal mesh
  !! to the covolume mesh. It interpolates to/from the finer level
  !! (if available) and calls itself on the next finer level.
  !! The projection is carried out by a L2 projection.
  subroutine atl_covolume_tocovolume( minlevel, maxlevel, iDir, state,         &
    & state_stab, mesh, filter, scheme, equation, tree, poly_proj,             &
    & poly_proj_pos, bc, boundary, general, commStateTimer, adaptive_orders    )
    ! --------------------------------------------------------------------------
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> The spatial direction for the projection:
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction \n
    integer , intent(in) :: iDir
    !> State to be filtered
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> The actual co-volume filter to be applied.
    type(atl_covolume_type), intent(in) :: filter
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    !> The filters orders, if adaptive filter is applied
    type(atl_adaptive_orders_type), intent(in), optional :: adaptive_orders(minlevel:maxlevel)
    ! --------------------------------------------------------------------------
    integer :: maxPolyDeg
    integer :: left_neighbor, right_neighbor
    integer :: iside, iLevel
    integer :: nElems_fluid, nDofs
    real(kind=rk), allocatable :: leftCovolume(:,:), rightCovolume(:,:), &
                                & prev_state(:,:), modalCoeff(:,:), &
                                & leftModalCoeff(:,:), rightModalCoeff(:,:)
    ! --------------------------------------------------------------------------

    ! Interpolate element states (output is in state_stab)
    call atl_interpolate_elemstate( minLevel       = minLevel,      &
      &                             maxLevel       = maxLevel,      &
      &                             iLevel         = maxLevel,      &
      &                             iDir           = iDir,          &
      &                             state          = state,         &
      &                             state_stab     = state_stab,    &
      &                             mesh           = mesh,          &
      &                             scheme         = scheme,        &
      &                             equation       = equation,      &
      &                             tree           = tree,          &
      &                             dimen          = 3,             &
      &                             poly_proj      = poly_proj,     &
      &                             poly_proj_pos  = poly_proj_pos, &
      &                             bc             = bc,            &
      &                             boundary       = boundary,      &
      &                             general        = general,       &
      &                             commStateTimer = commStateTimer )

    do iLevel = minLevel, maxLevel

      maxPolyDeg = scheme(iLevel)%modg%maxPolyDegree
      nDofs = (maxPolyDeg+1)**3

      allocate( leftCovolume( nDofs, equation%varSys%nScalars) )
      allocate( rightCovolume( nDofs, equation%varSys%nScalars) )
      allocate( prev_state( nDofs, equation%varSys%nScalars) )
      allocate( modalCoeff( nDofs, equation%varSys%nScalars) )
      allocate( leftModalCoeff( nDofs, equation%varSys%nScalars) )
      allocate( rightModalCoeff( nDofs, equation%varSys%nScalars) )

      ! Project onto the co-volume
      nElems_fluid = mesh(iLevel)%faces_stab         &
        &                        %dimByDimDesc(iDir) &
        &                        %elem               &
        &                        %nElems(eT_fluid)
      do iside = 1, nElems_fluid

        ! Get the neighbors for this face.
        left_neighbor = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%neigh(1)%nghElems(1,iside)
        right_neighbor = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%neigh(1)%nghElems(2,iside)


        ! Get modal coefficients of current element and left and right neighbor
        modalCoeff = state_stab(iLevel,iDir)%state(iside,:,:)
        leftModalCoeff = state_stab(iLevel,iDir)%state(left_neighbor,:,:)
        rightModalCoeff = state_stab(iLevel,iDir)%state(right_neighbor,:,:)


        ! Calculate projection onto co-volume (left and right)
        if(present(adaptive_orders)) then
          call atl_primal_to_covolume_projection(                 &
            & left       = leftModalCoeff,                        &
            & right      = modalCoeff,                            &
            & dir        = iDir,                                  &
            & filter     = filter,                                &
            & scheme     = scheme(iLevel),                        &
            & order      = adaptive_orders(iLevel)%orders(iside), &
            & maxPolyDeg = maxPolyDeg,                            &
            & covolume   = leftCovolume                           )
          call atl_primal_to_covolume_projection(                 &
            & left       = modalCoeff,                            &
            & right      = rightModalCoeff,                       &
            & dir        = iDir,                                  &
            & filter     = filter,                                &
            & scheme     = scheme(iLevel),                        &
            & order      = adaptive_orders(iLevel)%orders(iside), &
            & maxPolyDeg = maxPolyDeg,                            &
            & covolume   = rightCovolume                          )
        else
          call atl_primal_to_covolume_projection( &
            & left       = leftModalCoeff,        &
            & right      = modalCoeff,            &
            & dir        = iDir,                  &
            & filter     = filter,                &
            & scheme     = scheme(iLevel),        &
            & order      = filter%order,          &
            & maxPolyDeg = maxPolyDeg,            &
            & covolume   = leftCovolume           )
          call atl_primal_to_covolume_projection( &
            & left       = modalCoeff,            &
            & right      = rightModalCoeff,       &
            & dir        = iDir,                  &
            & filter     = filter,                &
            & scheme     = scheme(iLevel),        &
            & order      = filter%order,          &
            & maxPolyDeg = maxPolyDeg,            &
            & covolume   = rightCovolume          )
        end if

        ! Project back to the primal mesh.
        prev_state(:,:) = state(iLevel)%state(iside,:,:)
        state(iLevel)%state(iside,:,:) = atl_covolume_to_primal_projection( &
          & left       = leftCovolume,                                      &
          & right      = rightCovolume,                                     &
          & dir        = iDir,                                              &
          & filter     = filter,                                            &
          & scheme     = scheme(iLevel),                                    &
          & nScalars   = equation%varSys%nScalars,                          &
          & state      = prev_state,                                        &
          & maxPolyDeg = maxPolyDeg                                         )

      end do

      deallocate( prev_state )
      deallocate( leftCovolume )
      deallocate( rightCovolume )
      deallocate( modalCoeff )
      deallocate( leftModalCoeff )
      deallocate( rightModalCoeff )

    end do

  end subroutine atl_covolume_tocovolume


  !> Recursive interpolation of element states among the levels of the
  !! mesh.
  recursive subroutine atl_interpolate_elemstate(minlevel, maxlevel, iLevel,   &
    & iDir, state, state_stab, mesh, scheme, equation, tree, dimen,poly_proj,  &
    & poly_proj_pos, bc, boundary, general, commStateTimer                     )
    ! --------------------------------------------------------------------------!
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> The current level.
    integer, intent(in) :: iLevel
    !> The direction to interpolate
    integer, intent(in) :: iDir
    !> State to be filtered (input and output)
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The spatial dimension under consideration.
    integer, intent(in) :: dimen
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    ! --------------------------------------------------------------------------!
    integer :: nElems_fluid
    ! --------------------------------------------------------------------------!

    ! Assign state information to temporary array
    nElems_fluid = mesh(iLevel)%faces_stab         &
      &                        %dimByDimDesc(iDir) &
      &                        %elem               &
      &                        %nElems(eT_fluid)
    state_stab(iLevel,iDir)%state(:nElems_fluid,:,:) = &
                                        & state(iLevel)%state(:nElems_fluid,:,:)


    ! Assign boundary values for state_stab
    call atl_set_covolume_bnd(                                     &
      &    bc        = bc,                                         &
      &    boundary  = boundary(iLevel),                           &
      &    state     = state_stab(iLevel,iDir),                    &
      &    scheme    = scheme(iLevel),                             &
      &    poly_proj = poly_proj(boundary(iLevel)%poly_proj_pos),  &
      &    equation  = equation,                                   &
      &    tree      = tree,                                       &
      &    time      = state(iLevel)%local_time,                   &
      &    mesh      = mesh(iLevel),                               &
      &    nodalBnd  = equation%isNonlinear,                       &
      &    iDir      = iDir,                                       &
      &    iLevel    = iLevel                                      )

    ! Communicate co-volume states (on this level)
    call tem_startTimer( timerHandle = commStateTimer )
    call general%commPattern%exchange_real( &
    &    send = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%sendbuffer, &
    &    recv = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%recvbuffer, &
    &    state = state_stab(iLevel,iDir)%state, &
    &    message_flag = iLevel, &
    &    comm = general%proc%comm )
    call tem_stopTimer( timerHandle = commStateTimer )

    if(iLevel < maxlevel) then

      ! interpolate states from finer level to current level
      select case(dimen)
      case(3)
        call atl_modg_fineToCoarseElem( minLevel = minLevel, &
          &                             maxLevel = maxLevel, &
          &                             currentLevel = iLevel, &
          &                             iDir = iDir, &
          &                             mesh = mesh, &
          &                             state_stab = state_stab, &
          &                             scheme = scheme, &
          &                             nScalars = equation%varSys%nScalars )
      case(2)
        call atl_modg_fineToCoarseElem_2d( minLevel = minLevel, &
          &                                maxLevel = maxLevel, &
          &                                currentLevel = iLevel, &
          &                                iDir = iDir, &
          &                                mesh = mesh, &
          &                                state_stab = state_stab, &
          &                                scheme = scheme, &
          &                                nScalars = equation%varSys%nScalars )
      case(1)
        call atl_modg_fineToCoarseElem_1d( minLevel = minLevel, &
          &                                maxLevel = maxLevel, &
          &                                currentLevel = iLevel, &
          &                                iDir = iDir, &
          &                                mesh = mesh, &
          &                                state_stab = state_stab, &
          &                                scheme = scheme, &
          &                                nScalars = equation%varSys%nScalars )
      case default
        write(logUnit(1),*) 'ERROR in atl_interpolate_elemstate: '
        write(logUnit(1),*) ' Unknown spatial dimension, stopping ...'
        call tem_abort()
      end select

      ! Communicate co-volume states (only the from finer elements)
      call tem_startTimer( timerHandle = commStateTimer )
      call general%commPattern%exchange_real( &
      &    send = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%sendbufferfromfiner, &
      &    recv = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%recvbufferfromfiner, &
      &    state = state_stab(iLevel,iDir)%state, &
      &    message_flag = iLevel, &
      &    comm = general%proc%comm )
      call tem_stopTimer( timerHandle = commStateTimer )
    end if

    ! If there is a coarser level, then recurse up to the next level
    if(iLevel > minLevel) then

      call atl_interpolate_elemstate( minLevel       = minLevel,      &
        &                             maxLevel       = maxLevel,      &
        &                             iLevel         = iLevel-1,      &
        &                             iDir           = iDir,          &
        &                             state          = state,         &
        &                             state_stab     = state_stab,    &
        &                             mesh           = mesh,          &
        &                             scheme         = scheme,        &
        &                             equation       = equation,      &
        &                             tree           = tree,          &
        &                             dimen          = dimen,         &
        &                             poly_proj      = poly_proj,     &
        &                             poly_proj_pos  = poly_proj_pos, &
        &                             bc             = bc,            &
        &                             boundary       = boundary,      &
        &                             general        = general,       &
        &                             commStateTimer = commStateTimer )
      ! interpolate from the coarser level to the current level
      select case(dimen)
      case(3)
        call atl_modg_coarseToFineElem( minLevel = minLevel, &
          &                             maxLevel = maxLevel, &
          &                             currentLevel = iLevel, &
          &                             iDir = iDir, &
          &                             mesh = mesh, &
          &                             state_stab = state_stab, &
          &                             scheme = scheme, &
          &                             nScalars = equation%varSys%nScalars )
      case(2)
        call atl_modg_coarseToFineElem_2d( minLevel = minLevel, &
          &                                maxLevel = maxLevel, &
          &                                currentLevel = iLevel, &
          &                                iDir = iDir, &
          &                                mesh = mesh, &
          &                                state_stab = state_stab, &
          &                                scheme = scheme, &
          &                                nScalars = equation%varSys%nScalars )
      case(1)
        call atl_modg_coarseToFineElem_1d( minLevel = minLevel, &
          &                                maxLevel = maxLevel, &
          &                                currentLevel = iLevel, &
          &                                iDir = iDir, &
          &                                mesh = mesh, &
          &                                state_stab = state_stab, &
          &                                scheme = scheme, &
          &                                nScalars = equation%varSys%nScalars )
      case default
        write(logUnit(1),*) 'ERROR in atl_interpolate_elemstate: '
        write(logUnit(1),*) ' Unknown spatial dimension, stopping ...'
        call tem_abort()
      end select

      ! Communicate co-volume states (only the from coarser elements)
      call tem_startTimer( timerHandle = commStateTimer )
      call general%commPattern%exchange_real( &
      &    send = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%sendbufferfromcoarser, &
      &    recv = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%recvbufferfromcoarser, &
      &    state = state_stab(iLevel,iDir)%state, &
      &    message_flag = iLevel, &
      &    comm = general%proc%comm )
      call tem_stopTimer( timerHandle = commStateTimer )

    end if


  end subroutine atl_interpolate_elemstate


  !> Covolume filtering for 2D equation.
  subroutine atl_covolume_2d( minlevel, maxlevel, state, state_stab,  &
    & mesh, filter, scheme, equation, tree, poly_proj, poly_proj_pos, &
    & bc, boundary, general, commStateTimer, adaptive_orders          )
    ! --------------------------------------------------------------------------
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> State to be filtered (input and output)
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> The actual co-volume filter to be applied.
    type(atl_covolume_type), intent(in) :: filter
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    !> The filters orders, if adaptive filter is applied
    type(atl_adaptive_orders_type), intent(in), optional :: adaptive_orders(minlevel:maxlevel)
    ! --------------------------------------------------------------------------
    integer :: iDir
    ! --------------------------------------------------------------------------

    ! We compute the co-volume filter in a dim-by-dim procedure. Hence
    ! we loop over all the directions and use the results of the previous
    ! direction as an input for the next direction.
    do iDir = 1,2

       ! Project from primal mesh to covolume mesh
       call atl_covolume_tocovolume_2d( minlevel        = minlevel,       &
         &                              maxlevel        = maxLevel,       &
         &                              iDir            = iDir,           &
         &                              state           = state,          &
         &                              state_stab      = state_stab,     &
         &                              mesh            = mesh,           &
         &                              filter          = filter,         &
         &                              scheme          = scheme,         &
         &                              equation        = equation,       &
         &                              tree            = tree,           &
         &                              poly_proj       = poly_proj,      &
         &                              poly_proj_pos   = poly_proj_pos,  &
         &                              bc              = bc,             &
         &                              boundary        = boundary,       &
         &                              general         = general,        &
         &                              commStateTimer  = commStateTimer, &
         &                              adaptive_orders = adaptive_orders )

     end do

  end subroutine atl_covolume_2d


  !> Recursive routine to project the state from primal mesh to covolume mesh.
  !!
  !! This routine is a recursive subroutine to project from the primal mesh
  !! to the covolume mesh. It interpolates to/from the finer level
  !! (if available) and calls itself on the next finer level.
  !! The projection is carried out by a L2 projection.
  subroutine atl_covolume_tocovolume_2d( minlevel, maxlevel, iDir, state,   &
    & state_stab, mesh, filter, scheme, equation, tree, poly_proj,          &
    & poly_proj_pos, bc, boundary, general, commStateTimer, adaptive_orders )
    ! ---------------------------------------------------------------------------
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> The spatial direction for the projection:
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction \n
    integer , intent(in) :: iDir
    !> State to be filtered
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> The actual co-volume filter to be applied.
    type(atl_covolume_type), intent(in) :: filter
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    !> The filters orders, if adaptive filter is applied
    type(atl_adaptive_orders_type), intent(in), optional :: adaptive_orders(minlevel:maxlevel)
    ! --------------------------------------------------------------------------
    integer :: maxPolyDeg
    integer :: left_neighbor, right_neighbor
    integer :: iside, iLevel
    integer :: nElems_fluid
    real(kind=rk), allocatable :: leftCovolume(:,:), rightCovolume(:,:), &
                                & prev_state(:,:), modalCoeff(:,:), &
                                & leftModalCoeff(:,:), rightModalCoeff(:,:)
    ! --------------------------------------------------------------------------

    ! Interpolate element states (output is in state_stab)
    call atl_interpolate_elemstate( minLevel = minLevel,            &
                                  & maxLevel = maxLevel,            &
                                  & iLevel = maxLevel,              &
                                  & iDir = iDir,                    &
                                  & state = state,                  &
                                  & state_stab = state_stab,        &
                                  & mesh = mesh,                    &
                                  & scheme = scheme,                &
                                  & equation = equation,            &
                                  & tree = tree,                    &
                                  & dimen = 2,                      &
                                  & poly_proj = poly_proj,          &
                                  & poly_proj_pos = poly_proj_pos,  &
                                  & bc = bc,                        &
                                  & boundary = boundary,            &
                                  & general = general,              &
                                  & commStateTimer = commStateTimer )

    do iLevel = minLevel, maxLevel

      allocate( leftCovolume(scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( rightCovolume(scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( prev_state( scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( modalCoeff( scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( leftModalCoeff( scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( rightModalCoeff(scheme(iLevel)%nDofs, equation%varSys%nScalars) )

      maxPolyDeg = scheme(iLevel)%modg_2d%maxPolyDegree

      ! Project onto the co-volume
      nElems_fluid = mesh(iLevel)%faces_stab         &
        &                        %dimByDimDesc(iDir) &
        &                        %elem               &
        &                        %nElems(eT_fluid)
      do iside = 1, nElems_fluid
        ! Get the neighbors for this face.
        left_neighbor = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%neigh(1)%nghElems(1,iside)
        right_neighbor = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%neigh(1)%nghElems(2,iside)

        ! Get modal coefficients of current element and left and right neighbor
        modalCoeff = state_stab(iLevel,iDir)%state(iside,:,:)
        leftModalCoeff = state_stab(iLevel,iDir)%state(left_neighbor,:,:)
        rightModalCoeff = state_stab(iLevel,iDir)%state(right_neighbor,:,:)

        ! Calculate projection onto co-volume (left and right)
        if(present(adaptive_orders)) then
          call atl_primal_to_covolume_projection_2d(              &
            & left       = leftModalCoeff,                        &
            & right      = modalCoeff,                            &
            & dir        = iDir,                                  &
            & filter     = filter,                                &
            & scheme     = scheme(iLevel),                        &
            & order      = adaptive_orders(iLevel)%orders(iside), &
            & maxPolyDeg = maxPolyDeg,                            &
            & covolume   = leftCovolume                           )
          call atl_primal_to_covolume_projection_2d(              &
            & left       = modalCoeff,                            &
            & right      = rightModalCoeff,                       &
            & dir        = iDir,                                  &
            & filter     = filter,                                &
            & scheme     = scheme(iLevel),                        &
            & order      = adaptive_orders(iLevel)%orders(iside), &
            & maxPolyDeg = maxPolyDeg,                            &
            & covolume   = rightCovolume                          )
        else
          call atl_primal_to_covolume_projection_2d( &
            & left       = leftModalCoeff,           &
            & right      = modalCoeff,               &
            & dir        = iDir,                     &
            & filter     = filter,                   &
            & scheme     = scheme(iLevel),           &
            & order      = filter%order,             &
            & maxPolyDeg = maxPolyDeg,               &
            & covolume   = leftCovolume              )
          call atl_primal_to_covolume_projection_2d( &
            & left       = modalCoeff,               &
            & right      = rightModalCoeff,          &
            & dir        = iDir,                     &
            & filter     = filter,                   &
            & scheme     = scheme(iLevel),           &
            & order      = filter%order,             &
            & maxPolyDeg = maxPolyDeg,               &
            & covolume   = rightCovolume             )
        end if

        ! Project back to the primal mesh.
        prev_state(:,:) = state(iLevel)%state(iside,:,:)
        state(iLevel)%state(iside,:,:) = atl_covolume_to_primal_projection_2d( &
          & left       = leftCovolume,                                         &
          & right      = rightCovolume,                                        &
          & dir        = iDir,                                                 &
          & filter     = filter,                                               &
          & scheme     = scheme(iLevel),                                       &
          & nScalars   = equation%varSys%nScalars,                             &
          & state      = prev_state,                                           &
          & maxPolyDeg = maxPolyDeg                                            )

      end do

      deallocate( prev_state )
      deallocate( leftCovolume )
      deallocate( rightCovolume )
      deallocate( modalCoeff )
      deallocate( leftModalCoeff )
      deallocate( rightModalCoeff )

    end do

  end subroutine atl_covolume_tocovolume_2d

  !> Covolume filtering for 1D equation.
  subroutine atl_covolume_1d( minlevel, maxlevel, state, state_stab, mesh,     &
    & filter, scheme, equation, tree, poly_proj, poly_proj_pos, bc, boundary,  &
    & general, commStateTimer )
    ! --------------------------------------------------------------------------
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> State to be filtered (input and output)
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> The actual co-volume filter to be applied.
    type(atl_covolume_type), intent(in) :: filter
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Project from primal mesh to covolume mesh
    call atl_covolume_tocovolume_1d( minlevel = minlevel, maxlevel = maxLevel, &
                                & iDir = 1, &
                                & state = state, state_stab = state_stab, &
                                & mesh = mesh, filter = filter, &
                                & scheme = scheme, equation = equation, &
                                & tree = tree, &
                                & poly_proj = poly_proj, &
                                & poly_proj_pos = poly_proj_pos, &
                                & bc = bc, boundary = boundary, &
                                & general = general, &
                                & commStateTimer = commStateTimer  )


  end subroutine atl_covolume_1d

  !> Recursive routine to project the state from primal mesh to covolume mesh.
  !!
  !! This routine is a recursive subroutine to project from the primal mesh
  !! to the covolume mesh. It interpolates to/from the finer level
  !! (if available) and calls itself on the next finer level.
  !! The projection is carried out by a L2 projection.
  subroutine atl_covolume_tocovolume_1d( minlevel, maxlevel, iDir, state,      &
    & state_stab, mesh, filter, scheme, equation, tree, poly_proj,             &
    & poly_proj_pos, bc, boundary, general, commStateTimer )
    ! --------------------------------------------------------------------------
    !> The minimal refinement level of the mesh.
    integer, intent(in) :: minlevel
    !> The maximal refinement level of the mesh.
    integer, intent(in) :: maxlevel
    !> The spatial direction for the projection:
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction \n
    integer , intent(in) :: iDir
    !> State to be filtered
    type(atl_statedata_type), intent(inout) :: state(minlevel:maxlevel)
    type(atl_statedata_type), intent(inout) :: state_stab(minlevel:maxlevel,1:3)
    !> Mesh information for all the levels.
    type(atl_cube_elem_type), intent(inout) :: mesh(minlevel:maxlevel)
    !> The actual co-volume filter to be applied.
    type(atl_covolume_type), intent(in) :: filter
    !> List of numerical schemes for all the levels.
    type(atl_scheme_type), intent(inout) :: scheme(minlevel:maxlevel)
    !> Equation kind information
    type(atl_equations_type), intent(in) :: equation
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> The list of projections.
    type(ply_poly_project_type), intent(inout) :: poly_proj(:)
    !> The mapping from each level to the projections.
    integer, intent(inout) :: poly_proj_pos(minlevel:maxlevel)
    !> The global description of the boundaries.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> The levelwise collection of boundary elements and boundary faces.
    type(atl_level_boundary_type), intent(in) :: boundary(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    !> Timer for measuring the communication time inside this routine.
    integer,intent(inout) :: commStateTimer
    ! --------------------------------------------------------------------------
    integer :: maxPolyDeg
    integer :: left_neighbor, right_neighbor
    integer :: iside, iLevel
    integer :: nElems_fluid
    real(kind=rk), allocatable :: leftCovolume(:,:), rightCovolume(:,:), &
                                & prev_state(:,:), modalCoeff(:,:), &
                                & leftModalCoeff(:,:), rightModalCoeff(:,:)
    ! --------------------------------------------------------------------------

    ! Interpolate element states (output is in state_stab)
    call atl_interpolate_elemstate( minLevel = minLevel,            &
                                  & maxLevel = maxLevel,            &
                                  & iLevel = maxLevel,              &
                                  & iDir = iDir,                    &
                                  & state = state,                  &
                                  & state_stab = state_stab,        &
                                  & mesh = mesh,                    &
                                  & scheme = scheme,                &
                                  & equation = equation,            &
                                  & tree = tree,                    &
                                  & dimen = 1,                      &
                                  & poly_proj = poly_proj,          &
                                  & poly_proj_pos = poly_proj_pos,  &
                                  & bc = bc,                        &
                                  & boundary = boundary,            &
                                  & general = general,              &
                                  & commStateTimer = commStateTimer )

    do iLevel = minLevel, maxLevel

      allocate( leftCovolume(scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( rightCovolume(scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( prev_state( scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( modalCoeff( scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( leftModalCoeff( scheme(iLevel)%nDofs, equation%varSys%nScalars) )
      allocate( rightModalCoeff(scheme(iLevel)%nDofs, equation%varSys%nScalars) )

      maxPolyDeg = scheme(iLevel)%modg_1d%maxPolyDegree

      ! Project onto the co-volume
      nElems_fluid = mesh(iLevel)%faces_stab         &
        &                        %dimByDimDesc(iDir) &
        &                        %elem               &
        &                        %nElems(eT_fluid)
      do iside = 1, nElems_fluid
        ! Get the neighbors for this face.
        left_neighbor = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%neigh(1)%nghElems(1,iside)
        right_neighbor = mesh(iLevel)%faces_stab%dimByDimDesc(iDir)%neigh(1)%nghElems(2,iside)

        ! Get modal coefficients of current element and left and right neighbor
        modalCoeff = state_stab(iLevel,iDir)%state(iside,:,:)
        leftModalCoeff = state_stab(iLevel,iDir)%state(left_neighbor,:,:)
        rightModalCoeff = state_stab(iLevel,iDir)%state(right_neighbor,:,:)

        ! Calculate projection onto co-volume (left and right)
        call atl_primal_to_covolume_projection_1d( &
          & left       = leftModalCoeff,           &
          & right      = modalCoeff,               &
          & filter     = filter,                   &
          & scheme     = scheme(iLevel),           &
          & maxPolyDeg = maxPolyDeg,               &
          & covolume   = leftCovolume              )
        call atl_primal_to_covolume_projection_1d( &
          & left       = modalCoeff,               &
          & right      = rightModalCoeff,          &
          & filter     = filter,                   &
          & scheme     = scheme(iLevel),           &
          & maxPolyDeg = maxPolyDeg,               &
          & covolume   = rightCovolume             )

        ! Project back to the primal mesh.
        prev_state(:,:) = state(iLevel)%state(iside,:,:)
        state(iLevel)%state(iside,:,:) = atl_covolume_to_primal_projection_1d( &
          & left       = leftCovolume,                                         &
          & right      = rightCovolume,                                        &
          & filter     = filter,                                               &
          & scheme     = scheme(iLevel),                                       &
          & nScalars   = equation%varSys%nScalars,                             &
          & state      = prev_state,                                           &
          & maxPolyDeg = maxPolyDeg                                            )

      end do

      deallocate( prev_state )
      deallocate( leftCovolume )
      deallocate( rightCovolume )
      deallocate( modalCoeff )
      deallocate( leftModalCoeff )
      deallocate( rightModalCoeff )

    end do

  end subroutine atl_covolume_tocovolume_1d


  !> Damp )the modal coefficients of the state vector by a given spectral viscosity
  !! method.
  subroutine atl_spectral_visc_1d( state, mesh, filter, maxPolyDeg )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_spectral_visc_type), intent(in) :: filter
    integer, intent(in) :: maxPolyDeg
    ! --------------------------------------------------------------------------
    integer :: iElem, dof
    real(kind=rk) :: dofAbs, damping, cut
    integer :: nElems, mpd1
    integer :: iElemDeg
    ! --------------------------------------------------------------------------

    ! Check for single degree of freedom to avoid divisions by zeros in the filter.
    if(maxPolyDeg.le.0) then
      return
    end if

    if(filter%cut_order .le. 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    nElems = mesh%descriptor%elem%nElems(eT_fluid)
    mpd1 = maxPolyDeg+1

    select case(filter%kind)
    case(atl_exp_spectral_visc_prp)


      do iElemDeg=1,nElems*mpd1
        dof = mod(iElemDeg-1,mpd1) + 1
        iElem = (iElemDeg-1)/mpd1 + 1
        dofAbs = sqrt( ((real(dof,rk)-1.0_rk)**2) / (cut**2) )
        damping = exp( (-filter%alpha * (dofAbs**filter%order) ) )
        state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
      end do


    case default
      write(logUnit(1),*) 'ERROR in atl_spectral_visc_1d: Unknown filter function, stopping '
      call tem_abort()
    end select

  end subroutine atl_spectral_visc_1d

  !> Damp )the modal coefficients of the state vector by a given spectral
  !! viscosity method.
  subroutine atl_spectral_visc_2d( state, mesh, filter, maxPolyDeg, orders )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_spectral_visc_type), intent(in) :: filter
    integer, intent(in) :: maxPolyDeg
    real(kind=rk), intent(in), optional :: orders(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDegX, iDegY, dof
    real(kind=rk) :: dofAbs, damping, cut, dofAbsX, dofAbsY
    integer :: iEXY, mpd1, mpd1_square, nElems
    ! --------------------------------------------------------------------------

    ! Check for single degree of freedom to avoid divisions by zeros in the
    ! filter.
    if(maxPolyDeg.le.0) then
      return
    end if

    if(filter%cut_order .le. 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    mpd1 = maxPolyDeg+1
    mpd1_square = mpd1**2
    nElems = mesh%descriptor%elem%nElems(eT_fluid)

    select case(filter%kind)
    case(atl_exp_spectral_visc_prp)


      if(present(orders)) then

        do iEXY=1,nElems*mpd1_square
          iElem = (iEXY-1)/mpd1_square + 1
          dof = iEXY - (iElem-1)*mpd1_square
          iDegX = mod(dof-1, mpd1) + 1
          iDegY = (dof-1)/mpd1 + 1
          dofAbs = sqrt((((real(iDegX,rk) - 1.0_rk)**2) &
            & + ((real(iDegY,rk)-1.0_rk)**2)) / (cut**2) )
          damping = exp( (-filter%alpha * (dofAbs**orders(iElem)) ) )
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      else

        do iEXY=1,nElems*mpd1_square
          iElem = (iEXY-1)/mpd1_square + 1
          dof = iEXY - (iElem-1)*mpd1_square
          iDegX = mod(dof-1, mpd1) + 1
          iDegY = (dof-1)/mpd1 + 1
          dofAbs = sqrt((((real(iDegX,rk) - 1.0_rk)**2) &
            & + ((real(iDegY,rk)-1.0_rk)**2)) / (cut**2) )
          damping = exp( (-filter%alpha * (dofAbs**filter%order) ) )
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      end if

    case(atl_poly_spectral_visc_prp)

      if(present(orders)) then

        do iEXY=1,nElems*mpd1_square
          iElem = (iEXY-1)/mpd1_square + 1
          dof = iEXY - (iElem-1)*mpd1_square
          iDegX = mod(dof-1, mpd1) + 1
          iDegY = (dof-1)/mpd1 + 1
          dofAbsX = (real(iDegX,rk) - 1.0_rk) / cut
          dofAbsY = (real(iDegY,rk) - 1.0_rk) / cut
          damping = (1.0_rk - (dofAbsX**orders(iElem)))*(1.0_rk-(dofAbsY**orders(iElem)))
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      else

        do iEXY=1,nElems*mpd1_square
          iElem = (iEXY-1)/mpd1_square + 1
          dof = iEXY - (iElem-1)*mpd1_square
          iDegX = mod(dof-1, mpd1) + 1
          iDegY = (dof-1)/mpd1 + 1
          dofAbsX = (real(iDegX,rk) - 1.0_rk) / cut
          dofAbsY = (real(iDegY,rk) - 1.0_rk) / cut
          damping = (1.0_rk - (dofAbsX**filter%order))*(1.0_rk-(dofAbsY**filter%order))
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      end if

    case default
      write(logUnit(1),*) 'ERROR in atl_spectral_visc_2d: Unknown filter function, stopping '
      call tem_abort()
    end select

  end subroutine atl_spectral_visc_2d

  !> Damp )the modal coefficients of the state vector by a given spectral
  !! viscosity method.
  subroutine atl_spectral_visc_3d( state, mesh, filter, maxPolyDeg, orders )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_spectral_visc_type), intent(in) :: filter
    integer, intent(in) :: maxPolyDeg
    real(kind=rk), intent(in), optional :: orders(:)
    ! --------------------------------------------------------------------------
    integer :: iElem, iDegX, iDegY, iDegZ, dof
    real(kind=rk) :: dofAbs, damping, cut, dofAbsX, dofAbsY, dofAbsZ
    integer :: nElems, mpd1, mpd1_square, mpd1_cube
    integer :: iEXYZ
    ! --------------------------------------------------------------------------

    ! Check for single degree of freedom to avoid divisions by zeros in the
    ! filter.
    if(maxPolyDeg.le.0) then
      return
    end if

    if (filter%cut_order .le. 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    nElems = mesh%descriptor%elem%nElems(eT_fluid)
    mpd1 = maxPolyDeg+1
    mpd1_square = mpd1**2
    mpd1_cube = mpd1_square*mpd1

    select case(filter%kind)
    case(atl_exp_spectral_visc_prp)

      if(present(orders)) then

        do iEXYZ=1,nElems*mpd1_cube
          iElem = (iEXYZ-1)/mpd1_cube + 1
          dof = iEXYZ - (iElem-1)*mpd1_cube
          iDegX = mod(dof-1, mpd1) + 1
          iDegZ = (dof-1)/mpd1_square + 1
          iDegY = (dof - (iDegZ-1)*mpd1_square - 1)/mpd1 + 1
          dofAbs = sqrt( &
                         & ( &
                         &      ((real(iDegX,rk)-1.0_rk)**2) &
                         &    + ((real(iDegY,rk)-1.0_rk)**2) &
                         &    + ((real(iDegZ,rk)-1.0_rk)**2) &
                         & ) &
                         &  / (cut**2) &
                         &  )
          damping = exp( -filter%alpha * (real(dofAbs,rk)**orders(iElem)) )
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      else

        do iEXYZ=1,nElems*mpd1_cube
          iElem = (iEXYZ-1)/mpd1_cube + 1
          dof = iEXYZ - (iElem-1)*mpd1_cube
          iDegX = mod(dof-1, mpd1) + 1
          iDegZ = (dof-1)/mpd1_square + 1
          iDegY = (dof - (iDegZ-1)*mpd1_square - 1)/mpd1 + 1
          dofAbs = sqrt( &
                         & ( &
                         &      ((real(iDegX,rk)-1.0_rk)**2) &
                         &    + ((real(iDegY,rk)-1.0_rk)**2) &
                         &    + ((real(iDegZ,rk)-1.0_rk)**2) &
                         & ) &
                         &  / (cut**2) &
                         &  )
          damping = exp( -filter%alpha * (real(dofAbs,rk)**filter%order) )
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      end if

    case(atl_poly_spectral_visc_prp)

      if(present(orders)) then

        do iEXYZ=1,nElems*mpd1_cube
          iElem = (iEXYZ-1)/mpd1_cube + 1
          dof = iEXYZ - (iElem-1)*mpd1_cube
          iDegX = mod(dof-1, mpd1) + 1
          iDegZ = (dof-1)/mpd1_square + 1
          iDegY = (dof - (iDegZ-1)*mpd1_square - 1)/mpd1 + 1
          dofAbsX = (real(iDegX,rk) - 1.0_rk) / cut
          dofAbsY = (real(iDegY,rk) - 1.0_rk) / cut
          dofAbsZ = (real(iDegZ,rk) - 1.0_rk) / cut
          damping = (1.0_rk-(dofAbsX**orders(iElem))) &
                & * (1.0_rk-(dofAbsY**orders(iElem))) &
                & * (1.0_rk-(dofAbsZ**orders(iElem)))
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      else

        do iEXYZ=1,nElems*mpd1_cube
          iElem = (iEXYZ-1)/mpd1_cube + 1
          dof = iEXYZ - (iElem-1)*mpd1_cube
          iDegX = mod(dof-1, mpd1) + 1
          iDegZ = (dof-1)/mpd1_square + 1
          iDegY = (dof - (iDegZ-1)*mpd1_square - 1)/mpd1 + 1
          dofAbsX = (real(iDegX,rk) - 1.0_rk) / cut
          dofAbsY = (real(iDegY,rk) - 1.0_rk) / cut
          dofAbsZ = (real(iDegZ,rk) - 1.0_rk) / cut
          damping = (1.0_rk-(dofAbsX**filter%order)) &
                & * (1.0_rk-(dofAbsY**filter%order)) &
                & * (1.0_rk-(dofAbsZ**filter%order))
          state%state(iElem, dof, :) = damping*state%state(iElem, dof, :)
        end do

      end if

    case default
      write(logUnit(1),*) 'ERROR in atl_spectral_visc_3d: Unknown filter function, stopping '
      call tem_abort()
    end select

  end subroutine atl_spectral_visc_3d

  !> Damp the modal coefficients of the state vector by a given spectral
  !! viscosity method.
  subroutine atl_cheb_spectral_visc_1d( state, mesh, filter, poly_proj, &
    & maxPolyDeg )
    ! --------------------------------------------------------------------------
    type(atl_statedata_type), intent(inout) :: state
    type(atl_cube_elem_type), intent(in) :: mesh
    type(atl_spectral_visc_type), intent(in) :: filter
    integer, intent(in) :: maxPolyDeg
    type(ply_poly_project_type), intent(inout) :: poly_proj
    ! --------------------------------------------------------------------------
    integer :: iElem, iDof, iVar, nVars
    real(kind=rk) :: dofAbs, damping, cut
    real(kind=rk), allocatable :: legCoeffs(:), chebCoeffs(:)
    integer :: nElems, mpd1
    ! --------------------------------------------------------------------------

    ! Check for single degree of freedom to avoid divisions by zeros in the
    ! filter.
    if(maxPolyDeg.le.0) then
      return
    end if

    if(filter%cut_order .le. 0.0_rk) then
      cut = real(maxPolyDeg,rk)
    else
      cut = filter%cut_order
    end if

    if(.not.(poly_proj%kind.eq.'fpt')) then
      write(logUnit(1),*) 'ERROR in atl_cheb_spectral_visc_1d: &
        & Chebyshev-spectral viscosity requires initialized FPT, stopping ...'
      call tem_abort()
    end if


    nElems = mesh%descriptor%elem%nElems(eT_fluid)
    mpd1 = maxPolyDeg+1
    nVars = size(state%state,3)

    allocate(                                        &
            & chebCoeffs(mpd1), &
            & legCoeffs(mpd1)   &
            & )

    select case(filter%kind)
    case(atl_exp_spectral_visc_prp)


      do iElem=1,nElems
        do iVar = 1, nVars

          ! Copy the Legendre coefficients to dedicated array
          do iDof = 1, mpd1
            legCoeffs(iDof) = state%state(iElem, iDof, iVar)
          end do

          ! Transform Legendre expansion to Chebyshev expansion
          call ply_fpt_exec_striped( alph    = legCoeffs,                &
            &                        gam     = chebCoeffs,               &
            &                        nIndeps = 1,                        &
            &                        params  = poly_proj%body_1d         &
            &                                           %fpt             &
            &                                           %legToChebParams )

          do iDof = 1, mpd1
            dofAbs = sqrt( &
                           & ((real(iDof,rk)-1.0_rk)**2) / (cut**2) &
                           &  )
            damping = exp( (-filter%alpha * (dofAbs**filter%order) ) )
            chebCoeffs(iDof) = damping*chebCoeffs(iDof)
          end do

          ! Transform Chebyshev polynomials to Legendre polynomial
          call ply_fpt_exec_striped( alph    = chebCoeffs,               &
            &                        gam     = legCoeffs,                &
            &                        nIndeps = 1,                        &
            &                        params  = poly_proj%body_1d         &
            &                                           %fpt             &
            &                                           %chebToLegParams )

          ! copy back to the state array
          do iDof = 1, mpd1
            state%state(iElem, iDof,iVar) = legCoeffs(iDof)
          end do

        end do
      end do


    case default
      write(logUnit(1),*) 'ERROR in atl_cheb_spectral_visc_1d: Unknown filter function, stopping '
      call tem_abort()
    end select

  end subroutine atl_cheb_spectral_visc_1d

end module atl_stabilize_module
