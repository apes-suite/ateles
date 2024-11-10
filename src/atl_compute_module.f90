! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011, 2014 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2011-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Daniel Harlacher <daniel.harlacher@uni-siegen.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Michael Gaida  <michael.gaida@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2017-2018 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
! Copyright (c) 2018 Robin Weihe <robin.weihe@student.uni-siegen.de>
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

!> summary: module including the computation loops.
!! author: Jens Zudrop
!!
!! All routines, which are usually used inside the time step
!! iteration loop, are collected inside this module.
!! The computation of the right hand side is decomposed
!! into 3 main steps: pre-process, compute_rhs and post-process.
module atl_compute_module
  use env_module,                 only: rk

  ! Treelm modules
  use tem_aux_module,             only: tem_abort
  use tem_general_module,         only: tem_general_type
  use tem_element_module,         only: eT_fluid
  use tem_timer_module,           only: tem_startTimer, tem_stopTimer

  use atl_timer_module,           only: atl_timerHandles, atl_elemTimers
  use atl_eqn_linearEuler_module, only: atl_eqn_update_background
  use atl_materialPrp_module,     only: atl_ConstMatIdx, atl_VarMatIdx
  use atl_source_module,          only: atl_update_sourcedata


  !HK: WORKAROUND for Intel compiler, moved most used modules into the subroutines
  !    themselves, to avoid internal compiler error.

  implicit none

  private

  public :: atl_preprocess_rhs, atl_postprocess_rhs, atl_compute_rhs

  ! Intel 15 workaround:
  ! This exposure is needed to have the subroutine accessible for unit tests.
  public :: atl_initialize_state_der

  interface atl_preprocess_rhs
    module procedure preprocess_rhs_cubes
  end interface

  interface atl_postprocess_rhs
    module procedure postprocess_rhs_cubes
  end interface


  interface atl_compute_rhs
    module procedure compute_rhs_cubes
  end interface


contains


  ! ************************************************************************ !
  !> summary: compute the right hand side of your discrete equation.
  subroutine preprocess_rhs_cubes( minLevel, maxLevel, currentLevel,    &
    & mesh_list, tree, statedata_list, facedata_list, boundary_list,    &
    & bc, scheme_list, poly_proj_list, equation, material_list, general )

    use atl_cube_elem_module,           only: atl_cube_elem_type
    use atl_kerneldata_module,          only: atl_kerneldata_type, &
      &                                       atl_statedata_type
    use atl_boundary_module,            only: atl_level_boundary_type
    use atl_source_types_module,        only: atl_source_type
    use atl_scheme_module,              only: atl_scheme_type,        &
      &                                       atl_modg_scheme_prp,    &
      &                                       atl_modg_2D_scheme_prp, &
      &                                       atl_modg_1D_scheme_prp
    use atl_equation_module,            only: atl_equations_type
    use atl_bc_header_module,           only: atl_boundary_type
    use atl_modg_kernel_module,         only: atl_preprocess_modg_kernel, &
      &                                       atl_modg_modalVolToModalFace, &
      &                                       atl_modg_ensure_pos_facemean
    use atl_modg_multilevel_module,     only: atl_modg_coarseToFineFace
    use atl_modg_bnd_module,            only: atl_modg_set_bnd
    use atl_voltoface_module,           only: atl_modg_1d_modalVolToModalFace
    use atl_modg_1d_multilevel_module,  only: atl_modg_1d_coarseToFineFace
    use atl_modg_1d_kernel_module,      only: atl_modg_1d_ensure_pos_face, &
      &                                       atl_preprocess_modg_1d_kernel
    use atl_modg_1d_bnd_module,         only: atl_modg_1d_set_bnd
    use atl_modg_2d_kernel_module,      only: atl_preprocess_modg_2d_kernel, &
      &                                       atl_modg_2d_modalVolToModalFace, &
      &                                       atl_modg_2d_ensure_pos_facemean
    use atl_modg_2d_multilevel_module,  only: atl_modg_2d_coarseToFineFace
    use atl_modg_2d_bnd_module,         only: atl_modg_2d_set_bnd
    use atl_elemental_time_integration_module,                    &
      &                                 only: atl_timestep_type
    use atl_facedata_module,            only: atl_facedata_type
    use atl_materialPrp_module,         only: atl_material_type
    use ply_poly_project_module,        only: ply_poly_project_type, &
      &                                       assignment(=)
    use ply_dynArray_project_module,    only: dyn_projectionArray_type
    use treelmesh_module,               only: treelmesh_type
    use tem_comm_env_module,            only: tem_comm_env_type

    ! -------------------------------------------------------------------- !
    !> The minimum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> the level to compute on
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type),  intent(inout) :: mesh_list(minLevel:maxLevel)
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list(minLevel:maxLevel)
    !> List of face states you want to calc the rhs for. For each level we have
    !! one.
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of boundaries, for each level.
    type(atl_level_boundary_type), intent(inout) :: &
      & boundary_list(minLevel:maxLevel)
    !> The global boundary description.
    type(atl_boundary_type), intent(in) :: bc(:)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> The equation you are operating with.
    type(atl_equations_type),intent(inout) :: equation
    !> Information about the material parameters of the equation.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !> General treelm settings.
    type(tem_general_type), intent(inout) :: general
    ! -------------------------------------------------------------------- !
    integer :: iDir, iFace
    ! -------------------------------------------------------------------- !

    select case(scheme_list(currentLevel)%scheme)
    case(atl_modg_scheme_prp) ! MODG kernel
      select case(trim(equation%eq_kind))
      case('maxwell', 'maxwelldivcorrection', 'euler', 'navier_stokes', &
        &  'filtered_navier_stokes', 'acoustic', 'heat', 'lineareuler', &
        &  'loclineuler')

        call tem_startTimer( timerHandle = atl_timerHandles%preprocessKernel )

        ! For linear euler we need to update the background in each timestep and
        ! intermediate timestep
        select case(trim(equation%eq_kind))
        case('lineareuler')
          call tem_startTimer( timerHandle = atl_timerHandles%updateBackground)
          call atl_eqn_update_background(                            &
            & me          = equation%linearEuler,                    &
            & time        = statedata_list(currentLevel)%local_time, &
            & nDimensions = 3                                        )
          call tem_stopTimer( timerHandle = atl_timerHandles%updateBackground)
        end select

        ! Update the source terms.
        call atl_preprocess_modg_kernel(                          &
          &    equation           = equation,                     &
          &    statedata          = statedata_list(currentLevel), &
          &    mesh               = mesh_list(currentLevel),      &
          &    boundary           = boundary_list(currentLevel),  &
          &    scheme             = scheme_list(currentLevel),    &
          &    material           = material_list(currentLevel),  &
          &    poly_proj_material = poly_proj_list(               &
          &                           material_list(currentLevel) &
          &                             %poly_proj_pos),          &
          &    commPattern        = general%commPattern,          &
          &    proc               = general%proc                  )

        call tem_stopTimer( timerHandle = atl_timerHandles%preprocessKernel )

        call tem_startTimer( timerHandle = atl_timerHandles%projectToFace )

        ! Project modal representation to a modal representation on the faces.
        call atl_modg_modalVolToModalFace(                           &
          & nElems_fluid = mesh_list(currentLevel)%descriptor        &
          &                                       %elem              &
          &                                       %nElems(eT_fluid), &
          & length       = mesh_list(currentLevel)%length,           &
          & volState     = statedata_list(currentLevel)%state,       &
          & faceRep      = facedata_list(currentLevel)%faceRep,      &
          & nScalars     = equation%varSys%nScalars,                 &
          & nDerivatives = equation%nDerivatives,                    &
          & modg         = scheme_list(currentLevel)%modg            )

        call tem_stopTimer( timerHandle = atl_timerHandles%projectToFace )

        call tem_startTimer( timerHandle = atl_timerHandles%setBnd )

        ! Now, we can define the values for the boundaries. We do so, by diretly
        ! imposing face values on the boundaries in a weak sense (i.e. set modal
        ! face representation for the virtual boundary elements).
        select case(trim(equation%eq_kind))
        case('maxwell','maxwelldivcorrection')
          ! Set the boundaries for all faces with purely constant
          ! material parameters (done in non-nodal way)
          call atl_modg_set_bnd(                                               &
            & bc           = bc,                                               &
            & boundary     = material_list(currentLevel)%material_desc         &
            &                                           %bnd_faces(1)          &
            &                                           %boundary,             &
            & facedata     = facedata_list(currentLevel),                      &
            & poly_proj    = poly_proj_list(                                   &
            &                  boundary_list(currentLevel)%poly_proj_pos),     &
            & equation     = equation,                                         &
            & time         = statedata_list(currentLevel)%local_time,          &
            & mesh         = mesh_list(currentLevel),                          &
            & nodalBnd     = .false.,                                          &
            & material     = material_list(currentLevel)%material_dat          &
            &                                           %faceMaterialData(:,1),&
            & statedata    = statedata_list(currentLevel),                     &
            & currentLevel = currentLevel                                      )
          ! Set the boundaries for all faces with variable
          ! material parameters (done in non-nodal way)
          call atl_modg_set_bnd(                                               &
            & bc           = bc,                                               &
            & boundary     = material_list(currentLevel)%material_desc         &
            &                                           %bnd_faces(2)          &
            &                                           %boundary,             &
            & facedata     = facedata_list(currentLevel),                      &
            & poly_proj    = poly_proj_list(                                   &
            &                  boundary_list(currentLevel)%poly_proj_pos),     &
            & equation     = equation,                                         &
            & time         = statedata_list(currentLevel)%local_time,          &
            & mesh         = mesh_list(currentLevel),                          &
            & nodalBnd     = .true.,                                           &
            & material     = material_list(currentLevel)%material_dat          &
            &                                           %faceMaterialData(:,2),&
            & statedata    = statedata_list(currentLevel),                     &
            & currentLevel = currentLevel                                      )

        case('euler', 'navier_stokes', &
          & 'filtered_navier_stokes',  &
          & 'loclineuler'              )
          if (equation%euler%ensure_positivity) then
            call atl_modg_ensure_pos_facemean(                           &
              & nelems_fluid      = mesh_list(currentLevel)              &
              &                     %descriptor                          &
              &                     %elem                                &
              &                     %nElems(eT_fluid),                   &
              & volState          = statedata_list(currentLevel)%state,  &
              & faceRep           = facedata_list(currentLevel)%faceRep, &
              & nScalars          = equation%varSys%nScalars,            &
              & ensure_positivity = [.true.,                             &
              &                      .false., .false., .false.,          &
              &                      .true. ]                            )
          end if
          call atl_modg_set_bnd(                                           &
            & bc           = bc,                                           &
            & boundary     = boundary_list(currentLevel),                  &
            & facedata     = facedata_list(currentLevel),                  &
            & poly_proj    = poly_proj_list(                               &
            &                  boundary_list(currentLevel)%poly_proj_pos), &
            & equation     = equation,                                     &
            & time         = statedata_list(currentLevel)%local_time,      &
            & mesh         = mesh_list(currentLevel),                      &
            & nodalBnd     = equation%isNonlinear,                         &
            & statedata    = statedata_list(currentLevel),                 &
            & currentLevel = currentLevel                                  )

        case('heat', 'lineareuler', 'acoustic')
          call atl_modg_set_bnd(                                           &
            & bc           = bc,                                           &
            & boundary     = boundary_list(currentLevel),                  &
            & facedata     = facedata_list(currentLevel),                  &
            & poly_proj    = poly_proj_list(                               &
            &                  boundary_list(currentLevel)%poly_proj_pos), &
            & equation     = equation,                                     &
            & time         = statedata_list(currentLevel)%local_time,      &
            & mesh         = mesh_list(currentLevel),                      &
            & nodalBnd     = equation%isNonlinear,                         &
            & statedata    = statedata_list(currentLevel),                 &
            & currentLevel = currentLevel                                  )

        case default
          call tem_abort( 'ERROR in preprocess_rhs_cubes: unknown equation' &
            & // ' for MODG boundary handling, stopping ...'                )
        end select
        call tem_stopTimer( timerHandle = atl_timerHandles%setBnd )

        ! Now, we exchange all the face representations (for all 3 directions left and
        ! right face representations).
        call tem_startTimer( timerHandle = atl_timerHandles%commState )
        do iDir = 1,3
          do iFace = 1,2
            call general%commPattern%exchange_real(                           &
              & send         = mesh_list(currentLevel)                        &
              &                  %faces                                       &
              &                  %faces(iDir)                                 &
              &                  %sendBuffer_state(iFace),                    &
              & recv         = mesh_list(currentLevel)                        &
              &                  %faces                                       &
              &                  %faces(iDir)                                 &
              &                  %recvBuffer_state(iFace),                    &
              & state        = facedata_list(currentLevel)%faceRep(iDir)%dat, &
              & message_flag = currentLevel,                                  &
              & comm         = general%proc%comm                              )
          end do
        end do
        call tem_stopTimer( timerHandle = atl_timerHandles%commState )

        ! Interpolate face representation to next finer level
        if(currentLevel.lt.maxLevel) then
          call atl_modg_coarseToFineFace(                          &
            & minLevel     = minLevel,                             &
            & maxLevel     = maxLevel,                             &
            & currentLevel = currentLevel,                         &
            & mesh         = mesh_list,                            &
            & facedata     = facedata_list,                        &
            & scheme       = scheme_list,                          &
            & nScalars     = equation%varSys%nScalars              &
            &                  * ( 3 * equation%nDerivatives + 1 ) )
        end if

      case default
        call tem_abort( 'ERROR in preprocess_rhs_cubes: unknown equation for' &
          & // ' MODG, stopping...'                                           )
      end select

    case(atl_modg_2d_scheme_prp) ! 2D MODG kernel
      select case(trim(equation%eq_kind))
      case('maxwell_2d','pec_maxwell_2d', 'euler_2d', 'navier_stokes_2d', &
        & 'acoustic_2d', 'filtered_navier_stokes_2d', 'heat_2d',          &
        & 'lineareuler_2d'                                                )

        call tem_startTimer( timerHandle = atl_timerHandles%preprocessKernel )

        ! For linear euler we need to update the temporal background in
        ! every timestep
        select case(trim(equation%eq_kind))
        case('lineareuler_2d')
          call tem_startTimer( timerHandle = atl_timerHandles%updateBackground)
          call atl_eqn_update_background(                            &
            & me          = equation%linearEuler,                    &
            & time        = statedata_list(currentLevel)%local_time, &
            & nDimensions = 2                                        )
        call tem_stopTimer( timerHandle = atl_timerHandles%updateBackground)
        end select

        call atl_preprocess_modg_2d_kernel(                        &
          &    equation           = equation,                      &
          &    statedata          = statedata_list(currentLevel),  &
          &    mesh_list          = mesh_list,                     &
          &    boundary_list      = boundary_list,                 &
          &    scheme_list        = scheme_list ,                  &
          &    material_list      = material_list,                 &
          &    currentLevel       = currentLevel,                  &
          &    minLevel           = minLevel,                      &
          &    maxLevel           = maxLevel,                      &
          &    poly_proj_material = poly_proj_list(                &
          &                           material_list(currentLevel)  &
          &                             %poly_proj_pos),           &
          &    commPattern        = general%commPattern,           &
          &    Proc               = general%Proc                   )

        call tem_stopTimer( timerHandle = atl_timerHandles%preprocessKernel )

        call tem_startTimer( timerHandle = atl_timerHandles%projectToFace )

        ! Project modal representation to a modal representation on the faces.
        call atl_modg_2d_modalVolToModalFace(             &
          & mesh      = mesh_list(currentLevel),          &
          & statedata = statedata_list(currentLevel),     &
          & facedata  = facedata_list(currentLevel),      &
          & equation  = equation,                         &
          & modg_2d   = scheme_list(currentLevel)%modg_2d )

        call tem_stopTimer( timerHandle = atl_timerHandles%projectToFace )

        call tem_startTimer( timerHandle = atl_timerHandles%setBnd )

        select case(trim(equation%eq_kind))
        case('maxwell_2d')
          ! Set the boundaries for all faces with purely constant
          ! material parameters (done in non-nodal way)
          call atl_modg_2d_set_bnd(                                            &
            & bc           = bc,                                               &
            & boundary     = material_list(currentLevel)%material_desc         &
            &                                           %bnd_faces(1)          &
            &                                           %boundary,             &
            & facedata     = facedata_list(currentLevel),                      &
            & statedata    = statedata_list(currentLevel),                     &
            & modg         = scheme_list(currentLevel)%modg_2d,                &
            & poly_proj    = poly_proj_list(                                   &
            &                  boundary_list(currentLevel)%poly_proj_pos),     &
            & equation     = equation,                                         &
            & time         = statedata_list(currentLevel)%local_time,          &
            & mesh         = mesh_list(currentLevel),                          &
            & nodalBnd     = .false.,                                          &
            & material     = material_list(currentLevel)%material_dat          &
            &                                           %faceMaterialData(:,1),&
            & currentLevel = currentLevel                                      )
          ! Set the boundaries for all faces with variable
          ! material parameters (done in nodal way)
          call atl_modg_2d_set_bnd(                                            &
            & bc           = bc,                                               &
            & boundary     = material_list(currentLevel)%material_desc         &
            &                                           %bnd_faces(2)          &
            &                                           %boundary,             &
            & facedata     = facedata_list(currentLevel),                      &
            & statedata    = statedata_list(currentLevel),                     &
            & modg         = scheme_list(currentLevel)%modg_2d,                &
            & poly_proj    = poly_proj_list(                                   &
            &                  boundary_list(currentLevel)%poly_proj_pos),     &
            & equation     = equation,                                         &
            & time         = statedata_list(currentLevel)%local_time,          &
            & mesh         = mesh_list(currentLevel),                          &
            & nodalBnd     = .true.,                                           &
            & material     = material_list(currentLevel)%material_dat          &
            &                                           %faceMaterialData(:,2),&
            & currentLevel = currentLevel                                      )
        case('euler_2d', 'navier_stokes_2d', &
          &  'filtered_navier_stokes_2d'     )
          ! Now, we can define the values for the boundaries. We do so, by
          ! diretly imposing face values on the boundaries in a weak sense
          ! (i.e. set modal face representation for the virtual boundary
          ! elements).
          if (equation%euler%ensure_positivity) then
            call atl_modg_2d_ensure_pos_facemean(                        &
              & nelems_fluid      = mesh_list(currentLevel)              &
              &                     %descriptor                          &
              &                     %elem                                &
              &                     %nElems(eT_fluid),                   &
              & volState          = statedata_list(currentLevel)%state,  &
              & faceRep           = facedata_list(currentLevel)%faceRep, &
              & nScalars          = equation%varSys%nScalars,            &
              & ensure_positivity = [.true.,                             &
              &                      .false., .false.,                   &
              &                      .true. ]                            )
          end if
          call atl_modg_2d_set_bnd(                                        &
            & bc           = bc,                                           &
            & boundary     = boundary_list(currentLevel),                  &
            & facedata     = facedata_list(currentLevel),                  &
            & statedata    = statedata_list(currentLevel),                 &
            & modg         = scheme_list(currentLevel)%modg_2d,            &
            & poly_proj    = poly_proj_list(                               &
            &                  boundary_list(currentLevel)%poly_proj_pos), &
            & equation     = equation,                                     &
            & time         = statedata_list(currentLevel)%local_time,      &
            & mesh         = mesh_list(currentLevel),                      &
            & nodalBnd     = equation%isNonlinear,                         &
            & currentLevel = currentLevel                                  )
        case('acoustic_2d','heat_2d', &
          &  'lineareuler_2d'         )
          ! Now, we can define the values for the boundaries. We do so, by
          ! diretly imposing face values on the boundaries in a weak sense
          ! (i.e. set modal face representation for the virtual boundary
          ! elements).
          call atl_modg_2d_set_bnd(                                        &
            & bc           = bc,                                           &
            & boundary     = boundary_list(currentLevel),                  &
            & facedata     = facedata_list(currentLevel),                  &
            & statedata    = statedata_list(currentLevel),                 &
            & modg         = scheme_list(currentLevel)%modg_2d,            &
            & poly_proj    = poly_proj_list(                               &
            &                  boundary_list(currentLevel)%poly_proj_pos), &
            & equation     = equation,                                     &
            & time         = statedata_list(currentLevel)%local_time,      &
            & mesh         = mesh_list(currentLevel),                      &
            & nodalBnd     = equation%isNonlinear,                         &
            & currentLevel = currentLevel                                  )
        case default
          call tem_abort( 'ERROR in preprocess_rhs_cubes: unknown equation' &
            & // ' for 2D MODG boundary handling, stopping ...'             )
        end select
        call tem_stopTimer( timerHandle = atl_timerHandles%setBnd )

        ! Now, we exchange all the face representations (only for 2 directions,
        ! i.e. x and y)
        call tem_startTimer( timerHandle = atl_timerHandles%commState )
        do iDir = 1,2
          do iFace = 1,2
            call general%commPattern%exchange_real(                           &
              & send         = mesh_list(currentLevel)                        &
              &                  %faces                                       &
              &                  %faces(iDir)                                 &
              &                  %sendBuffer_state(iFace),                    &
              & recv         = mesh_list(currentLevel)                        &
              &                  %faces                                       &
              &                  %faces(iDir)                                 &
              &                  %recvBuffer_state(iFace),                    &
              & state        = facedata_list(currentLevel)%faceRep(iDir)%dat, &
              & message_flag = currentLevel,                                  &
              & comm         = general%proc%comm                              )
          end do
        end do
        call tem_stopTimer( timerHandle = atl_timerHandles%commState )

        ! Interpolate face representation to next finer level
        if(currentLevel.lt.maxLevel) then
          call atl_modg_2d_coarseToFineFace(                       &
            & minLevel     = minLevel,                             &
            & maxLevel     = maxLevel,                             &
            & currentLevel = currentLevel,                         &
            & mesh         = mesh_list,                            &
            & facedata     = facedata_list,                        &
            & scheme       = scheme_list,                          &
            & nScalars     = equation%varSys%nScalars              &
            &                  * ( 2 * equation%nDerivatives + 1 ) )
        end if

      case default
        call tem_abort( 'ERROR in preprocess_rhs_cubes: unknown equation for' &
          & // ' 2d-MODG, stopping...'                                        )
      end select

    case(atl_modg_1d_scheme_prp) ! 1D MODG kernel
      select case(equation%eq_kind)
      case('euler_1d','advection_1d', 'heat_1d','loclineuler_1d')

        call tem_startTimer( timerHandle = atl_timerHandles%preprocessKernel )

        call tem_stopTimer( timerHandle = atl_timerHandles%preprocessKernel )

        call tem_startTimer( timerHandle = atl_timerHandles%projectToFace )

        call atl_preprocess_modg_1d_kernel(                       &
          &    equation           = equation,                     &
          &    statedata          = statedata_list(currentLevel), &
          &    mesh_list          = mesh_list,                    &
          &    boundary_list      = boundary_list,                &
          &    scheme_list        = scheme_list ,                 &
          &    material_list      = material_list,                &
          &    currentLevel       = currentLevel,                 &
          &    minLevel           = minLevel,                     &
          &    maxLevel           = maxLevel,                     &
          &    poly_proj_material = poly_proj_list(               &
          &                           material_list(currentLevel) &
          &                             %poly_proj_pos),          &
          &    commPattern        = general%commPattern,          &
          &    proc               = general%proc                  )
        ! Project modal representation to a modal representation on the faces.
        call atl_modg_1d_modalVolToModalFace(                                &
          & mesh          = mesh_list(currentLevel),                         &
          & statedata     = statedata_list(currentLevel),                    &
          & facedata      = facedata_list(currentLevel),                     &
          & nScalars      = equation%varSys%nScalars,                        &
          & maxPolyDegree = scheme_list(currentLevel)%modg_1d%maxPolyDegree, &
          & basisType     = scheme_list(currentLevel)%modg_1d%basisType,     &
          & equation      = equation                                         )

        if (equation%eq_kind == 'euler_1d') then
          if (equation%euler%ensure_positivity) then
            call atl_modg_1d_ensure_pos_face(                            &
              & nelems_fluid      = mesh_list(currentLevel)              &
              &                     %descriptor                          &
              &                     %elem                                &
              &                     %nElems(eT_fluid),                   &
              & volState          = statedata_list(currentLevel)%state,  &
              & faceRep           = facedata_list(currentLevel)%faceRep, &
              & nScalars          = equation%varSys%nScalars,            &
              & ensure_positivity = [ .true., .false., .true. ]          )
          end if
        end if

        call tem_stopTimer( timerHandle = atl_timerHandles%projectToFace )

        call tem_startTimer( timerHandle = atl_timerHandles%setBnd )

        ! Now, we can define the values for the boundaries. We do so, by diretly
        ! imposing face values on the boundaries in a weak sense (i.e. set modal
        ! face representation for the virtual boundary elements).
        call atl_modg_1d_set_bnd(                                     &
          & bc        = bc,                                           &
          & boundary  = boundary_list(currentLevel),                  &
          & facedata  = facedata_list(currentLevel),                  &
          & statedata = statedata_list(currentLevel),                 &
          & poly_proj = poly_proj_list(                               &
          &               boundary_list(currentLevel)%poly_proj_pos), &
          & equation  = equation,                                     &
          & tree      = tree,                                         &
          & time      = statedata_list(currentLevel)%local_time,      &
          & mesh      = mesh_list(currentLevel)                       )

        call tem_stopTimer( timerHandle = atl_timerHandles%setBnd )

        ! Now, we exchange all the face representations (only for 1 directions,
        ! i.e. x and y)
        call tem_startTimer( timerHandle = atl_timerHandles%commState )
        do iDir = 1,1
          do iFace = 1,2
            call general%commPattern%exchange_real(                           &
              & send         = mesh_list(currentLevel)                        &
              &                  %faces                                       &
              &                  %faces(iDir)                                 &
              &                  %sendBuffer_state(iFace),                    &
              & recv         = mesh_list(currentLevel)                        &
              &                  %faces                                       &
              &                  %faces(iDir)                                 &
              &                  %recvBuffer_state(iFace),                    &
              & state        = facedata_list(currentLevel)%faceRep(iDir)%dat, &
              & message_flag = currentLevel,                                  &
              & comm         = general%proc%comm                              )
          end do
        end do
        call tem_stopTimer( timerHandle = atl_timerHandles%commState )

        ! Interpolate face representation to next finer level
        if(currentLevel.lt.maxLevel) then
          call atl_modg_1d_coarseToFineFace(          &
            & minLevel     = minLevel,                &
            & maxLevel     = maxLevel,                &
            & currentLevel = currentLevel,            &
            & mesh         = mesh_list,               &
            & facedata     = facedata_list,           &
            & nScalars     = equation%varSys%nScalars )
        end if

      case default
        call tem_abort( 'ERROR in preprocess_rhs_cubes: unknown equation for' &
          & // ' 1d-MODG, stopping...'                                        )
      end select

    case default
      call tem_abort( 'ERROR in preprocess_rhs_cubes: unknown scheme,' &
        & // ' stopping...'                                            )
    end select

  end subroutine preprocess_rhs_cubes
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: Applies the postprocessing step of the compute step.
  subroutine postprocess_rhs_cubes( mesh, kerneldata, statedata, scheme, &
    &                               timestep, equation                   )
    use ply_poly_project_module,              only: ply_poly_project_type, &
      &                                             assignment(=)
    use atl_cube_elem_module,                 only: atl_cube_elem_type
    use atl_kerneldata_module,                only: atl_kerneldata_type, &
      &                                             atl_statedata_type
    use atl_scheme_module,                    only: atl_scheme_type,        &
      &                                             atl_modg_scheme_prp,    &
      &                                             atl_modg_2D_scheme_prp, &
      &                                             atl_modg_1d_scheme_prp
    use atl_equation_module,                  only: atl_equations_type
    use atl_modg_kernel_module,               only: atl_modg_invMassMatrix
    use atl_modg_2d_kernel_module,            only: atl_modg_2d_invMassMatrix
    use atl_modg_1d_kernel_module,            only: atl_modg_1d_invMassMatrix
    use atl_elemental_time_integration_module, only: atl_timestep_type
    ! -------------------------------------------------------------------- !
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme
    !> List of levelwise timestepping algorihtms
    type(atl_timestep_type), intent(inout) :: timestep
    !> The equation you are operating with.
    type(atl_equations_type),intent(in) :: equation
    ! -------------------------------------------------------------------- !

    call tem_startTimer( timerHandle = atl_timerHandles%invMassMatrix )

    select case(scheme%scheme)
    case(atl_modg_scheme_prp) ! MODG kernel
      call atl_modg_invMassMatrix( mesh              = mesh,                  &
        &                          kerneldata        = kerneldata,            &
        &                          statedata         = statedata,             &
        &                          elementalTimestep = timestep%elemStep_vec, &
        &                          timestep          = timestep,              &
        &                          scheme            = scheme%modg            )

    case(atl_modg_2d_scheme_prp) ! 2D MODG kernel
      call atl_modg_2d_invMassMatrix(                &
        & mesh              = mesh,                  &
        & equation          = equation,              &
        & kerneldata        = kerneldata,            &
        & statedata         = statedata,             &
        & elementalTimestep = timestep%elemStep_vec, &
        & timestep          = timestep,              &
        & scheme            = scheme%modg_2d         )

    case(atl_modg_1d_scheme_prp) ! 1D MODG kernel
      call atl_modg_1d_invMassMatrix(                &
        & mesh              = mesh,                  &
        & kerneldata        = kerneldata,            &
        & statedata         = statedata,             &
        & elementalTimestep = timestep%elemStep_vec, &
        & timestep          = timestep,              &
        & scheme            = scheme%modg_1d         )

    case default
      call tem_abort( 'ERROR in postprocess_rhs_cubes: postrocessing for this' &
        & // ' scheme is not supported yet!'                                   )
    end select

    call tem_stopTimer( timerHandle = atl_timerHandles%invMassMatrix )

  end subroutine postprocess_rhs_cubes
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> summary: compute the right hand side of your discrete equation.
  subroutine compute_rhs_cubes( minLevel, maxLevel, currentLevel, mesh_list, &
    & tree, kerneldata_list, statedata_list, facedata_list, source,          &
    & penalizationdata_list, scheme_list, poly_proj_pos, poly_proj_list,     &
    & equation, material_list, general, computePenalization                 )
    use atl_cube_elem_module,                   only: atl_cube_elem_type
    use atl_source_types_module,                only: atl_source_type
    use atl_source_module,                      only: atl_deallocate_sourceData
    use atl_scheme_module,                      only: atl_scheme_type,        &
      &                                               atl_modg_scheme_prp,    &
      &                                               atl_modg_2D_scheme_prp, &
      &                                               atl_modg_1d_scheme_prp
    use atl_equation_module,                    only: atl_equations_type
    use atl_elemental_time_integration_module,  only: atl_timestep_type
    use atl_facedata_module,                    only: atl_facedata_type
    use atl_kerneldata_module,                  only: atl_kerneldata_type, &
      &                                               atl_statedata_type
    use atl_materialPrp_module,                 only: atl_material_type
    use ply_poly_project_module,                only: ply_poly_project_type, &
      &                                               assignment(=)
    use ply_dynArray_project_module,            only: dyn_projectionArray_type
    use atl_penalization_module,                only: atl_penalizationData_type
    use treelmesh_module,                       only: treelmesh_type
    ! -------------------------------------------------------------------- !
    !> The minimum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The level of the mesh you are computing the rhs for.
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh_list(minLevel:maxLevel)
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: kerneldata_list(minLevel:maxLevel)
    !> List of facedatas. For each level we have one
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list(minLevel:maxLevel)
    !> Levelwise list of sources
    type(atl_source_type), intent(inout) :: source
    !> Levelwise list of penalization data
    type(atl_penalizationData_type), intent(inout) :: &
      & penalizationdata_list(minLevel:maxLevel)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> List of levelwise position of  projection method in unique
    !! projection_list
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> The equation you are operating with.
    type(atl_equations_type),intent(in) :: equation
    !> Material parameter description.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> Flag to indicate whether penalization terms should be computed or not.
    !!
    !! This is used to switch off the application of penalizing terms from
    !! the convective computation and instead compute it in an implicit
    !! update, see the atl_imexrk_module.
    !! Default is .true.!
    logical, intent(in), optional :: computePenalization
    ! -------------------------------------------------------------------- !

    select case(scheme_list(currentLevel)%scheme)
    case(atl_modg_scheme_prp, atl_modg_2d_scheme_prp, atl_modg_1d_scheme_prp)
      ! modg scheme
      ! No reconstruction necessary, so we simply skip the reconstruction.
      ! Furhtermore, for the modal Discontinuous Galerkin scheme we need the
      ! full modal information for each cell, so we do not communicate any cell
      ! values. We just use the cell values we got from the preprocess step.
    case default
      call tem_abort( 'ERROR in compute_rhs_cubes: this scheme is not' &
        & // ' supported yet (reconstruction)!'                        )
    end select

    select case(scheme_list(currentLevel)%scheme)
    case(atl_modg_scheme_prp) ! MODG scheme
      call compute_rhs_cubes_modg( minLevel, maxLevel, currentLevel,         &
        & mesh_list, kerneldata_list, statedata_list, facedata_list, source, &
        & penalizationdata_list, scheme_list, poly_proj_pos, poly_proj_list, &
        & equation, material_list, general, computePenalization              )

    case(atl_modg_2d_scheme_prp) ! 2D MODG scheme
      call compute_rhs_cubes_modg_2d( minLevel, maxLevel, currentLevel, &
        & mesh_list, kerneldata_list, statedata_list, facedata_list,    &
        & source, penalizationdata_list, scheme_list, poly_proj_pos,    &
        & poly_proj_list, equation, material_list, general,             &
        & computePenalization                                           )

    case(atl_modg_1d_scheme_prp) ! 1D MODG scheme
      call compute_rhs_cubes_modg_1d( minLevel, maxLevel, currentLevel,    &
        & mesh_list, tree, kerneldata_list, statedata_list, facedata_list, &
        & source, penalizationdata_list, scheme_list, poly_proj_pos,       &
        & poly_proj_list, equation, material_list, general,                &
        & computePenalization                                              )
    case default
      call tem_abort( 'ERROR in compute_rhs_cubes: this scheme is not' &
        & // ' supported yet (kernel)!'                                )
    end select

  end subroutine compute_rhs_cubes
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Computes the right hand side for cubical elements and MODG scheme.
  subroutine compute_rhs_cubes_modg( minLevel, maxLevel, currentLevel,   &
    & mesh_list, kerneldata_list, statedata_list, facedata_list, source, &
    & penalizationdata_list, scheme_list, poly_proj_pos,poly_proj_list,  &
    & equation, material_list, general, computePenalization              )
    use atl_cube_elem_module, only: &
      & atl_cube_elem_type
    use atl_source_types_module, only: &
      & atl_source_type
    use atl_scheme_module, only: &
      & atl_scheme_type
    use atl_equation_module, only: &
      & atl_equations_type
    use atl_project_physflux_module, only: &
      & atl_modg_project_PhysFlux_testFunc
    use atl_modg_kernel_module, only:     &
      & atl_modg_project_source,          &
      & atl_modg_project_NumFlux
    use atl_modg_multilevel_module, only: &
      & atl_modg_fineToCoarseFace
    use atl_modg_maxwell_kernel_module, only: &
      & atl_modg_maxwell_numflux,             &
      & atl_modg_maxwell_physflux_const,      &
      & atl_modg_maxwell_physflux_NonConst
    use atl_modg_maxwellDivCor_kernel_module, only: &
      & atl_modg_maxwellDivCor_numflux,             &
      & atl_modg_maxwellDivCor_physflux_const,      &
      & atl_modg_maxwellDivCor_physflux_NonConst
    use atl_modg_euler_kernel_module, only: &
      & atl_modg_euler_numflux,             &
      & atl_modg_euler_physFlux_const,      &
      & atl_modg_euler_physFlux_NonConst,   &
      & atl_modg_euler_penalization_const,  &
      & atl_modg_euler_penalization_Nonconst
    use atl_modg_navierstokes_kernel_module, only: &
      & atl_modg_navierstokes_numflux,             &
      & atl_modg_navierstokes_physFlux_const,      &
      & atl_modg_navierstokes_physFlux_Nonconst,   &
      & atl_modg_navierstokes_penalization_const,  &
      & atl_modg_navierstokes_penalization_Nonconst
    use atl_modg_filNvrStk_kernel_module, only: &
      & atl_modg_filNvrStk_numflux,             &
      & atl_modg_filNvrStk_physFlux_const,      &
      & atl_modg_filNvrStk_physFlux_Nonconst
    use atl_elemental_time_integration_module, only: &
      & atl_timestep_type
    use atl_facedata_module, only: &
      & atl_facedata_type
    use atl_kerneldata_module, only: &
      & atl_kerneldata_type,         &
      & atl_statedata_type
    use atl_materialPrp_module, only: &
      & atl_material_type
    use ply_poly_project_module, only: &
      & ply_poly_project_type,         &
      & assignment(=)
    use ply_dynArray_project_module, only: &
      & dyn_projectionArray_type
    use atl_modg_acoustic_kernel_module, only: &
      & atl_modg_acoustic_numflux,             &
      & atl_modg_acoustic_physFlux
    use atl_modg_linearEuler_kernel_module, only: &
      & atl_modg_linearEuler_numflux,             &
      & atl_modg_linearEuler_physFlux
    use atl_modg_LoclinEuler_kernel_module, only: &
      & atl_modg_LoclinEuler_physFlux
    use atl_modg_heat_kernel_module, only: &
      & atl_modg_heat_numflux,             &
      & atl_modg_heat_physFlux
    use atl_penalization_module, only: &
      & atl_penalizationData_type
    use atl_physFlux_module, only:  &
      & atl_physFlux_pointer_type,      &
      & atl_penalization_pointer_type
    use treelmesh_module, only: &
      & treelmesh_type
    ! -------------------------------------------------------------------- !
    !> The minimum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The level of the mesh you are computing the rhs for.
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type),  intent(inout) :: mesh_list(minLevel:maxLevel)
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: &
      & kerneldata_list(minLevel:maxLevel)
    !> List of facedatas. For each level we have one
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list(minLevel:maxLevel)
    !> Levelwise list of sources
    type(atl_source_type), intent(inout) :: source
    !> Levelwise list of penalization data
    type(atl_penalizationData_type), intent(inout) :: &
      & penalizationdata_list(minLevel:maxLevel)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> List of levelwise position of  projection method in unique list
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> The equation you are operating with.
    type(atl_equations_type),intent(in) :: equation
    !> Information about the material parameters of the equation.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> Flag to indicate whether penalization terms should be computed or not.
    !!
    !! This is used to switch off the application of penalizing terms from
    !! the convective computation and instead compute it in an implicit
    !! update, see the atl_imexrk_module.
    !! Default is .true.!
    logical, intent(in), optional :: computePenalization
    ! -------------------------------------------------------------------- !
    integer :: iDir, iFace
    integer :: dirVec(3,3)
    type(atl_physFlux_pointer_type) :: eval(2)
    type(atl_penalization_pointer_type) :: apply(2)
    logical :: usePenalization
    ! -------------------------------------------------------------------- !

    usePenalization = .true.
    if (present(computePenalization)) usePenalization = computePenalization

    ! calculate rhs of the PDE from source terms
    call atl_update_sourcedata(                                            &
      & equation     = equation,                                           &
      & time         = statedata_list(currentLevel)%local_time,            &
      & mesh         = mesh_list(currentLevel),                            &
      & poly_proj    = poly_proj_list(source%poly_proj_pos(currentLevel)), &
      & currentLevel = currentLevel,                                       &
      & state        = statedata_list(currentLevel)%state,                 &
      & material     = material_list(currentLevel),                        &
      & source       = source,                                             &
      & scheme       = scheme_list(currentLevel)                           )

    ! Interpolate fluxes from finer level to current level. Please keep in mind
    ! that the timestepping routine is called recursively for the levels.
    ! Therefore the compute rhs subroutine (i.e. the current routine) is
    ! executed on the finer levels before it is executed on the coarser level.
    ! So, the flux on the finer level is already available when a coarser level
    ! enters this routine.
    if( currentLevel .lt. maxLevel ) then
      call atl_modg_fineToCoarseFace(                                         &
        & minLevel     = minLevel,                                            &
        & maxLevel     = maxLevel,                                            &
        & currentLevel = currentLevel,                                        &
        & mesh         = mesh_list,                                           &
        & facedata     = facedata_list,                                       &
        & scheme       = scheme_list,                                         &
        & nScalars     = equation%varSys%nScalars * (equation%nDerivatives+1) )
    end if

    call tem_startTimer( timerHandle = atl_timerHandles%numFlux )
    ! Compute all the numerical fluxes, by the modal representations of the
    ! states on the faces. The modal representation on the faces is provided by
    ! the preprocess step.
    select case(trim(equation%eq_kind))
    case('maxwell')
      call atl_modg_maxwell_numFlux(                               &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg,              &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )
      !Set up the pointer for the physical fluxes
      eval(1)%physFlux  => atl_modg_maxwell_physflux_const
      eval(2)%physFlux  => atl_modg_maxwell_physflux_NonConst

    case('maxwelldivcorrection')
      call atl_modg_maxwellDivCor_numFlux(                         &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg,              &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )
      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_maxwellDivCor_physflux_const
      eval(2)%physFlux => atl_modg_maxwellDivCor_physflux_NonConst

    case('euler')
      call atl_modg_euler_numFlux(                                 &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux  => atl_modg_euler_physFlux_const
      eval(2)%physFlux  => atl_modg_euler_physFlux_NonConst

      !Set up the pointer for the penalization routines
      apply(1)%pen => atl_modg_euler_penalization_const
      apply(2)%pen => atl_modg_euler_penalization_Nonconst

    case('navier_stokes')
      call atl_modg_navierstokes_numFlux(                          &
        & mesh      = mesh_list(currentLevel),                     &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg,              &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_navierstokes_physflux_const
      eval(2)%physFlux => atl_modg_navierstokes_physflux_NonConst

      !Set up the pointer for the penalization routines
      apply(1)%pen => atl_modg_navierstokes_penalization_Const
      apply(2)%pen => atl_modg_navierstokes_penalization_NonConst

    case('filtered_navier_stokes')
      call atl_modg_filNvrStk_numFlux(                             &
        & mesh      = mesh_list(currentLevel),                     &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg,              &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_filNvrStk_physflux_const
      eval(2)%physFlux => atl_modg_filNvrStk_physflux_NonConst

    case('acoustic')
      call atl_modg_acoustic_numFlux(               &
        & mesh     = mesh_list(currentLevel),       &
        & equation = equation,                      &
        & facedata = facedata_list(currentLevel),   &
        & scheme   = scheme_list(currentLevel)%modg )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_acoustic_physFlux
      !eval(2)%physFlux

    case('lineareuler')
      call atl_modg_linearEuler_numFlux(            &
        & mesh     = mesh_list(currentLevel),       &
        & equation = equation,                      &
        & facedata = facedata_list(currentLevel),   &
        & scheme   = scheme_list(currentLevel)%modg )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_linearEuler_physFlux
      !eval(2)%physFlux

    case('loclineuler')
      call atl_modg_euler_numFlux(                                 &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_LoclinEuler_physFlux
      !eval(2)%physFlux

      !Set up the pointer for the penalization routines
      apply(1)%pen => atl_modg_euler_penalization_const
      apply(2)%pen => atl_modg_euler_penalization_Nonconst

    case('heat')
      call atl_modg_heat_numFlux(                                 &
        & mesh      = mesh_list(currentLevel),                    &
        & equation  = equation,                                   &
        & facedata  = facedata_list(currentLevel),                &
        & scheme    = scheme_list(currentLevel)%modg,             &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)) )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux  => atl_modg_heat_physFlux
      !eval(2)%physFlux

    case default
      call tem_abort( 'ERROR in compute_rhs_cubes: '                   &
        & // 'modg is not supporting this PDE (num flux calculation)!' )
    end select
    call tem_stopTimer( timerHandle = atl_timerHandles%numFlux )

    ! Exchange all the fluxes on the faces (for all 3 directions).
    ! We send data about the fluxes at those faces where we received
    ! information about the states (in preprocess step). Therefore,
    ! we use the send buffer for receiving and the receive buffer for
    ! sending.
    call tem_startTimer( timerHandle = atl_timerHandles%commState )
    do iDir = 1,3
      do iFace = 1,2
        call general%commPattern%exchange_real(                            &
          & send         = mesh_list(currentLevel)                         &
          &                  %faces                                        &
          &                  %faces(iDir)                                  &
          &                  %recvBuffer_flux(iFace),                      &
          & recv         = mesh_list(currentLevel)                         &
          &                  %faces                                        &
          &                  %faces(iDir)                                  &
          &                  %sendBuffer_flux(iFace),                      &
          & state        = facedata_list(currentLevel)%faceFlux(iDir)%dat, &
          & message_flag = currentLevel,                                   &
          & comm         = general%proc%comm                               )

      end do
    end do
    call tem_stopTimer( timerHandle = atl_timerHandles%commState )

    ! order of test function loops for the different directions
    ! The one indicates the current normal direction, while the
    ! other two refer to the tangentials.
    ! They have to be kept in order (2 before 3) in order to
    ! match the surface modes with the right volume modes.
    dirVec(:,1) = [ 1,2,3 ]
    dirVec(:,2) = [ 2,1,3 ]
    dirVec(:,3) = [ 2,3,1 ]

    call tem_startTimer( timerHandle = atl_timerHandles%physFlux )
    call modg_compute_project_physFlux(                                 &
      & mesh             = mesh_list(currentLevel),                     &
      & equation         = equation,                                    &
      & kerneldata       = kerneldata_list(currentLevel),               &
      & statedata        = statedata_list(currentLevel),                &
      & scheme           = scheme_list,                                 &
      & dl_prod          = scheme_list(currentLevel)%dl_prod,           &
      & dirVec           = dirVec,                                      &
      & eval_phy         = eval,                                        &
      & usePenalization  = usePenalization,                             &
      & apply_pen        = apply,                                       &
      & minLevel         = minLevel,                                    &
      & currentLevel     = currentLevel,                                &
      & maxLevel         = maxLevel,                                    &
      & poly_proj        = poly_proj_list(poly_proj_pos(currentLevel)), &
      & penalizationdata = penalizationdata_list(currentLevel),         &
      & material         = material_list(currentLevel)                  )
    call tem_stopTimer( timerHandle = atl_timerHandles%physFlux )


    call tem_startTimer( timerHandle = atl_timerHandles%projectTestFunc )
    call atl_modg_project_NumFlux(                                      &
      & mesh             = mesh_list(currentLevel),                     &
      & equation         = equation,                                    &
      & kerneldata       = kerneldata_list(currentLevel),               &
      & facedata         = facedata_list(currentLevel),                 &
      & scheme           = scheme_list(currentLevel)%modg,              &
      & dl_prod          = scheme_list(currentLevel)%dl_prod,           &
      & dl_prodDiff      = scheme_list(currentLevel)%dl_prodDiff,       &
      & poly_proj        = poly_proj_list(poly_proj_pos(currentLevel)), &
      & usePenalization  = usePenalization,                             &
      & penalizationdata = penalizationdata_list(currentLevel),         &
      & dirvec           = dirvec                                       )

    call atl_modg_project_source(                       &
      & sourcedata    = source,                         &
      & nScalars      = equation%varSys%nScalars,       &
      & mesh          = mesh_list(currentLevel),        &
      & scheme        = scheme_list(currentLevel)%modg, &
      & kerneldata    = kerneldata_list(currentLevel),  &
      & currentLevel  = currentLevel                    )

    call tem_stopTimer( timerHandle = atl_timerHandles%projectTestFunc )

  end subroutine compute_rhs_cubes_modg
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This routine is used to initialize an array in an OpenMP PARALLEL region.
  !! Usually this is done using a WORKSHARE directive, but due to a bug in
  !! Intel 15 we cannot make use of WORKSHARE.
  !!
  !! This routine is specifically made to initialize the state_der array of the
  !! atl_kerneldata_type, which is a three-dimensional array.
  subroutine atl_initialize_state_der(state_der)
    ! -------------------------------------------------------------------- !
    !> The state derivates of the kerneldata.
    real(kind=rk), intent(inout) :: state_der(:,:,:)
    ! -------------------------------------------------------------------- !
    integer :: lb(3), ub(3), length, n, x, y, z
    ! -------------------------------------------------------------------- !

    !! First we get the extents for the different dimensions.
    lb = lbound(state_der)
    ub = ubound(state_der)
    !! With these dimension we then can calculate the number of elements in
    !! the array.
    length = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*(ub(3)-lb(3)+1)

    !! This number is then used in a collapsed loop to initialize the array
    !! elements.

    do n = 1, length
      x = modulo(n - 1, ub(1) - lb(1) + 1) + lb(1)
      y = modulo((n - 1) / (ub(1) - lb(1) + 1), ub(2) - lb(2) + 1) + lb(2)
      z = ((n - 1) / ((ub(1) - lb(1) + 1) * (ub(2) - lb(2) + 1))) + lb(3)
      state_der(x, y, z) = 0.0_rk
    end do


  end subroutine atl_initialize_state_der
  ! ************************************************************************ !


  ! ************************************************************************ !
  !>TODO NA - Move this routine to the atl_modg_kernel_module
  subroutine modg_compute_project_physFlux( mesh, equation,kerneldata,   &
    & statedata, dirVec, poly_proj, material,scheme, dl_prod, apply_pen, &
    & penalizationdata, minLevel, currentLevel, maxLevel, eval_phy,      &
    & usePenalization                                                    )
    ! -------------------------------------------------------------------- !
    use atl_cube_elem_module,               only: atl_cube_elem_type
    use atl_equation_module,                only: atl_equations_type
    use atl_kerneldata_module,              only: atl_kerneldata_type, &
      &                                           atl_statedata_type
    use atl_scheme_module,                  only: atl_scheme_type
    use ply_poly_project_module,            only: ply_poly_project_type, &
      &                                           assignment(=)
    use atl_project_physflux_module,        only: &
      & atl_modg_project_PhysFlux_testFunc
    use atl_modg_LoclinEuler_kernel_module, only: &
      & atl_modg_LoclinEuler_physFlux
    use atl_modg_euler_kernel_module, only: &
      & atl_modg_euler_physFlux_const
    use atl_materialPrp_module,             only: atl_material_type
    use ply_oversample_module,              only: ply_convert2oversample,   &
      &                                           ply_convertFromoversample
    use ply_poly_project_module,            only: ply_poly_project_type, &
      &                                           assignment(=),         &
      &                                           ply_poly_project_n2m,  &
      &                                           ply_poly_project_m2n
    use atl_penalization_module,            only: atl_penalizationData_type
    use ply_leg_diff_module,                only: ply_calcDiff_leg
    use atl_physFlux_module,                only: atl_physFlux_pointer_type, &
      &                                           physflux_interface,        &
      &                                           atl_penalization_pointer_type
    use atl_physFluxEuler_module,           only: atl_physFluxEuler
    ! -------------------------------------------------------------------- !
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The representation on the face + representation of the flux.
    type(atl_statedata_type), intent(inout) :: statedata
    integer, intent(in) :: minLevel
    integer, intent(in) :: currentLevel
    integer, intent(in) :: maxLevel
    !> The parameters of the MODG scheme
    type(atl_scheme_type), intent(inout) :: scheme(minLevel:maxLevel)
    !> stored scalar products of the testfunction and anstaz function
    real(kind=rk), intent(in) :: dl_prod(:,:)
    !> vector for direction indicators
    integer, intent(in) :: dirVec(3,3)
    !> Material parameters (mu, epsilon) for all elements
    type(atl_material_type), intent(inout) :: material
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
    type(atl_penalizationData_type), intent(inout) :: penalizationdata
    type(atl_physFlux_pointer_type) :: eval_phy(2)
    type(atl_penalization_pointer_type) :: apply_pen(2)
    !> Flag indicating whether to apply the penalization or not.
    !!
    !! When a implicit scheme is used to integrate the penalized parts, this
    !! can be used to switch it off here.
    logical, intent(in) :: usePenalization
    ! -------------------------------------------------------------------- !
    !> The direction
    integer :: iDir
    !> Indicates whether the element has contant or variable material
    integer :: cov
    !> The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:), modalCoeffs_gradient(:,:,:)
    !> Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:), pointVal_gradient(:,:,:)
    !> The nodal representation of the physical flux along the 3 spatial
    !! directions.
    real(kind=rk), allocatable :: nodal_res(:,:)
    real(kind=rk), allocatable :: tmp_state_der(:,:)
    integer :: nquadpoints, oversamp_dofs, iElem, ndofs, elems, elemPos
    integer :: rot(5)
    logical :: use_linear_flux
    logical :: use_inviscid_flux
    procedure(physFlux_interface), pointer :: physFlux => null()
    ! -------------------------------------------------------------------- !
    ! get correct amount of quadrature points and degree due to projection
    ! method. oversamp_dof and oversamp_degree is used for the oversampling
    ! loop
    nquadpoints = poly_proj%body_3D%nquadpoints
    oversamp_dofs = poly_proj%body_3D%oversamp_dofs
    nDofs = poly_proj%body_3D%ndofs


    allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars) )
    allocate( pointVal(nQuadPoints,equation%varSys%nScalars) )
    allocate( nodal_Res(nQuadPoints, equation%varSys%nScalars) )
    allocate( tmp_state_der(nDofs, equation%varSys%nScalars) )
    ! Temp arrays for the Navierstokes. Should be removed
    allocate( modalCoeffs_gradient(oversamp_dofs,equation%varSys%nScalars,3))
    allocate( pointVal_gradient(nQuadPoints,equation%varSys%nScalars,3) )


!ICE Intel 15 !!!
!ICE!    kerneldata%state_der = 0.0_rk
    ! Intel 15 workaround:
    call tem_startTimer( timerHandle = atl_timerHandles%pF_initStateDer )
    call atl_initialize_state_der(kerneldata%state_der)
    call tem_stopTimer( timerHandle = atl_timerHandles%pF_initStateDer )

    ! Loop over the elements with material parameter
    ! cov = 1 : Constant Material Parameter
    ! cov = 2 : NonConstant Material Parameter
    do cov = 1,2

      use_linear_flux = .false.
      use_inviscid_flux = .false.

      ! begin the element loop for these elements
      do elems = 1, material%material_desc%computeElems(cov)%nElems
        iElem = material%material_desc%computeElems(cov)%totElemIndices(elems)

        elemPos = mesh%descriptor%PntTID(iElem)
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )

        physFlux => eval_phy(cov)%physflux

        ! In some cases we need the nodal values for the physflux evaluation!
        ! NOTE: when penalization is added to the equation system, it would also
        ! require nodal values
        call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv)
        select case(trim(equation%eq_kind))
        case('euler','navier_stokes','filtered_navier_stokes')
          if ( (cov == atl_ConstMatIdx) ) then
            ! Check whether it is allowed to use the linearized flux function
            ! to avoid projections between modal and nodal space.
            ! But we do this only for elements with constant material, in
            ! non-constant material elements, we always use the nonlinear
            ! flux computation.
            use_linear_flux = equation%euler%linear(      &
              & mean      = statedata%state(iElem,1,:),   &
              & deviation = kerneldata%deviation(iElem,:) )
            if (equation%nDerivatives == 1) then
              use_inviscid_flux = equation%navierstokes%inviscous( &
                & mean      = statedata%state(iElem,1,:),          &
                & deviation = kerneldata%deviation(iElem,:),       &
                & grad      = kerneldata%maxgrad(iElem,:)          )
              use_inviscid_flux = use_inviscid_flux         &
                &                .or. material%material_dat &
                &                             %mode_reducable(ielem)
            end if
          end if

          linearity: if (use_linear_flux) then
            ! It's deemed sufficient to compute the linearized equations
            ! locally.
            physflux => atl_modg_LoclinEuler_physFlux
          else linearity
            ! get the modal coefficients of the current cell (for all variables
            ! of the Euler equation, therefore we use ":" for the third index).
            ! ATTENTION: have to be duplicated as the FPT is modifying the
            ! input vector.
            modalCoeffs = 0.0
            ! --> modal space
            if (equation%euler%ensure_positivity) then
              call ply_convert2oversample(                        &
                & state             = statedata%state(iElem,:,:), &
                & poly_proj         = poly_proj,                  &
                & nDim              = 3,                          &
                & modalCoeffs       = modalCoeffs,                &
                & ensure_positivity = [.true.,                    &
                &                      .false., .false., .false., &
                &                      .true.  ]                  )
            else
              call ply_convert2oversample(                  &
                & state       = statedata%state(iElem,:,:), &
                & poly_proj   = poly_proj,                  &
                & nDim        = 3,                          &
                & modalCoeffs = modalCoeffs                 )
            end if
            ! --> oversamp modal space

            ! Now, we transform the modal representation of this element to
            ! nodal space by making use of fast polynomial transformations (FPT)
            call ply_poly_project_m2n( me         = poly_proj,                &
              &                        dim        = 3 ,                       &
              &                        nVars      = equation%varSys%nScalars, &
              &                        nodal_data = pointVal,                 &
              &                        modal_data = modalCoeffs               )
            ! --> oversamp nodal space

            ! for the navier stokes equation
            ! calculates the nodal gradient of modal representation
            navier: if (equation%nDerivatives == 1) then
              inviscosity: if (use_inviscid_flux) then
                ! It's deemed sufficient to compute the inviscid fluxes.
                ! Falling back to the Euler flux computation.
                physFlux => atl_modg_euler_physFlux_const
              else inviscosity
                call ply_convert2oversample(                  &
                  & state       = statedata%state(iElem,:,:), &
                  & poly_proj   = poly_proj,                  &
                  & nDim        = 3,                          &
                  & modalCoeffs = modalCoeffs                 )

                ! Calculate the gradient of modal Coeffs
                call ply_calcDiff_leg(                            &
                  &    legCoeffs     = modalCoeffs,               &
                  &    legCoeffsDiff = modalCoeffs_gradient,      &
                  &    maxPolyDegree = poly_proj%oversamp_degree, &
                  &    nVars         = equation%varSys%nScalars,  &
                  &    elemLength    = mesh%length                )

                ! Now, transform modal gradient representation of this element
                ! to nodal space by making use of fast polynomial
                ! transformations (FPT)
                do iDir = 1, 3
                  call ply_poly_project_m2n(                      &
                    & me         = poly_proj,                     &
                    & dim        = 3,                             &
                    & nVars      = equation%varSys%nScalars,      &
                    & nodal_data = pointVal_gradient(:,:,iDir),   &
                    & modal_data = modalCoeffs_gradient(:,:,iDir) )
                end do
              end if inviscosity
            end if navier

          end if linearity

        end select
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv)

        ! Perform Penalization if active and not computed elsewhere (imex)
        call tem_startTimer( timerHandle = atl_timerHandles%pF_pen )
        if (penalizationData%isActive .and. usePenalization) then

          call apply_pen(cov)%pen( equation         = equation,         &
            &                      poly_proj        = poly_proj,        &
            &                      nodal_data       = pointVal,         &
            &                      scheme_min       = scheme(minlevel), &
            &                      penalizationData = penalizationData, &
            &                      iElem            = elems,            &
            &                      material         = material          )

        end if
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_pen )
        do iDir = 1,3
          ! Evaluate the physical flux
          ! If the element is covered completely by the material, than just
          ! the first mode is considered for the flux computation.
          ! This flag can only be true for Euler equations (and Navier-Stokes,
          ! but for Navier-Stokes the velocity gradients are assumed 0, and
          ! thus, the viscous fluxes are ignored).
          if (  material%material_dat%mode_reducable(iElem) ) then

            rot = equation%varRotation(iDir)%varTransformIndices(1:5)
            call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
            tmp_state_der = 0.0_rk
            tmp_state_der(1,rot) = atl_physFluxEuler(                        &
              & state        = statedata%state(iElem,1,rot),                 &
              & isenCoeff    = equation%euler%isen_coef,                     &
              & penalty_char = material%material_dat%elemMaterialData(cov)   &
              &                                     %materialDat(elems,1,1), &
              & U_o          = material%material_dat                         &
              &                        %elemMaterialData(cov)                &
              &                        %materialDat(elems,1,iDir+1),         &
              & porosity     = equation%euler%porosity                       )
            call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
          else
            call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
            call physflux(                                     &
              & equation         = equation ,                  &
              & res              = tmp_state_der,              &
              & state            = statedata%state(iElem,:,:), &
              & iElem            = elems,                      &
              & iDir             = iDir,                       &
              & penalizationData = penalizationData,           &
              & poly_proj        = poly_proj,                  &
              & material         = material,                   &
              & nodal_data       = pointVal,                   &
              & nodal_gradData   = pointVal_gradient,          &
              & nodal_res        = nodal_res,                  &
              & elemLength       = mesh%length,                &
              & scheme_min       = scheme(minlevel),           &
              & scheme_current   = scheme(currentLevel)        )
            call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
            call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv )
            select case(trim(equation%eq_kind))
            case('euler','navier_stokes','filtered_navier_stokes')
              if (.not. use_linear_flux) then
                tmp_state_der = 0.0_rk
                ! Transform the nodal physical flux back to modal space
                ! --> oversamp modal space
                call ply_poly_project_n2m(                 &
                  & me         = poly_proj,                &
                  & dim        = 3 ,                       &
                  & nVars      = equation%varSys%nScalars, &
                  & nodal_data = nodal_res,                &
                  & modal_data = modalCoeffs               )
                ! --> oversamp modal space
                call ply_convertFromOversample( modalCoeffs = modalCoeffs,  &
                  &                             poly_proj   = poly_proj,    &
                  &                             nDim        = 3,            &
                  &                             state       = tmp_state_der )
              end if

            end select
            ! --> modal space
            call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv )
          end if

          call tem_startTimer(                                  &
            & timerHandle = atl_timerHandles%pF_projectTestFunc )
          call atl_modg_project_PhysFlux_testFunc(    &
            & mesh       = mesh,                      &
            & equation   = equation,                  &
            & kerneldata = kerneldata,                &
            & iDir       = iDir,                      &
            & scheme     = scheme(currentLevel)%modg, &
            & dl_prod    = dl_prod,                   &
            & dirvec     = dirvec,                    &
            & iElem      = iElem,                     &
            & state_der  = tmp_state_der              )
          call tem_stopTimer(                                   &
            & timerHandle = atl_timerHandles%pF_projectTestFunc )

        end do ! Direction Loop
        call tem_stopTimer( me = atl_elemTimers , timerHandle = elemPos )
      end do ! Element loop
    end do   ! end loop over materials

  end subroutine modg_compute_project_physFlux
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Computes the right hand side for cubical elements and 2D MODG scheme.
  subroutine compute_rhs_cubes_modg_2d( minLevel, maxLevel, currentLevel, &
    & mesh_list, kerneldata_list, statedata_list, facedata_list, source,  &
    & penalizationdata_list, scheme_list, poly_proj_pos, poly_proj_list,  &
    & equation, material_list, general, computePenalization               )
    ! -------------------------------------------------------------------- !
    use atl_cube_elem_module, only: &
      & atl_cube_elem_type
    use atl_source_types_module, only: &
      & atl_source_type
    use atl_scheme_module, only: &
      & atl_scheme_type
    use atl_equation_module, only: &
      & atl_equations_type
    use atl_modg_2d_kernel_module, only: &
      & atl_modg_2d_project_NumFlux,     &
      & atl_modg_2d_project_source
    use atl_modg_2d_multilevel_module, only: &
      & atl_modg_2d_fineToCoarseFace
    use atl_modg_2d_maxwell_kernel_module, only: &
      & atl_modg_maxwell_2d_numflux,             &
      & atl_modg_maxwell_2d_physFlux_const,      &
      & atl_modg_maxwell_2d_physFlux_Nonconst,   &
      & atl_modg_maxwell_2d_penalization_Const,  &
      & atl_modg_maxwell_2d_penalization_NonConst
    use atl_modg_2d_lineareuler_kernel_module,  only: &
      & atl_modg_2d_lineareuler_numflux,              &
      & atl_modg_2d_lineareuler_physflux
    use atl_modg_2d_euler_kernel_module, only: &
      & atl_modg_2d_euler_numflux,             &
      & atl_modg_2d_euler_physFlux_const,      &
      & atl_modg_2d_euler_physFlux_Nonconst,   &
      & atl_modg_2d_euler_penalization_const,  &
      & atl_modg_2d_euler_penalization_Nonconst
    use atl_modg_2d_navierstokes_kernel_module, only: &
      & atl_modg_2d_navierstokes_numflux,             &
      & atl_modg_2d_navierstokes_physFlux_const,      &
      & atl_modg_2d_navierstokes_physFlux_Nonconst,   &
      & atl_modg_2d_navierstokes_penalization_const,  &
      & atl_modg_2d_navierstokes_penalization_Nonconst
    use atl_modg_2d_filNvrStk_kernel_module, only: &
      & atl_modg_2d_filNvrStk_numflux,             &
      & atl_modg_2d_filNvrStk_physFlux_const,      &
      & atl_modg_2d_filNvrStk_physFlux_Nonconst
    use atl_elemental_time_integration_module, only: &
      & atl_timestep_type
    use atl_facedata_module, only: &
      & atl_facedata_type
    use atl_kerneldata_module, only: &
      & atl_kerneldata_type,         &
      & atl_statedata_type
    use ply_poly_project_module, only: &
      & ply_poly_project_type, &
      & assignment(=)
    use atl_materialPrp_module, only: &
      & atl_material_type
    use ply_dynArray_project_module, only: &
      & dyn_projectionArray_type
    use atl_penalization_module, only: &
      & atl_penalizationData_type
    use atl_modg_2d_heat_kernel_module, only: &
      & atl_modg_2d_heat_numflux,             &
      & atl_modg_2d_heat_physFlux
    use atl_modg_2d_acoustic_kernel_module, only: &
      & atl_modg_2d_acoustic_numflux, &
      & atl_modg_2d_acoustic_physFlux
    use atl_physFlux_module, only: &
      & atl_physFlux_pointer_type, &
      & atl_penalization_pointer_type
    use treelmesh_module, only: &
      & treelmesh_type
    ! -------------------------------------------------------------------- !
    !> The minimum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The level of the mesh you are computing the rhs for.
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type),  intent(inout) :: mesh_list(minLevel:maxLevel)
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: &
      & kerneldata_list(minLevel:maxLevel)
    !> List of facedatas. For each level we have one
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list(minLevel:maxLevel)
    !> Levelwise list of sources
    type(atl_source_type), intent(inout) :: source
    !> Levelwise list of penalization data
    type(atl_penalizationData_type), intent(inout) :: &
      & penalizationdata_list(minLevel:maxLevel)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> List of levelwise position of  projection method in unique
    !! projection_list
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> The equation you are operating with.
    type(atl_equations_type),intent(in) :: equation
    !> Information about the material parameters of the equation.
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxlevel)
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> Flag to indicate whether penalization terms should be computed or not.
    !!
    !! This is used to switch off the application of penalizing terms from
    !! the convective computation and instead compute it in an implicit
    !! update, see the atl_imexrk_module.
    !! Default is .true.!
    logical, intent(in), optional :: computePenalization
    ! -------------------------------------------------------------------- !
    integer :: iDir, iFace
    type(atl_physFlux_pointer_type) :: eval(2)
    type(atl_penalization_pointer_type) :: apply(2)
    logical :: usePenalization
    ! -------------------------------------------------------------------- !

    usePenalization = .true.
    if (present(computePenalization)) usePenalization = computePenalization

    ! calculate rhs of the PDE from source terms
    call atl_update_sourcedata(                                            &
      & equation     = equation,                                           &
      & time         = statedata_list(currentLevel)%local_time,            &
      & mesh         = mesh_list(currentLevel),                            &
      & poly_proj    = poly_proj_list(source%poly_proj_pos(currentLevel)), &
      & currentLevel = currentLevel,                                       &
      & state        = statedata_list(currentLevel)%state,                 &
      & material     = material_list(currentLevel),                        &
      & source       = source,                                             &
      & scheme       = scheme_list(currentLevel)                           )

    ! Interpolate fluxes from finer level to current level. Please keep in mind
    ! that the timestepping routine is called recursively for the levels.
    ! Therefore the compute rhs subroutine (i.e. the current routine) is
    ! executed on the finer levels before it is executed on the coarser level.
    ! So, the flux on the finer level is already available when a coarser level
    ! enters this routine.
    if( currentLevel .lt. maxLevel ) then
      call atl_modg_2d_fineToCoarseFace(                                    &
        & minLevel     = minLevel,                                          &
        & maxLevel     = maxLevel,                                          &
        & currentLevel = currentLevel,                                      &
        & mesh         = mesh_list,                                         &
        & facedata     = facedata_list,                                     &
        & scheme       = scheme_list,                                       &
        & nScalars     = equation%varSys%nScalars*(equation%nDerivatives+1) )
    end if

    ! Compute all the numerical fluxes, by the modal representations of the
    ! states on the faces. The modal representation on the faces is provided by
    ! the preprocess step.
    call tem_startTimer( timerHandle = atl_timerHandles%numFlux )
    select case(equation%eq_kind)
    case('maxwell_2d')
      call atl_modg_maxwell_2d_numFlux(                            &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg_2d,           &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux  => atl_modg_maxwell_2d_physflux_const
      eval(2)%physFlux  => atl_modg_maxwell_2d_physflux_NonConst

      !Set up the pointer for the penalization routines
      apply(1)%pen => atl_modg_maxwell_2d_penalization_Const
      apply(2)%pen => atl_modg_maxwell_2d_penalization_NonConst

    case('acoustic_2d')
      call atl_modg_2d_acoustic_numFlux(                &
        & mesh      = mesh_list(currentLevel),          &
        & equation  = equation,                         &
        & facedata  = facedata_list(currentLevel),      &
        & scheme    = scheme_list(currentLevel)%modg_2d )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_2d_acoustic_physflux

    case('lineareuler_2d')
      call atl_modg_2d_linearEuler_numFlux(            &
        &   mesh = mesh_list(currentLevel),            &
        &   equation = equation,                       &
        &   facedata = facedata_list(currentLevel),    &
        &   scheme = scheme_list(currentLevel)%modg_2d )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_2d_linearEuler_physFlux

    case('euler_2d')
      call atl_modg_2d_euler_numFlux(                              &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_2d_euler_physflux_const
      eval(2)%physFlux => atl_modg_2d_euler_physflux_NonConst

      !Set up the pointer for the penalization routines
      apply(1)%pen => atl_modg_2d_euler_penalization_Const
      apply(2)%pen => atl_modg_2d_euler_penalization_NonConst

    case('navier_stokes_2d')
      call atl_modg_2d_navierstokes_numFlux(                       &
        & mesh      = mesh_list(currentLevel),                     &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg_2d,           &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_2d_navierstokes_physflux_const
      eval(2)%physFlux => atl_modg_2d_navierstokes_physflux_NonConst

      !Set up the pointer for the penalization routines
      apply(1)%pen => atl_modg_2d_navierstokes_penalization_Const
      apply(2)%pen => atl_modg_2d_navierstokes_penalization_NonConst

    case('filtered_navier_stokes_2d')
      call atl_modg_2d_filNvrStk_numFlux(                          &
        & mesh      = mesh_list(currentLevel),                     &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & scheme    = scheme_list(currentLevel)%modg_2d,           &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & material  = material_list(currentLevel)                  )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux => atl_modg_2d_filNvrStk_physflux_const
      eval(2)%physFlux => atl_modg_2d_filNvrStk_physflux_NonConst

    case('heat_2d')
      call atl_modg_2d_heat_numFlux(                               &
        & mesh      = mesh_list(currentLevel),                     &
        & equation  = equation,                                    &
        & facedata  = facedata_list(currentLevel),                 &
        & poly_proj = poly_proj_list(poly_proj_pos(currentLevel)), &
        & scheme    = scheme_list(currentLevel)%modg_2d            )

      !Set up the pointer for the physical fluxes
      eval(1)%physFlux=> atl_modg_2d_heat_physflux

    case default
      call tem_abort( 'ERROR in compute_rhs_cubes: modg is not supporting' &
        & // ' this PDE (num flux calculation), stopping...'               )
    end select
    call tem_stopTimer( timerHandle = atl_timerHandles%numFlux )

    ! Exchange all the fluxes on the faces (for 2 directions, x and y dir)
    ! We send data about the fluxes at those faces where we received
    ! information about the states (in preprocess step). Therefore,
    ! we use the send buffer for receiving and the receive buffer for
    ! sending.
    call tem_startTimer( timerHandle = atl_timerHandles%commState )
    do iDir = 1,2
      do iFace = 1,2
        call general%commPattern%exchange_real(                            &
          & send         = mesh_list(currentLevel)                         &
          &                  %faces                                        &
          &                  %faces(iDir)                                  &
          &                  %recvBuffer_flux(iFace),                      &
          & recv         = mesh_list(currentLevel)                         &
          &                  %faces                                        &
          &                  %faces(iDir)                                  &
          &                  %sendBuffer_flux(iFace),                      &
          & state        = facedata_list(currentLevel)%faceFlux(iDir)%dat, &
          & message_flag = currentLevel,                                   &
          & comm         = general%proc%comm                               )

      end do
    end do
    call tem_stopTimer( timerHandle = atl_timerHandles%commState )

    call tem_startTimer( timerHandle = atl_timerHandles%physFlux )
    call modg_2d_compute_project_physFlux(                              &
      & mesh             = mesh_list(currentLevel),                     &
      & equation         = equation,                                    &
      & kerneldata       = kerneldata_list(currentLevel),               &
      & statedata        = statedata_list(currentLevel),                &
      & scheme           = scheme_list,                                 &
      & poly_proj        = poly_proj_list(poly_proj_pos(currentLevel)), &
      & dl_prod          = scheme_list(currentLevel)%dl_prod,           &
      & penalizationdata = penalizationdata_list(currentLevel),         &
      & material         = material_list(currentLevel),                 &
      & minLevel         = minLevel,                                    &
      & maxLevel         = maxLevel,                                    &
      & currentLevel     = currentLevel,                                &
      & eval_phy         = eval,                                        &
      & usePenalization  = usePenalization,                             &
      & apply_pen        = apply                                        )
    call tem_stopTimer( timerHandle = atl_timerHandles%physFlux )

    ! now project the numerical fluxes and the source terms to the test func
    call tem_startTimer( timerHandle = atl_timerHandles%projectTestFunc )
    call atl_modg_2d_project_NumFlux(                                   &
      & mesh             = mesh_list(currentLevel),                     &
      & equation         = equation,                                    &
      & kerneldata       = kerneldata_list(currentLevel),               &
      & facedata         = facedata_list(currentLevel),                 &
      & penalizationdata = penalizationdata_list(currentLevel),         &
      & usePenalization  = usePenalization,                             &
      & scheme           = scheme_list(currentLevel)%modg_2d,           &
      & poly_proj        = poly_proj_list(poly_proj_pos(currentLevel)), &
      & dl_prod          = scheme_list(currentLevel)%dl_prod,           &
      & dl_prodDiff      = scheme_list(currentLevel)%dl_prodDiff        )

    call atl_modg_2d_project_source(                       &
      & sourcedata    = source,                            &
      & nScalars      = equation%varSys%nScalars,          &
      & mesh          = mesh_list(currentLevel),           &
      & scheme        = scheme_list(currentLevel)%modg_2d, &
      & kerneldata    = kerneldata_list(currentLevel),     &
      & currentLevel  = currentLevel                       )

    call tem_stopTimer( timerHandle = atl_timerHandles%projectTestFunc )

  end subroutine compute_rhs_cubes_modg_2d
  ! *********************************************************************** !


  ! *********************************************************************** !
  !> This subroutine computes the physical fluxes for various equation system
  ! and projects it onto the test function directionwise
  subroutine modg_2d_compute_project_physFlux( mesh, equation, kerneldata,   &
    & statedata, scheme, poly_proj, dl_prod, penalizationdata, material,     &
    & minLevel, maxLevel, currentLevel, eval_phy, apply_pen, usePenalization )

    use atl_cube_elem_module,               only: atl_cube_elem_type
    use atl_scheme_module,                  only: atl_scheme_type
    use atl_equation_module,                only: atl_equations_type
    use atl_project_physflux_module,        only: &
      & atl_modg_2d_project_physFlux_testFunc
    use atl_kerneldata_module,              only: atl_kerneldata_type, &
      &                                           atl_statedata_type
    use ply_poly_project_module,            only: ply_poly_project_type, &
      &                                           ply_poly_project_n2m,  &
      &                                           assignment(=),         &
      &                                           ply_poly_project_m2n
    use atl_materialPrp_module,             only: atl_material_type
    use ply_oversample_module,              only: ply_convert2oversample,   &
      &                                           ply_convertFromoversample
    use atl_penalization_module,            only: atl_penalizationData_type
    use atl_materialPrp_module,             only: atl_material_type
    use atl_modg_2d_heat_kernel_module,     only: atl_modg_2d_heat_physFlux
    use ply_leg_diff_module,                only: ply_calcDiff_leg_2d
    use atl_modg_2d_euler_kernel_module, only: &
      & atl_modg_2d_euler_physFlux_const
    use atl_modg_2d_acoustic_kernel_module, only: &
      & atl_modg_2d_acoustic_physFlux
    use atl_modg_2d_loclineuler_kernel_module, only: &
      &   atl_modg_2d_loclineuler_physflux
    use atl_physFlux_module,                only: &
      &                                     atl_physFlux_pointer_type,  &
      &                                     physFlux_interface,  &
      &                                     atl_penalization_pointer_type
    use atl_physFluxEuler_2d_module,        only: atl_physFluxEuler_2d

    !> the levels of the geometry
    integer, intent(in) :: minLevel, maxLevel, currentLevel
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The representation on the face + representation of the flux.
    type(atl_statedata_type), intent(in) :: statedata
    !> The parameters of the MODG scheme
    type(atl_scheme_type), intent(inout) :: scheme(minLevel:maxLevel)
    !> stored scalar products of the testfunction and anstaz function
    real(kind=rk), intent(in) :: &
      & dl_prod(2, scheme(currentLevel)%modg_2d%maxPolyDegree+1)
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
    type(atl_penalizationData_type), intent(inout) :: penalizationdata
    !> Material parameters (mu, epsilon) for all elements
    type(atl_material_type), intent(inout) :: material
    type(atl_physFlux_pointer_type) :: eval_phy(2)
    type(atl_penalization_pointer_type) :: apply_pen(2)
    !> Flag indicating whether to apply the penalization or not.
    !!
    !! When a implicit scheme is used to integrate the penalized parts, this
    !! can be used to switch it off here.
    logical, intent(in) :: usePenalization
    ! -------------------------------------------------------------------- !
    !> The direction
    integer :: iDir, matType
    !> The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    !> Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:)
    !> The nodal representation of the physical flux along the 3 spatial
    !! directions.
    real(kind=rk), allocatable :: nodal_res(:,:)
    real(kind=rk), allocatable :: tmp_state_der(:,:)
    integer :: nquadpoints, oversamp_dofs, iElem, ndofs, elems, elemPos
    logical :: use_linear_flux
    logical :: use_inviscid_flux
    procedure(physFlux_interface), pointer :: physFlux => null()
    ! -------------------------------------------------------------------- !
    integer :: rot(4)
    real(kind=rk), allocatable :: modalCoeffs_gradient(:,:,:)
    real(kind=rk), allocatable :: pointVal_gradient(:,:,:)
    ! -------------------------------------------------------------------- !

    oversamp_dofs = poly_proj%body_2D%oversamp_dofs
    nquadpoints = poly_proj%body_2D%nquadpoints
    nDofs = poly_proj%body_2D%ndofs

    allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars) )
    allocate( pointVal(nQuadPoints,equation%varSys%nScalars) )
    allocate( nodal_Res(nQuadPoints, equation%varSys%nScalars) )
    allocate( tmp_state_der(nDofs, equation%varSys%nScalars) )
    ! Temp arrays a.t.m. only useful for the Navierstokes and rans
    allocate( modalCoeffs_gradient(oversamp_dofs,equation%varSys%nScalars,2))
    allocate( pointVal_gradient(nQuadPoints,equation%varSys%nScalars,2) )


!ICE Intel 15 !!!
!ICE!    kerneldata%state_der = 0.0_rk
    ! Intel 15 workaround:
    call tem_startTimer( timerHandle = atl_timerHandles%pF_initStateDer )
    call atl_initialize_state_der(kerneldata%state_der)
    call tem_stopTimer( timerHandle = atl_timerHandles%pF_initStateDer )


    ! Loop over the elements with material parameter
    do matType = atl_ConstMatIdx, atl_VarMatIdx

      use_linear_flux = .false.
      use_inviscid_flux = .false.

      ! begin the element loop for these elements
      do elems = 1, material%material_desc%computeElems(matType)%nElems
        iElem = material%material_desc%computeElems(matType) &
          &                           %totElemIndices(elems)

        physFlux => eval_phy(matType)%physflux

        elemPos = mesh%descriptor%PntTID(iElem)
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )

        call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv)
        select case(trim(equation%eq_kind))
        case('euler_2d', 'navier_stokes_2d', 'filtered_navier_stokes_2d')
          if ( (matType == atl_ConstMatIdx) ) then
            ! Check whether it is allowed to use the linearized flux function
            ! to avoid projections between modal and nodal space.
            ! But we do this only for elements with constant material, in
            ! non-constant material elements, we always use the nonlinear
            ! flux computation.
            use_linear_flux = equation%euler%linear(                      &
              &                 mean      = statedata%state(iElem,1,:),   &
              &                 deviation = kerneldata%deviation(iElem,:) )
            if (equation%nDerivatives == 1) then
              use_inviscid_flux = equation%navierstokes%inviscous( &
                & mean      = statedata%state(iElem,1,:),          &
                & deviation = kerneldata%deviation(iElem,:),       &
                & grad      = kerneldata%maxgrad(iElem,:)          )
              use_inviscid_flux = use_inviscid_flux         &
                &                .or. material%material_dat &
                &                             %mode_reducable(ielem)
            end if
          end if

          linearity: if (use_linear_flux) then
            ! It's deemed sufficient to compute the linearized equations
            ! locally.
            physFlux => atl_modg_2d_LoclinEuler_physFlux
          else linearity
            ! get the modal coefficients of the current cell (for all variables
            ! of the Euler equation, therefore we use ":" for the third index).
            ! ATTENTION: have to be duplicated as the FPT is modifying
            ! the input vector.
            modalCoeffs = 0.0
            ! --> modal space
            if (equation%euler%ensure_positivity) then
              call ply_convert2oversample(                  &
                & state       = statedata%state(iElem,:,:), &
                & poly_proj   = poly_proj,                  &
                & nDim        = 2,                          &
                & modalCoeffs = modalCoeffs,                &
                & ensure_positivity = [.true.,              &
                &                      .false., .false.,    &
                &                      .true.]              )
            else
              call ply_convert2oversample(                  &
                & state       = statedata%state(iElem,:,:), &
                & poly_proj   = poly_proj,                  &
                & nDim        = 2,                          &
                & modalCoeffs = modalCoeffs                 )
            end if
            ! --> oversamp modal space

            ! Now, we transform the modal representation of this element to nodal
            ! space by making use of fast polynomial transformations (FPT)
            call ply_poly_project_m2n( me         = poly_proj,                &
             &                         dim        = 2,                        &
             &                         nVars      = equation%varSys%nScalars, &
             &                         nodal_data = pointVal,                 &
             &                         modal_data = modalCoeffs               )

            ! --> oversamp nodal space
            ! for the navier stokes equation
            ! calculates the nodal gradient of modal representation
            navier: if (equation%nDerivatives ==1) then
              inviscosity: if (use_inviscid_flux) then
                ! It's deemed sufficient to compute the inviscid fluxes.
                ! Falling back to the Euler flux computation.
                physFlux => atl_modg_2d_euler_physFlux_const
              else inviscosity
                call ply_convert2oversample(                  &
                  & state       = statedata%state(iElem,:,:), &
                  & poly_proj   = poly_proj,                  &
                  & nDim        = 2,                          &
                  & modalCoeffs = modalCoeffs                 )

                ! Calculate the gradient of modal Coeffs
                call ply_calcDiff_leg_2d(                      &
                  & legCoeffs     = modalCoeffs,               &
                  & legCoeffsDiff = modalCoeffs_gradient,      &
                  & maxPolyDegree = poly_proj%oversamp_degree, &
                  & nVars         = equation%varSys%nScalars,  &
                  & elemLength    = mesh%length                )

                ! Now, transform modal gradient representation of this element to
                ! nodal space by making use of fast polynomial transformations
                ! (FPT)
                do iDir = 1, 2
                  call ply_poly_project_m2n(                      &
                    & me         = poly_proj,                     &
                    & dim        = 2,                             &
                    & nVars      = equation%varSys%nScalars,      &
                    & nodal_data = pointVal_gradient(:,:,iDir),   &
                    & modal_data = modalCoeffs_gradient(:,:,iDir) )
                end do
              end if inviscosity
            end if navier

          end if linearity

        case('maxwell_2d')
          if (penalizationData%isActive .and. usePenalization) then
            ! get the modal coefficients of the current cell (for all variables
            ! of the Euler equation, therefore we use ":" for the third index).
            ! ATTENTION: have to be duplicated as the FPT is modifying
            ! the input vector.
            modalCoeffs = 0.0
            ! --> modal space
            call ply_convert2oversample(                  &
              & state       = statedata%state(iElem,:,:), &
              & poly_proj   = poly_proj,                  &
              & nDim        = 2,                          &
              & modalCoeffs = modalCoeffs                 )
            ! --> oversamp modal space

            ! Now, we transform the modal representation of this element to nodal
            ! space by making use of fast polynomial transformations (FPT)
            call ply_poly_project_m2n( me         = poly_proj,                &
             &                         dim        = 2,                        &
             &                         nVars      = equation%varSys%nScalars, &
             &                         nodal_data = pointVal,                 &
             &                         modal_data = modalCoeffs               )

            ! --> oversamp nodal space
            ! for the navier stokes equation
            ! calculates the nodal gradient of modal representation
            if (equation%nDerivatives ==1) then
              ! it needs to be done again as the modalCoeffs above is modified
              ! when converting it to nodal values
              call ply_convert2oversample(                  &
                & state       = statedata%state(iElem,:,:), &
                & poly_proj   = poly_proj,                  &
                & nDim        = 2,                          &
                & modalCoeffs = modalCoeffs                 )

              ! Calculate the gradient of modal Coeffs
              call ply_calcDiff_leg_2d(                      &
                & legCoeffs     = modalCoeffs,               &
                & legCoeffsDiff = modalCoeffs_gradient,      &
                & maxPolyDegree = poly_proj%oversamp_degree, &
                & nVars         = equation%varSys%nScalars,  &
                & elemLength    = mesh%length                )

              ! Now, transform modal gradient representation of this element to
              ! nodal space by making use of fast polynomial transformations
              ! (FPT)
              do iDir = 1, 2
                call ply_poly_project_m2n(                      &
                  & me         = poly_proj,                     &
                  & dim        = 2,                             &
                  & nVars      = equation%varSys%nScalars,      &
                  & nodal_data = pointVal_gradient(:,:,iDir),   &
                  & modal_data = modalCoeffs_gradient(:,:,iDir) )
              end do
            end if
          end if
        end select
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv)


        ! Perform Penalization if active and not computed elsewhere (imex)
        call tem_startTimer( timerHandle = atl_timerHandles%pF_pen )
        if (penalizationData%isActive .and. usePenalization) then

          call apply_pen(matType)%pen( equation         = equation,         &
            &                          poly_proj        = poly_proj,        &
            &                          nodal_data       = pointVal,         &
            &                          scheme_min       = scheme(minlevel), &
            &                          penalizationData = penalizationData, &
            &                          iElem            = elems,            &
            &                          material         = material          )

        endif
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_pen )

        do iDir = 1,2

          if (use_linear_flux) then
            call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
            call atl_modg_2d_LoclinEuler_physFlux(                &
              & equation         = equation ,                  &
              & res              = tmp_state_der,              &
              & state            = statedata%state(iElem,:,:), &
              & iElem            = elems,                      &
              & iDir             = iDir,                       &
              & penalizationData = penalizationData,           &
              & poly_proj        = poly_proj,                  &
              & material         = material,                   &
              & nodal_data       = pointVal,                   &
              & nodal_gradData   = pointVal_gradient,          &
              & nodal_res        = nodal_res,                  &
              & elemLength       = mesh%length,                &
              & scheme_min       = scheme(minlevel),           &
              & scheme_current   = scheme(currentLevel)        )
            call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
          else
            !> If the element is covered completely by the material, than just the first
            !! mode is considered for the flux computation.
            if (  material%material_dat%mode_reducable(iElem) ) then
              rot = equation%varRotation(iDir)%varTransformIndices(1:4)

             call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
             tmp_state_der(:,:) = 0.0_rk
             tmp_state_der(1,rot) = atl_physFluxEuler_2d(            &
               & state        = statedata%state(iElem,1,rot),        &
               & isenCoeff    = equation%euler%isen_coef,            &
               & penalty_char = material%material_dat                &
               &                        %elemMaterialData(matType)   &
               &                        %materialDat(elems,1,1),     &
               & porosity     = equation%euler%porosity,             &
               & U_o          = material%material_dat                &
               &                        %elemMaterialData(matType)   &
               &                        %materialDat(elems,1,iDir+1) )
             call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
            else
              call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
              ! Evaluate the physical flux
              !! write(*,*) 'Elements with chi 0 using nonlinear flux '
              call physFlux(                                     &
                & equation         = equation,                   &
                & res              = tmp_state_der,              &
                & state            = statedata%state(iElem,:,:), &
                & iElem            = elems,                      &
                & iDir             = iDir,                       &
                & penalizationData = penalizationData,           &
                & poly_proj        = poly_proj,                  &
                & material         = material,                   &
                & nodal_data       = pointVal,                   &
                & nodal_gradData   = pointVal_gradient,          &
                & nodal_res        = nodal_res,                  &
                & elemLength       = mesh%length,                &
                & scheme_min       = scheme(minlevel),           &
                & scheme_current   = scheme(currentLevel)        )
              call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
              call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv )
              select case(trim(equation%eq_kind))
              case('euler_2d', 'navier_stokes_2d', 'filtered_navier_stokes_2d')

                ! Transform the nodal physical flux back to modal space
                ! --> oversamp modal space
                call ply_poly_project_n2m(                 &
                  & me         = poly_proj,                &
                  & dim        = 2,                        &
                  & nVars      = equation%varSys%nScalars, &
                  & nodal_data = nodal_res,                &
                  & modal_data = modalCoeffs               )

                ! --> oversamp modal space
                call ply_convertFromOversample( modalCoeffs = modalCoeffs,  &
                  &                             poly_proj   = poly_proj,    &
                  &                             nDim        = 2,            &
                  &                             state       = tmp_state_der )

              end select
              ! modal space
              call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv )
            end if
          end if

          call tem_startTimer(timerHandle = atl_timerHandles%pF_projectTestFunc)
          call atl_modg_2d_project_PhysFlux_testFunc(    &
            & mesh       = mesh,                         &
            & equation   = equation,                     &
            & kerneldata = kerneldata,                   &
            & iDir       = iDir,                         &
            & scheme     = scheme(currentLevel)%modg_2d, &
            & dl_prod    = dl_prod,                      &
            & iElem      = iElem,                        &
            & nDofs      = nDofs,                        &
            & state_data = tmp_state_der                 )
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_projectTestFunc)

        end do ! The direction loop
        call tem_stopTimer ( me = atl_elemTimers , timerHandle = elemPos  )
      end do  ! Loop over elements
    end do


  end subroutine modg_2d_compute_project_physFlux
  ! *********************************************************************** !


  ! *********************************************************************** !
  !> Computes the right hand side for cubical elements and 1D MODG scheme.
  subroutine compute_rhs_cubes_modg_1d(                                    &
    & minLevel, maxLevel, currentLevel, mesh_list, tree, kerneldata_list,  &
    & statedata_list, facedata_list, source, penalizationdata_list,        &
    & scheme_list, poly_proj_pos, poly_proj_list, equation, material_list, &
    & general, computePenalization                                         )
    use atl_cube_elem_module,             only: atl_cube_elem_type
    use atl_source_types_module,          only: atl_source_type
    use atl_scheme_module,                only: atl_scheme_type
    use atl_equation_module,              only: atl_equations_type
    use atl_materialPrp_module,           only: atl_material_type
    use atl_modg_1d_kernel_module,        only: atl_modg_1d_project_testFunc
    use atl_modg_1d_multilevel_module,    only: atl_modg_1d_fineToCoarseFace
    use atl_modg_1d_euler_kernel_module,  only: atl_modg_1d_euler_numflux
    use atl_modg_1d_LoclinEuler_kernel_module,                                 &
      &                                   only: atl_modg_1d_LoclinEuler_physFlux
    use atl_modg_1d_advection_kernel_module, only: &
      & atl_modg_1d_advection_numflux
    use atl_modg_1d_heat_kernel_module,   only: atl_modg_1d_heat_numflux
    use atl_elemental_time_integration_module,              &
      &                                   only: atl_timestep_type
    use atl_facedata_module,              only: atl_facedata_type
    use atl_kerneldata_module,            only: atl_kerneldata_type, &
      &                                         atl_statedata_type
    use atl_penalization_module,          only: atl_penalizationData_type
    use ply_poly_project_module,          only: ply_poly_project_type, &
      &                                         assignment(=)
    use ply_dynArray_project_module,      only: dyn_projectionArray_type
    use treelmesh_module,                 only: treelmesh_type
    ! -------------------------------------------------------------------- !
    !> The minimum level of the mesh.
    integer, intent(in) :: minLevel
    !> The maximum level of the mesh.
    integer, intent(in) :: maxLevel
    !> The level of the mesh you are computing the rhs for.
    integer, intent(in) :: currentLevel
    !> List of mesh parts. For each level we have one.
    type(atl_cube_elem_type), intent(inout) :: mesh_list(minLevel:maxLevel)
    !> treelm mesh
    type(treelmesh_type), intent(in) :: tree
    !> List of kerneldatas. For each level we have one
    type(atl_kerneldata_type), intent(inout) :: &
      & kerneldata_list(minLevel:maxLevel)
    !> List of facedatas. For each level we have one
    type(atl_facedata_type), intent(inout) :: facedata_list(minLevel:maxLevel)
    !> List of states you want to calc the rhs for. For each level we have one.
    type(atl_statedata_type), intent(inout) :: statedata_list(minLevel:maxLevel)
    !> Levelwise list of sources
    type(atl_source_type), intent(inout) :: source
    !> Levelwise list of penalization data
    type(atl_penalizationData_type), intent(inout) :: &
      & penalizationdata_list(minLevel:maxLevel)
    !> List of schemes, for each level.
    type(atl_scheme_type), intent(inout) :: scheme_list(minLevel:maxLevel)
    !> List of levelwise position of  projection method in unique
    !! projection_list
    integer, intent(in) :: poly_proj_pos(minLevel:maxLevel)
    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    !> The equation you are operating with.
    type(atl_equations_type),intent(in) :: equation
    !> Material parameter description.
    type(atl_material_type), intent(in) :: material_list(minlevel:maxlevel)
    !> General treelm settings
    type(tem_general_type), intent(inout) :: general
    !> Flag to indicate whether penalization terms should be computed or not.
    !!
    !! This is used to switch off the application of penalizing terms from
    !! the convective computation and instead compute it in an implicit
    !! update, see the atl_imexrk_module.
    !! Default is .true.!
    logical, intent(in), optional :: computePenalization
    ! -------------------------------------------------------------------- !
    logical :: usePenalization
    integer :: iDir, iFace
    ! -------------------------------------------------------------------- !

    usePenalization = .true.
    if (present(computePenalization)) usePenalization = computePenalization

    ! Interpolate fluxes from finer level to current level. Please keep in mind
    ! that the timestepping routine is called recursively for the levels.
    ! Therefore the compute rhs subroutine (i.e. the current routine) is
    ! executed on the finer levels before it is executed on the coarser level.
    ! So, the flux on the finer level is already available when a coarser level
    ! enters this routine.
    if( currentLevel .lt. maxLevel ) then
      call atl_modg_1d_fineToCoarseFace(          &
        & minLevel     = minLevel,                &
        & maxLevel     = maxLevel,                &
        & currentLevel = currentLevel,            &
        & mesh         = mesh_list,               &
        & facedata     = facedata_list,           &
        & nScalars     = equation%varSys%nScalars )
    end if

    ! Compute all the numerical fluxes, by the modal representations of the
    ! states on the faces. The modal representation on the faces is provided by
    ! the preprocess step.
    call tem_startTimer( timerHandle = atl_timerHandles%numFlux )
    select case(equation%eq_kind)
    case('euler_1d')
     call atl_modg_1d_euler_numFlux(             &
       & equation = equation,                    &
       & material = material_list(currentLevel), &
       & facedata = facedata_list(currentLevel)  )
    case('advection_1d')
     call atl_modg_1d_advection_numFlux(        &
       & mesh     = mesh_list(currentLevel),    &
       & equation = equation,                   &
       & facedata = facedata_list(currentLevel) )
    case('heat_1d')
     call atl_modg_1d_heat_numFlux(                    &
       & mesh      = mesh_list(currentLevel),          &
       & equation  = equation,                         &
       & facedata  = facedata_list(currentLevel),      &
       & scheme    = scheme_list(currentLevel)%modg_1d )
    case('loclineuler_1d')
     call atl_modg_1d_euler_numFlux(             &
       & equation = equation,                    &
       & material = material_list(currentLevel), &
       & facedata = facedata_list(currentLevel)  )
    case default
      call tem_abort( 'ERROR in compute_rhs_cubes: modg 1D is not supporting' &
        & // ' this PDE (num flux calculation), stopping...'                  )
    end select
    call tem_stopTimer( timerHandle = atl_timerHandles%numFlux )

    ! Exchange all the fluxes on the faces (for 1 direction, x)
    ! We send data about the fluxes at those faces where we received
    ! information about the states (in preprocess step). Therefore,
    ! we use the send buffer for receiving and the receive buffer for
    ! sending.
    call tem_startTimer( timerHandle = atl_timerHandles%commState )
    do iDir = 1,1
      do iFace = 1,2
        call general%commPattern%exchange_real(                            &
          & send         = mesh_list(currentLevel)                         &
          &                  %faces                                        &
          &                  %faces(iDir)                                  &
          &                  %recvBuffer_flux(iFace),                      &
          & recv         = mesh_list(currentLevel)                         &
          &                  %faces                                        &
          &                  %faces(iDir)                                  &
          &                  %sendBuffer_flux(iFace),                      &
          & state        = facedata_list(currentLevel)%faceFlux(iDir)%dat, &
          & message_flag = currentLevel,                                   &
          & comm         = general%proc%comm                               )
      end do
    end do
    call tem_stopTimer( timerHandle = atl_timerHandles%commState )

    call modg_1d_compute_project_physFlux(                              &
      & mesh             = mesh_list(currentLevel),                     &
      & equation         = equation,                                    &
      & material         = material_list(currentLevel),                 &
      & kerneldata       = kerneldata_list(currentLevel),               &
      & statedata        = statedata_list(currentLevel),                &
      & scheme           = scheme_list,                                 &
      & minLevel         = minLevel,                                    &
      & maxLevel         = maxLevel,                                    &
      & currentLevel     = CurrentLevel,                                &
      & poly_proj        = poly_proj_list(poly_proj_pos(currentLevel)), &
      & penalizationdata = penalizationdata_list(currentlevel),         &
      & usePenalization  = usePenalization                              )


    ! Fourth, project to test functions. This projects the numerical flux
    ! and the source terms to the test functions.
    call atl_modg_1d_project_testFunc(                          &
      & mesh             = mesh_list(currentLevel),             &
      & equation         = equation,                            &
      & kerneldata       = kerneldata_list(currentLevel),       &
      & penalizationdata = penalizationdata_list(currentLevel), &
      & usePenalization  = usePenalization,                     &
      & facedata         = facedata_list(currentLevel),         &
      & sourcedata       = source,                              &
      & level            = currentLevel,                        &
      & scheme           = scheme_list(currentLevel)%modg_1d    )

  end subroutine compute_rhs_cubes_modg_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine modg_1d_compute_project_physFlux( mesh, equation, material,   &
    &          kerneldata, statedata, scheme, poly_proj, penalizationdata, &
    &          minLevel, maxLevel, currentLevel, usePenalization           )
    use atl_cube_elem_module,                 only: atl_cube_elem_type
    use atl_scheme_module,                    only: atl_scheme_type
    use atl_equation_module,                  only: atl_equations_type
    use atl_penalization_module,              only: atl_penalizationData_type
    use atl_modg_1d_heat_kernel_module,       only: atl_modg_1d_heat_physFlux
    use atl_modg_1d_kernel_module,            only: &
      & atl_modg_1d_project_physFlux_testFunc
    use atl_modg_1d_multilevel_module,        only: atl_modg_1d_fineToCoarseFace
    use atl_modg_1d_euler_kernel_module,      only: &
      &   atl_modg_1d_euler_physFlux_const,         &
      &   atl_modg_1d_euler_physFlux_nonConst,      &
      &   atl_modg_1d_euler_penalization_const,     &
      &   atl_modg_1d_euler_penalization_nonConst
    use atl_physFluxEuler_1d_module,          only: atl_physFluxEuler_1d
    use atl_modg_1d_LoclinEuler_kernel_module,       &
      & only: atl_modg_1d_LoclinEuler_physFlux
    use atl_modg_1d_advection_kernel_module,  only: &
      & atl_modg_1d_advection_physFlux
    use atl_elemental_time_integration_module, only: &
      & atl_timestep_type
    use ply_oversample_module,                only: ply_convert2oversample,   &
      &                                             ply_convertFromoversample
    use atl_materialPrp_module,               only: atl_material_type
    use atl_facedata_module,                  only: atl_facedata_type
    use atl_kerneldata_module,                only: atl_kerneldata_type, &
      &                                             atl_statedata_type
    use ply_poly_project_module,              only: ply_poly_project_type, &
      &                                             assignment(=)
    use ply_dynArray_project_module,          only: dyn_projectionArray_type
    use ply_poly_project_module,              only: ply_poly_project_type, &
      &                                             ply_poly_project_n2m,  &
      &                                             assignment(=),         &
      &                                             ply_poly_project_m2n
    ! -------------------------------------------------------------------- !
    !> the levels of the geometry
    integer, intent(in) :: minLevel, maxLevel, currentLevel
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> Description of the material parameters.
    type(atl_material_type), intent(in) :: material
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The representation on the face + representation of the flux.
    type(atl_statedata_type), intent(in) :: statedata
    !> The parameters of the MODG scheme
    type(atl_scheme_type), intent(inout) :: scheme(minLevel:maxLevel)
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
    !> Penalization terms
    type(atl_penalizationData_type), intent(inout) :: penalizationdata
    !> Flag indicating whether to apply the penalization or not.
    !!
    !! When a implicit scheme is used to integrate the penalized parts, this
    !! can be used to switch it off here.
    logical, intent(in) :: usePenalization
    ! -------------------------------------------------------------------- !
    !> The direction
    integer :: iDir
    !> The modal coefficients of the current element in the loop.
    real(kind=rk), allocatable :: modalCoeffs(:,:)
    !> Nodal representation of the polynomial with in each cell.
    real(kind=rk), allocatable :: pointVal(:,:)
    !> The nodal representation of the physical flux along the 3 spatial
    !! directions.
    real(kind=rk), allocatable :: Res(:,:)
    real(kind=rk), allocatable :: tmp_state_der(:,:)
    integer :: nquadpoints, oversamp_dofs, iElem, ndofs, elemPos
    integer :: iVar
    integer :: elems
    logical :: use_linear_flux
    ! -------------------------------------------------------------------- !
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs
    nquadpoints = poly_proj%body_1D%nquadpoints
    nDofs = poly_proj%body_1D%ndofs


    allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars) )
    allocate( pointVal(nQuadPoints,equation%varSys%nScalars) )
    allocate( Res(nQuadPoints, equation%varSys%nScalars) )
    allocate( tmp_state_der(nDofs, equation%varSys%nScalars) )

    modalCoeffs = 0.0_rk

!ICE Intel 15
!ICE!    kerneldata%state_der = 0.0_rk
    ! Intel 15 workaround:
    call atl_initialize_state_der(kerneldata%state_der)

    ! Only one direction in 1D
    iDir = 1

    ! Compute the physical fluxes, by the modal representation in each cell.
    ! This operations is cell local.
    select case(trim(equation%eq_kind))
    case('euler_1d')
      use_linear_flux = .false.
      do elems = 1, material%material_desc%computeElems(atl_ConstMatIdx)%nElems
        iElem = material%material_desc%computeElems(atl_ConstMatIdx) &
          &                           %totElemIndices(elems)

        elemPos = mesh%descriptor%PntTID(iElem)
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )

        use_linear_flux = equation%euler%linear(                      &
          &                 mean      = statedata%state(iElem,1,:),   &
          &                 deviation = kerneldata%deviation(iElem,:) )

        if (use_linear_flux) then

          call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
          call atl_modg_1d_LoclinEuler_physFlux(       &
            & equation  = equation,                    &
            & state     = statedata%state(iElem,:,:),  &
            & poly_proj = poly_proj,                   &
            & res       = tmp_state_der                )
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )

        else
          if ( material%material_dat%mode_reducable(iElem) ) then
            call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
            tmp_state_der(:,:) = 0.0_rk
            tmp_state_der(1,:) = atl_physfluxEuler_1d(                &
              &  state           = statedata%state(iElem,1,:),        &
              &  isenCoeff       = equation%euler%isen_coef,          &
              &  penalty_scaling = material                           &
              &                      %material_dat                    &
              &                      %elemMaterialData(atl_ConstMatIdx) &
              &                      %materialDat(elems,1,1),         &
              &  U_o             = material                           &
              &                      %material_dat                    &
              &                      %elemMaterialData(atl_ConstMatIdx) &
              &                      %materialDat(elems,1,2)          )
            call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
            if ( penalizationData%isActive .and. usePenalization ) then
              do iVar=1,equation%varSys%nScalars
                pointval(:,iVar) = statedata%state(iElem,1,iVar)
              end do
            end if

          else
            ! get the modal coefficients of the current cell (for all variables
            ! of the Euler equation, therefore we use ":" for the third index).
            call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv)
            if (equation%euler%ensure_positivity) then
              call ply_convert2oversample(                      &
                & state       = statedata%state(iElem,:,:),     &
                & poly_proj   = poly_proj,                      &
                & nDim        = 1,                              &
                & modalCoeffs = modalCoeffs,                    &
                & ensure_positivity = [.true., .false., .true.] )
            else
              call ply_convert2oversample(                  &
                & state       = statedata%state(iElem,:,:), &
                & poly_proj   = poly_proj,                  &
                & nDim        = 1,                          &
                & modalCoeffs = modalCoeffs                 )
            end if

            ! Now, we transform the modal representation of this element to nodal
            ! space
            call ply_poly_project_m2n( me         = poly_proj,                &
             &                         dim        = 1,                        &
             &                         nVars      = equation%varSys%nScalars, &
             &                         nodal_data = pointVal,                 &
             &                         modal_data = modalCoeffs               )

            call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv)
            ! Pass the nodal values obtained above for evaluation of phys flux
            ! in pointwise fasion.
            ! Note the result is also in pointwise space so needs to be
            ! converted to modal values in next steps
            call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
            call atl_modg_1d_euler_physFlux_const(                          &
              & equation      = equation,                                   &
              & res           = Res,                                        &
              & PointVal      = pointval,                                   &
              & penalty_char  = material%material_dat                       &
              &                         %elemMaterialData(atl_ConstMatIdx)  &
              &                         %materialDat(elems,1,1),            &
              & U_o           =  material%material_dat                      &
              &                          %elemMaterialData(atl_ConstMatIdx) &
              &                          %materialDat(elems,1,2),           &
              & poly_proj     = poly_proj                                   )
          end if
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )

          ! If penalization is active then penalize the momentum
          ! and density also passing the pointwise data
          ! The penalized terms calculated are converted to modal values
          ! and stored in penalizationdata
          if (penalizationData%isActive .and. usePenalization) then
            call tem_startTimer( timerHandle = atl_timerHandles%pF_pen )
            call atl_modg_1d_euler_penalization_Const( &
              & equation         = equation,           &
              & poly_proj        = poly_proj,          &
              & nodal_data       = pointVal,           &
              & scheme_min       = scheme(minLevel),   &
              & penalizationData = penalizationData,   &
              & iElem            = elems,              &
              & material         = material            )
          end if
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_pen )

          if ( .not. material%material_dat%mode_reducable(iElem) ) then
            ! Convert the phys fluxes to modal data before projecting
            ! it onto test function
            ! --> oversamp modal space
            call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv )
            call ply_poly_project_n2m(                 &
              & me         = poly_proj,                &
              & dim        = 1,                        &
              & nVars      = equation%varSys%nScalars, &
              & nodal_data = Res,                      &
              & modal_data = modalCoeffs               )

            ! --> oversamp modal space
            call ply_convertFromOversample( modalCoeffs = modalCoeffs,  &
              &                             poly_proj   = poly_proj,    &
              &                             nDim        = 1,            &
              &                             state       = tmp_state_der )
          end if
            call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv )

        end if

        ! Now project the physical fluxes stored in temp_state_der
        call tem_startTimer(timerHandle = atl_timerHandles%pF_projectTestFunc)
        call atl_modg_1d_project_PhysFlux_testFunc(    &
          & equation   = equation,                     &
          & kerneldata = kerneldata,                   &
          & scheme     = scheme(currentLevel)%modg_1d, &
          & iElem      = iElem,                        &
          & nDofs      = nDofs,                        &
          & state_data = tmp_state_der                 )
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_projectTestFunc)

        call tem_stopTimer( me = atl_elemTimers , timerHandle = elemPos )
      end do ! elem loop
      do elems = 1, material%material_desc%computeElems(atl_VarMatIdx)%nElems
        iElem = material%material_desc%computeElems(atl_VarMatIdx) &
          &                           %totElemIndices(elems)

        elemPos = mesh%descriptor%totalPnt(iElem)
        call tem_startTimer( me = atl_elemTimers , TimerHandle = elemPos )

        ! Pass the nodal values obtained above for evaluation of phys flux
        ! in pointwise fasion.
        ! Note the result is also in pointwise space so needs to be
        ! converted to modal values in next steps
        if ( material%material_dat%mode_reducable(iElem) ) then
          call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
          tmp_state_der(:,:) = 0.0_rk
          tmp_state_der(1,:) = atl_physfluxEuler_1d(                &
            &  state           = statedata%state(iElem,1,:),        &
            &  isenCoeff       = equation%euler%isen_coef,          &
            &  penalty_scaling = material                           &
            &                      %material_dat                    &
            &                      %elemMaterialData(atl_VarMatIdx) &
            &                      %materialDat(elems,1,1),         &
            &  U_o             = material                           &
            &                      %material_dat                    &
            &                      %elemMaterialData(atl_VarMatIdx) &
            &                      %materialDat(elems,1,2)          )
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )
          if ( penalizationData%isActive .and. usePenalization ) then
            do iVar=1,equation%varSys%nScalars
              pointval(:,iVar) = statedata%state(iElem,1,iVar)
            end do
          end if
        else
          ! get the modal coefficients of the current cell (for all variables
          ! of the Euler equation, therefore we use ":" for the third index).
          call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv)
          if (equation%euler%ensure_positivity) then
            call ply_convert2oversample(                      &
              & state       = statedata%state(iElem,:,:),     &
              & poly_proj   = poly_proj,                      &
              & nDim        = 1,                              &
              & modalCoeffs = modalCoeffs,                    &
              & ensure_positivity = [.true., .false., .true.] )
          else
            call ply_convert2oversample(                  &
              & state       = statedata%state(iElem,:,:), &
              & poly_proj   = poly_proj,                  &
              & nDim        = 1,                          &
              & modalCoeffs = modalCoeffs                 )
          end if

          ! Now, we transform the modal representation of this element to
          ! nodal space
          call ply_poly_project_m2n( me         = poly_proj,                &
            &                        dim        = 1,                        &
            &                        nVars      = equation%varSys%nScalars, &
            &                        nodal_data = pointVal,                 &
            &                        modal_data = modalCoeffs               )
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv )

          call tem_startTimer( timerHandle = atl_timerHandles%pF_eval)
          call atl_modg_1d_euler_physFlux_nonConst(                       &
            &    equation      = equation,                                &
            &    res           = Res,                                     &
            &    PointVal      = pointval,                                &
            &    penalty_char  = material%material_dat                    &
            &                            %elemMaterialData(atl_VarMatIdx) &
            &                            %materialDat(elems,:,1),         &
            &    U_o           =  material%material_dat                   &
            &                            %elemMaterialData(atl_VarMatIdx) &
            &                            %materialDat(elems,:,2),         &
            &    poly_proj     = poly_proj                                )
        end if
          call tem_stopTimer( timerHandle = atl_timerHandles%pF_eval )

        ! If penalization is active then penalize the momentum
        ! and density also passing the pointwise data
        ! The penalized terms calculated are converted to modal values
        ! and stored in penalizationdata
        if ( penalizationData%isActive .and. usePenalization ) then
        call tem_startTimer( timerHandle = atl_timerHandles%pF_pen )
          call atl_modg_1d_euler_penalization_nonConst( &
            &    equation         = equation,           &
            &    poly_proj        = poly_proj,          &
            &    nodal_data       = pointVal,           &
            &    scheme_min       = scheme(minLevel),   &
            &    penalizationData = penalizationData,   &
            &    iElem            = elems,              &
            &    material         = material            )
        end if
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_pen )

        call tem_startTimer( timerHandle = atl_timerHandles%pF_projConv)
        if ( .not. material%material_dat%mode_reducable(iElem) ) then
          ! Convert the phys fluxes to modal data before projecting
          ! it onto test function
          ! --> oversamp modal space
          call ply_poly_project_n2m(                    &
            &    me         = poly_proj,                &
            &    dim        = 1,                        &
            &    nVars      = equation%varSys%nScalars, &
            &    nodal_data = Res,                      &
            &    modal_data = modalCoeffs               )

          ! --> oversamp modal space
          call ply_convertFromOversample( modalCoeffs = modalCoeffs,  &
            &                             poly_proj   = poly_proj,    &
            &                             nDim        = 1,            &
            &                             state       = tmp_state_der )
        end if
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_projConv )

        ! Now project the physical fluxes stored in temp_state_der
        call tem_startTimer(timerHandle = atl_timerHandles%pF_projectTestFunc)
        call atl_modg_1d_project_PhysFlux_testFunc(       &
          &    equation   = equation,                     &
          &    kerneldata = kerneldata,                   &
          &    scheme     = scheme(currentLevel)%modg_1d, &
          &    iElem      = iElem,                        &
          &    nDofs      = nDofs,                        &
          &    state_data = tmp_state_der                 )
        call tem_stopTimer( timerHandle = atl_timerHandles%pF_projectTestFunc)
        call tem_stopTimer( me = atl_elemTimers , timerHandle = elemPos )
      end do ! elem loop
    case('heat_1d')
      !Loop over Elements
      do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

        elemPos = mesh%descriptor%PntTID(iElem)
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )

        call atl_modg_1d_heat_physFlux(                           &
          & mesh            = mesh,                               &
          & equation        = equation,                           &
          & state           = statedata%state(iElem,:,:),         &
          & res             = tmp_state_der,                      &
          & modalCoeffs     = scheme(minLevel)%temp_modal(:,:,1), &
          & modalCoeffsDiff = scheme(minLevel)%temp_modal(:,:,2), &
          & poly_proj       = poly_proj                           )
        call atl_modg_1d_project_PhysFlux_testFunc(    &
          & equation   = equation,                     &
          & kerneldata = kerneldata,                   &
          & scheme     = scheme(currentLevel)%modg_1d, &
          & iElem      = iElem,                        &
          & nDofs      = nDofs,                        &
          & state_data = tmp_state_der                 )
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )
      end do ! elem loop

    case('advection_1d')
      !Loop over Elements
      do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

        elemPos = mesh%descriptor%PntTID(iElem)
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )

        call atl_modg_1d_advection_physFlux(       &
          & equation = equation,                   &
          & state    = statedata%state(iElem,:,:), &
          & res      = tmp_state_der               )
        call atl_modg_1d_project_PhysFlux_testFunc(    &
          & equation   = equation,                     &
          & kerneldata = kerneldata,                   &
          & scheme     = scheme(currentLevel)%modg_1d, &
          & iElem      = iElem,                        &
          & nDofs      = nDofs,                        &
          & state_data = tmp_state_der                 )

        call tem_stopTimer( me = atl_elemTimers , timerHandle = elemPos )
      end do ! elem loop

    case('loclineuler_1d')
      !Loop over Elements
      do iElem = 1, mesh%descriptor%elem%nElems(eT_fluid)

        elemPos = mesh%descriptor%PntTID(iElem)
        call tem_startTimer( me = atl_elemTimers , timerHandle = elemPos )

        call atl_modg_1d_LoclinEuler_physFlux(       &
          & equation  = equation,                    &
          & state     = statedata%state(iElem,:,:),  &
          & poly_proj = poly_proj,                   &
          & res       = tmp_state_der                )
        call atl_modg_1d_project_PhysFlux_testFunc(      &
          & equation     = equation,                     &
          & kerneldata   = kerneldata,                   &
          & scheme       = scheme(currentLevel)%modg_1d, &
          & iElem        = iElem,                        &
          & nDofs        = nDofs,                        &
          & state_data   = tmp_state_der                 )

        call tem_stopTimer( me = atl_elemTimers , timerHandle = elemPos )
      end do ! elem loop

    case default
      call tem_abort( 'ERROR in compute_rhs_cubes: modg is not supporting' &
        & // ' this PDE (phys flux calculation), stopping...'              )
    end select

  end subroutine modg_1d_compute_project_physFlux
  ! ************************************************************************ !


end module atl_compute_module
