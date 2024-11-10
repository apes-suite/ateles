! Copyright (c) 2014-2016, 2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
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

!> Configuring the projection (modal to nodal transformation)
!!
!! Ateles uses the Legendre polynomial basis to represent the solution in the
!! DG elements. For nonlinear operations this incurs the necessity to transform
!! this representation into nodal values, such that a pointwise evaluation can
!! be performed for those nonlinear operations.
!! The method to be used for this projection between the two representations is
!! configured in the `projection` table.
!!
!! There is a general projection method to be configured according to the
!! described in [[ply_poly_project_module]] and in addition individual
!! projection parameters can be defined for various tasks. These override the
!! general parameters for those tasks.
!!
!! The general parameters will be used for all tasks for which no specific
!! projection definitions are provided and are the ones that will be used in the
!! main computational operations.
!!
!! For the following tasks individual projection definitions can be provided
!! within the `projection` table (each needs to be a table describing the
!! projection method, as layed out in [[ply_poly_project_module]]):
!!
!! * `source_terms`: evaluation of source terms
!! * `boundary_condition`: evaluation of boundary condition
!! * `initial_condition`: evaluation of initial condition
!! * `material`: evaluation of material parameters
!!
!! The projection definition for ateles could for example look like this:
!!
!!```lua
!!  projection = {
!!    kind = 'fpt',
!!    factor = 1.5,
!!    initial_condition = {
!!      kind = 'fpt',
!!      factor = 3
!!    }
!!  }
!!```
!!
!! Please refer to the [[ply_poly_project_module]] for details on the
!! configuration of the projection method.
module atl_load_project_module
  use env_module,                  only: rk
  use flu_binding,                 only: flu_state
  use aot_table_module,            only: aot_table_open, &
    &                                    aot_table_close
  use tem_aux_module,              only: tem_abort
  use tem_logging_module,          only: logUnit
  use ply_dynArray_project_module, only: dyn_ProjectionArray_type, &
    &                                    ply_fill_dynProjectArray, &
    &                                    ply_prj_init_type, append

  use atl_materialPrp_module,      only: atl_material_type
  use atl_scheme_module,           only: atl_scheme_type,        &
    &                                    atl_modg_scheme_prp,    &
    &                                    atl_modg_2d_scheme_prp, &
    &                                    atl_modg_1d_scheme_prp

  use atl_boundary_module,         only: atl_level_boundary_type
  use atl_equation_module,         only: atl_equations_type

  implicit none

  private

  public :: atl_load_projection


contains


  ! ------------------------------------------------------------------------ !
  subroutine atl_load_projection( minLevel, maxLevel, equation, poly_proj_pos, &
    &                             dynprojectArray, conf, scheme_list,          &
    &                             source_projPos, boundary_list,               &
    &                             boundary_list_stab, ic_pos_list,             &
    &                             material_list                                )
    ! -------------------------------------------------------------------- !
    integer, intent(in)                 :: minLevel
    integer, intent(in)                 :: maxLevel
    type(atl_equations_type), intent(in) :: equation
    integer, intent(inout)              :: poly_proj_pos(minlevel:maxLevel)
    type(dyn_ProjectionArray_type), intent(inout) :: dynprojectArray
    type(flu_State), intent(in)         :: conf
    type(atl_scheme_type), intent(in)   :: scheme_list(minlevel:maxLevel)
    type(atl_level_boundary_type), intent(inout) &
      & :: boundary_list(minlevel:maxLevel), &
      &    boundary_list_stab(minlevel:maxLevel)
    !> The levelwise positions of the projections that has to be used for the
    !! sources
    integer, intent(inout)              :: source_projPos(minlevel:maxLevel)
    integer, intent(inout), allocatable :: ic_pos_list(:)
    type(atl_material_type), intent(inout) :: material_list(minlevel:maxLevel)
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    integer :: thandle
    integer :: mathandle
    integer :: matProjectPos
    integer :: genProjectPos
    real(kind=rk) :: matProjFact, State2MatFact
    type(ply_prj_init_type) :: proj_state2Mat
    logical :: enable_mat_project
    ! -------------------------------------------------------------------- !

    allocate (ic_pos_list(minLevel:maxLevel))
    write(logUnit(1),*) 'Loading projection tables...'
    ! open the projection table
    call aot_table_open(L = conf,  thandle = thandle,  key = 'projection')

    ! load general projection
    call atl_load_general_projection( poly_proj_pos   = poly_proj_pos,   &
      &                               dynprojectArray = dynprojectArray, &
      &                               minLevel        = minLevel,        &
      &                               maxLevel        = maxLevel,        &
      &                               scheme_list     = scheme_list,     &
      &                               conf            = conf,            &
      &                               parent          = thandle          )

    prjtable: if (thandle /= 0) then

      ! load individual projection method for source terms
      call atl_load_subprojection( minLevel         = minLevel,         &
        &                          maxLevel         = maxLevel,         &
        &                          subkey           = 'source_terms',   &
        &                          proj_pos         = source_projPos,   &
        &                          dynprojectArray  = dynprojectArray,  &
        &                          general_pos_list = poly_proj_pos,    &
        &                          scheme_list      = scheme_list,      &
        &                          conf             = conf,             &
        &                          parent           = thandle           )

      ! load individual projection method for initial condition
      call atl_load_subprojection( minLevel         = minLevel,            &
        &                          maxLevel         = maxLevel,            &
        &                          subkey           = 'initial_condition', &
        &                          proj_pos         = ic_pos_list,         &
        &                          dynprojectArray  = dynprojectArray,     &
        &                          general_pos_list = poly_proj_pos,       &
        &                          scheme_list      = scheme_list,         &
        &                          conf             = conf,                &
        &                          parent           = thandle              )

      ! load individual projection method for boundary conditions
      call atl_load_subprojection( minLevel         = minLevel,             &
        &                          maxLevel         = maxLevel,             &
        &                          subkey           = 'boundary_condition', &
        &                          proj_pos         = boundary_list(:)      &
        &                                               %poly_proj_pos,     &
        &                          dynprojectArray  = dynprojectArray,      &
        &                          general_pos_list = poly_proj_pos,        &
        &                          scheme_list      = scheme_list,          &
        &                          conf             = conf,                 &
        &                          parent           = thandle               )
      boundary_list_stab(:)%poly_proj_pos = boundary_list(:)%poly_proj_pos

      enable_mat_project = .true.
      select case(trim(equation%eq_kind))
      case( 'euler', 'euler_2d', 'euler_1d',                      &
        &   'navier_stokes', 'navier_stokes_2d',                  &
        &   'filtered_navier_stokes', 'filtered_navier_stokes_2d' )
        enable_mat_project = ( equation%euler%porosity     &
          &                    >= (1.0_rk-epsilon(1.0_rk)) )
        call aot_table_open( L       = conf,      &
          &                  parent  = thandle,   &
          &                  thandle = mathandle, &
          &                  key     = 'material' )
        if (.not. enable_mat_project .and. mathandle /= 0) then
          write(logUnit(1),*) 'WARNING: material projection provided but' &
            &                 // ' ignored!'
          write(logUnit(1),*) 'Not supported for porosity < 1.'
          write(logUnit(1),*) 'The general projection setting will be used'
          write(logUnit(1),*) 'for the material projection.'
        end if
        call aot_table_close(L = conf, thandle = mathandle)
      end select

      act_mat: if (enable_mat_project) then
        ! load individual projection method for material
        call atl_load_subprojection( minLevel         = minLevel,         &
          &                          maxLevel         = maxLevel,         &
          &                          subkey           = 'material',       &
          &                          proj_pos         = material_list(:)  &
          &                                               %poly_proj_pos, &
          &                          dynprojectArray  = dynprojectArray,  &
          &                          general_pos_list = poly_proj_pos,    &
          &                          scheme_list      = scheme_list,      &
          &                          conf             = conf,             &
          &                          parent           = thandle           )
        do iLevel = minLevel, maxLevel
          ! Copy the material transformation and change the factor and the
          ! polynomial degree to match with the parameter of the state
          ! transformation.
          write(logUnit(1),*) 'Preparing special state to material' &
            &                 // ' transformation'
          matProjectPos = material_list(iLevel)%poly_proj_pos
          genProjectPos = poly_proj_pos(iLevel)
          proj_state2Mat = dynProjectArray%val( matProjectPos )
          if (proj_state2Mat%header%kind == 'fpt') then
            matProjFact = dynProjectArray%val(matProjectpos)%header &
              &                          %fpt_header%factor
          else
            matProjFact = dynProjectArray%val(matProjectpos)%header &
              &                          %l2p_header%factor
          end if
          State2MatFact = matProjFact &
            &             * real(dynProjectArray%val(matProjectPos) &
            &                                   %maxPolyDegree+1,   &
            &                    kind=rk)                           &
            &             / real(dynProjectArray%val(genProjectPos) &
            &                                   %maxPolyDegree+1,   &
            &                    kind=rk)
          if (proj_state2Mat%header%kind == 'fpt') then
            proj_state2Mat%header%fpt_header%factor = State2MatFact
          else
            proj_state2Mat%header%l2p_header%factor = State2MatFact
          end if

          proj_state2Mat%maxPolyDegree = dynProjectArray%val(genProjectPos) &
            &                                           %maxPolyDegree

          call append( me  = dynProjectArray,                              &
            &          val = proj_state2Mat,                               &
            &          pos = material_list(iLevel)%poly_proj_pos_state2Mat )

          write(logUnit(1),*) ' * polydegree: ', proj_state2Mat%maxPolyDegree
          if (proj_state2Mat%header%kind == 'fpt') then
            write(logUnit(1),*) ' * FPT Oversamp factor: ', &
              &                 proj_state2Mat%header%fpt_header%factor
          else
            write(logUnit(1),*) ' * L2P Oversamp factor: ', &
              &                 proj_state2Mat%header%l2p_header%factor
          end if
        end do

      else act_mat

        ! Ignore material projection and fall back to general projection.
        proj_state2Mat = dynProjectArray%val( poly_proj_pos(iLevel) )
        if (proj_state2Mat%header%kind == 'fpt') then
          proj_state2Mat%header%fpt_header%factor = 1.0_rk
        else
          proj_state2Mat%header%l2p_header%factor = 1.0_rk
        end if
        call append( me  = dynProjectArray,                              &
          &          val = proj_state2Mat,                               &
          &          pos = material_list(iLevel)%poly_proj_pos_state2Mat )
      end if act_mat

    else prjtable

      ! only load the general one and copy all the other position pointers
      do iLevel = minLevel, maxLevel
        boundary_list(iLevel)%poly_proj_pos = poly_proj_pos(iLevel)
        boundary_list_stab(iLevel)%poly_proj_pos = poly_proj_pos(iLevel)
        source_projPos(iLevel) = poly_proj_pos(iLevel)
        ic_pos_list(iLevel) = poly_proj_pos(iLevel)
        material_list(iLevel)%poly_proj_pos = poly_proj_pos(iLevel)
        proj_state2Mat = dynProjectArray%val( poly_proj_pos(iLevel) )
        if (proj_state2Mat%header%kind == 'fpt') then
          proj_state2Mat%header%fpt_header%factor = 1.0_rk
        else
          proj_state2Mat%header%l2p_header%factor = 1.0_rk
        end if
        call append( me  = dynProjectArray,                              &
          &          val = proj_state2Mat,                               &
          &          pos = material_list(iLevel)%poly_proj_pos_state2Mat )
      end do

    end if prjtable

    call aot_table_close(L = conf, thandle = thandle)

  end subroutine atl_load_projection
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine atl_load_general_projection( poly_proj_pos, dynprojectArray,  &
    &                                     minLevel, maxLevel, scheme_list, &
    &                                     conf, parent                     )
    ! -------------------------------------------------------------------- !
    integer, intent(in)           :: minLevel
    integer, intent(in)           :: maxLevel
    integer, intent(inout) :: poly_proj_pos(minlevel:maxLevel)
    type(dyn_ProjectionArray_type), intent(inout) :: dynprojectArray
    type(atl_scheme_type), intent(in)   :: scheme_list(minlevel:maxLevel)
    type(flu_State), intent(in)   :: conf
    integer, intent(in)           :: parent
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    ! -------------------------------------------------------------------- !

    do iLevel= minLevel, maxLevel
      write(logUnit(1),*) ''
      write(logUnit(1),*) '...for Level', iLevel

      select case(scheme_list(iLevel)%scheme)
      case (atl_modg_scheme_prp)
        call ply_fill_dynProjectArray(                                    &
          & proj_pos            = poly_proj_pos(iLevel),                  &
          & dyn_projectionArray = dynProjectArray,                        &
          & maxPolyDegree       = scheme_list(iLevel)%modg%maxPolyDegree, &
          & basisType           = scheme_list(iLevel)%modg%basisType,     &
          & conf                = conf,                                   &
          & parent              = parent                                  )

      case (atl_modg_2d_scheme_prp)
        call ply_fill_dynProjectArray(                                       &
          & proj_pos            = poly_proj_pos(iLevel),                     &
          & dyn_projectionArray = dynProjectArray,                           &
          & maxPolyDegree       = scheme_list(iLevel)%modg_2d%maxPolyDegree, &
          & basisType           = scheme_list(iLevel)%modg_2d%basisType,     &
          & conf                = conf,                                      &
          & parent              = parent                                     )

      case (atl_modg_1d_scheme_prp)
        call ply_fill_dynProjectArray(                                       &
          & proj_pos            = poly_proj_pos(iLevel),                     &
          & dyn_projectionArray = dynProjectArray,                           &
          & maxPolyDegree       = scheme_list(iLevel)%modg_1d%maxPolyDegree, &
          & basisType           = scheme_list(iLevel)%modg_1d%basisType,     &
          & conf                = conf,                                      &
          & parent              = parent                                     )

      case default
        write(logUnit(1),*) 'ERROR in init_projection: unknown scheme name '
        write(logUnit(1),*) 'Stopping....'
        call tem_abort()

      end select

    end do
  end subroutine atl_load_general_projection
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine atl_load_subprojection( minLevel, maxLevel, subkey, proj_pos,    &
    &                                dynprojectArray, general_pos_list, conf, &
    &                                scheme_list, parent                      )
    ! -------------------------------------------------------------------- !
    integer, intent(in)           :: minLevel
    integer, intent(in)           :: maxLevel
    character(len=*), intent(in)  :: subkey
    integer, intent(out)          :: proj_pos(minlevel:maxlevel)
    type(dyn_ProjectionArray_type), intent(inout) :: dynprojectArray
    integer, intent(in)           :: general_pos_list(minLevel:maxLevel)
    type(atl_scheme_type), intent(in) :: scheme_list(minlevel:maxLevel)
    type(flu_State), intent(in)   :: conf
    integer, intent(in)           :: parent
    ! -------------------------------------------------------------------- !
    integer :: iLevel
    integer :: thandle
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) ''
    write(logUnit(1),*) 'Load the individual projection method for ' &
      &                 // 'polynomials for ', trim(subkey), '...'

    call aot_table_open( L       = conf,        &
      &                  parent  = parent,      &
      &                  thandle = thandle,     &
      &                  key     = trim(subkey) )

    if (thandle /= 0) then ! --> there is a projection table for the subkey
      ! levelwise we open the projection table and fill the projection_init_type
      ! to append the DA for projection.

      do iLevel= minLevel, maxLevel
        select case(scheme_list(iLevel)%scheme)
        case (atl_modg_scheme_prp)
          call ply_fill_dynProjectArray(                                   &
            &    proj_pos            = proj_pos(iLevel),                   &
            &    dyn_projectionArray = dynProjectArray,                    &
            &    maxPolyDegree       = scheme_list(iLevel)                 &
            &                            %modg                             &
            &                            %maxPolyDegree,                   &
            &    basisType           = scheme_list(iLevel)%modg%basisType, &
            &    conf                = conf,                               &
            &    parent              = thandle                             )

        case (atl_modg_2d_scheme_prp)
          call ply_fill_dynProjectArray(                   &
            &    proj_pos            = proj_pos(iLevel),   &
            &    dyn_projectionArray = dynProjectArray,    &
            &    maxPolyDegree       = scheme_list(iLevel) &
            &                            %modg_2d          &
            &                            %maxPolyDegree,   &
            &    basisType           = scheme_list(iLevel) &
            &                            %modg_2d          &
            &                            %basisType,       &
            &    conf                = conf,               &
            &    parent              = thandle             )

        case (atl_modg_1d_scheme_prp)
          call ply_fill_dynProjectArray(                   &
            &    proj_pos            = proj_pos(iLevel),   &
            &    dyn_projectionArray = dynProjectArray,    &
            &    maxPolyDegree       = scheme_list(iLevel) &
            &                            %modg_1d          &
            &                            %maxPolyDegree,   &
            &    basisType           = scheme_list(iLevel) &
            &                            %modg_1d          &
            &                            %basisType,       &
            &    conf                = conf,               &
            &    parent              = thandle             )

        end select
      end do

    else

      ! There is no projection table, means the position of the general
      ! projection for all level is used
      ! position should be stored in the source datatype
      write(logUnit(1),*) ' No individual projection table for ' &
        &                 //trim(subkey) &
        &                 //' is provided, using general projection.'
      ! store the position of the general projection method
      proj_pos = general_pos_list

    end if

    call aot_table_close(L=conf, thandle=thandle)

  end subroutine atl_load_subprojection
  ! ------------------------------------------------------------------------ !


end module atl_load_project_module
