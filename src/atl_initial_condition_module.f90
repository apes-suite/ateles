! Copyright (c) 2011-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2016,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2011 Gaurang Phadke <g.phadke@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

!> Initial conditions describe all field variables at the beginning of the
!! simulation.
!!
!! They take the form of a spatial function.
!! See [[tem_spatial_module]] for the definition of spatial functions.
!! In their simplest form they can just be a constant value.
!!
!! Usually, the variables of the field in the equation system need to be
!! defined, but for nonlinear flow equations the primitive variables are
!! required instead.
!! Variables with multiple components need to be provided as individual
!! spatial functions, and the names have an X, Y or Z appended to signify
!! the respective directions.
!!
!! An example initial condition with only constant values would look for
!! flows as follows:
!!
!!```lua
!!  initial_condition = {
!!    density = 1.0,
!!    velocityX = 0.0,
!!    velocityY = 0.0,
!!    velocityZ = 0.0,
!!    pressure = 1.0
!!  }
!!```
!!
!! The initial condition will be evaluated pointwise and then transformed
!! to a modal representation in Legendre polynomials.
module atl_initial_condition_module
  use env_module,                   only: rk, labelLen
  use flu_binding,                  only: flu_state
  use aot_path_module,              only: aot_path_type, &
    &                                     assignment(=), &
    &                                     aot_init_path

  use tem_aux_module,               only: tem_abort
  use tem_element_module,           only: eT_fluid
  use tem_ini_condition_module,     only: tem_load_ic, tem_ini_condition_type
  use tem_spatial_module,           only: tem_spatial_for
  use tem_logging_module,           only: logUnit
  use tem_element_module,           only: eT_fluid
  use treelmesh_module,             only: treelmesh_type

  use atl_equation_module,          only: atl_equations_type
  use atl_cube_container_module,    only: atl_cube_container_type
  use atl_scheme_module,            only: atl_modg_scheme_prp,    &
    &                                     atl_modg_2d_scheme_prp, &
    &                                     atl_modg_1d_scheme_prp
  use atl_reference_element_module, only: atl_refToPhysCoord
  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=),         &
    &                                     ply_poly_project_n2m
  use ply_oversample_module,        only: ply_convertFromoversample

  implicit none

  private

  public :: atl_load_initial_condition


contains


  ! ------------------------------------------------------------------------ !
  !> subroutine to load the initial conditions from a lua configuration file.
  subroutine atl_load_initial_condition( cube_container, equation, prj_pos,  &
    &                                    poly_proj_list, conf, luaFile, tree )
    ! -------------------------------------------------------------------- !
    !> Tree representation of your mesh.
    type(treelmesh_type), intent(in) :: tree

    !> Description of the equations to solve.
    type(atl_Equations_type), intent(inout) :: equation

    !> Handle, providing access to the configuration script
    type(flu_State),intent(inout) :: conf

    !> Name of Lua configuration file.
    character(len=*), intent(in) :: LuaFile

    !> Container which holds all cubic elements for the mesh.
    !!
    !! This parameter has to be initialzed already, since we need information
    !! about the physical coordinates of the cells.
    type(atl_cube_container_type), intent(inout) :: cube_container

    !> Levelwise position pointer for the projection method used for IC
    integer, intent(in) :: prj_pos( tree%global%minLevel  &
      &                             :tree%global%maxLevel )

    !> unique list for projection methods
    type(ply_poly_project_type), intent(inout) :: poly_proj_list(:)
    ! -------------------------------------------------------------------- !
    type(tem_ini_condition_type)          :: ic
    type(aot_path_type)                   :: path
    integer                               :: nVars
    integer                               :: iLevel, iVar, n_fluids
    integer                               :: nComponents, iState
    real(kind=rk)                         :: bcoord(3)
    real(kind=rk)                         :: dummycoord(1,3)
    character(len=labelLen), allocatable  :: stateName(:)
    !for oversampling
    integer                               :: nScalars
    integer                               :: iElem
    integer                               :: nquadpoints
    integer                               :: oversamp_dofs
    ! The current quadrature point on the physical element.
    real(kind=rk), allocatable :: physQuadLoc(:,:)
    ! The coordiantes of reference quadrature points
    real(kind=rk), allocatable :: refpoints(:,:)
    ! Array for the modal data
    real(kind=rk), allocatable :: modalData(:,:)
    ! Array for the nodal data
    real(kind=rk), allocatable :: nodalData(:,:,:)
    integer, allocatable :: errCode(:)
    ! -------------------------------------------------------------------- !
    !> @todo Do not set input here, but handle generic case
    !! (use variable system)
    ! read lua config file
    call aot_init_path(me = path, Filename = LuaFile)

    nScalars = equation%varSys%nScalars
    ! get variable names with X,Y,Z at the end for vector variables
    ! if the equation-system has primitive variables, use them
    allocate(StateName(nScalars))
    allocate(errCode(nSCalars))

    if (equation%hasPrimitiveVariables) then
      nVars = size(equation%primVar)
    else
      nVars = equation%varSys%nStateVars
    end if

    iState = 0
    do iVar = 1, nVars
      if( equation%hasPrimitiveVariables ) then
        nComponents = equation%varSys%method%val( &
          & equation%primVar(iVar) )%nComponents
        stateName(iState+1:iState+nComponents) = &
          & equation%varSys%Varname%val( &
          & equation%primVar(iVar) )
      else
        nComponents = equation%varSys%method%val( iVar )%nComponents
        stateName(iState+1:iState+nComponents) = &
          & equation%varSys%varName%val( iVar )
      end if

      select case(nComponents)
      case(1)
        iState = iState + 1
      case(2)
        iState = iState + 1
        stateName(iState) = trim(stateName(iState)) // 'X'
        iState = iState + 1
        stateName(iState) = trim(stateName(iState)) // 'Y'
      case(3)
        iState = iState + 1
        stateName(iState) = trim(stateName(iState)) // 'X'
        iState = iState + 1
        stateName(iState) = trim(stateName(iState)) // 'Y'
        iState = iState + 1
        stateName(iState) = trim(stateName(iState)) // 'Z'
      case default
        write(logUnit(1),*) 'unknown number of variable components &
          &                 (must be 1 or 3)'
        call tem_abort()
      end select
    end do

    ! Load initial condition specification
    call tem_load_ic( me        = ic,        &
      &               conf      = conf,      &
      &               stateName = stateName, &
      &               errCode   = errCode    )
     ! check if all variables are loaded correctly or an error occured
     do iState = 1, nScalars
       if ( errCode(iState)  /= 0 ) then
         write(logUnit(1),*) 'With Initial condition variable ', &
           &                 trim(stateName(iState))
         write(logUnit(1),*) ' something went wrong, maybe it is not defined.'
         call tem_abort()
       end if
     end do
     deallocate(errCode)


    ! Now we read out the initial conditions and set the right state variables
    ! in our kernel states. So we loop over all levels.
    levelLoop: do iLevel = tree%global%minLevel, tree%global%maxLevel
      !>\remark ATTENTION: we set the initial conditions only in fluid cells!
      !! so we have no initial condition for ghost and halo cells here!
      n_fluids = cube_container%mesh_list(iLevel) &
        &                      %descriptor        &
        &                      %elem              &
        &                      %nElems(eT_fluid)

      anyfluids: if (n_fluids > 0) then
        ! preset initial data
        cube_container%statedata_list(iLevel)%state = 0.0_rk

        ! Now we iterate over the different levels of the cubic elements and
        ! write the initial values to a temporary array.

        ! Initialization is different for the different schemes. So,
        ! we check which type of scheme we have and to the right
        ! initialization.
        select case(cube_container%scheme_list(iLevel)%scheme)
        case(atl_modg_scheme_prp)
          ! get the correct amount of number of quadpoints and dofs due to
          ! the proejction
          ! refpoints are the quad volume points, where the inital condition
          ! need to be evaluated
          ! oversamo_dofs and oversamp_degree is used for the oversampling
          ! loops
          nquadpoints = poly_proj_list(prj_pos(iLevel))%body_3D%nquadpoints
          oversamp_dofs= poly_proj_list(prj_pos(iLevel))%body_3D%oversamp_dofs


          if (all(ic%ini_state(:)%kind == 'const')) then

            ! If all states are constant, there is no need for pointwise
            ! evaluations.
            ! (Just filling the first mode = integral mean is sufficient)
            allocate(physQuadloc(1,3))
            allocate(NodalData(1,1,nScalars))
            physQuadloc = 0.0_rk
            cube_container%statedata_list(iLevel)%state = 0.0_rk
            do iVar=1,nScalars
              nodalData(:,1,iVar) &
                &  = tem_spatial_for( me    = ic%ini_state(iVar), &
                &                     n     = 1,                  &
                &                     coord = physQuadloc         )
            end do
            if (equation%hasPrimitiveVariables) then
              call equation%prim2cons(instate = NodalData, nElems = 1)
            end if
            do iVar=1,nScalars
              cube_container%statedata_list(iLevel)%state(:,1,iVar) &
                &  = nodalData(1,1,iVar)
            end do
            deallocate(physQuadloc)
            deallocate(NodalData)

          else

            ! Not all states are constant. If we need to convert primitive to
            ! conservative variables we need to fill all point data anyway.
            ! Otherwise, we can still fill the constant modes directly and only
            ! do the transformation for varying quantities.
            ! allocate tmpArray
            allocate(physQuadloc(nquadpoints,3))
            if (equation%hasPrimitiveVariables) then
              allocate(NodalData(1,nquadpoints,nScalars))
              allocate(modaldata(oversamp_dofs, nScalars))
            else
              allocate(NodalData(1,nquadpoints,1))
              allocate(modaldata(oversamp_dofs,1))
            end if

            allocate(refpoints(nquadpoints,3))

            refpoints = poly_proj_list(prj_pos(iLevel))%body_3d%nodes

            ! If the equation has primitive vars, we use them to load the
            ! initial conditions. And we need to transform them to
            ! conservative variables.
            if (equation%hasPrimitiveVariables) then
              ! Iterate over all the elements and project the spatial data to
              ! all the basis functions.
              do iElem = 1, n_fluids

                bcoord = cube_container%mesh_list(iLevel)%bary_coord(iElem,:)

                call  atl_refToPhysCoord(                                 &
                  &               refpoints  = refpoints,                 &
                  &               nPoints    = nquadpoints,               &
                  &               baryCoord  = bCoord,                    &
                  &               elemLength = cube_container             &
                  &                            %mesh_list(iLevel)%length, &
                  &               physPoints = physQuadLoc                )

                ! Evaluate spatial data at these points.
                do iVar=1,nScalars
                  NodalData(1,:,iVar) = tem_spatial_for(             &
                    &                       me = ic%ini_state(iVar), &
                    &                        n = nquadpoints,        &
                    &                    coord = physQuadLoc         )
                end do
                ! Convert to conservative (inplace, overwriting the original
                ! nodal data)
                call equation%prim2cons(instate = NodalData, nElems = 1)
                ! ---> nodal values



                call ply_poly_project_n2m( &
                  &      me = poly_proj_list(prj_pos(iLevel)), &
                  &      dim = 3 ,                             &
                  &      nVars = nScalars,                     &
                  &      nodal_data = NodalData(1,:,:),        &
                  &      modal_data = modalData                )
                ! ---> modal oversampling space
                ! --> oversamp modal space
                call ply_convertFromOversample(                         &
                  & modalCoeffs = modalData,                            &
                  & poly_proj   = poly_proj_list(prj_pos(iLevel)),      &
                  & nDim        = 3,                                    &
                  & state       = cube_container%statedata_list(iLevel) &
                  &                             %state(iElem,:,:)       )
                ! --> modal space

              end do !elemLoop
            else
              ! Iterate over all the elements and project the spatial data to
              ! all the basis functions.
              do iElem = 1, n_fluids

                bcoord = cube_container%mesh_list(iLevel)%bary_coord(iElem,:)

                call  atl_refToPhysCoord(                                 &
                  &               refpoints  = refpoints,                 &
                  &               nPoints    = nquadpoints,               &
                  &               baryCoord  = bCoord,                    &
                  &               elemLength = cube_container             &
                  &                            %mesh_list(iLevel)%length, &
                  &               physPoints = physQuadLoc                )

                  ! No conversion from primitive to conservative, process each
                  ! Variable on its own.
                dummycoord(1,:) = bcoord
                do iVar=1,nScalars
                  if (ic%ini_state(iVar)%kind == 'const') then
                    cube_container%statedata_list(iLevel) &
                      &           %state(iElem,:,iVar)   = 0.0_rk
                    cube_container%statedata_list(iLevel) &
                      &           %state(iElem,1:1,iVar)  &
                      &  = tem_spatial_for( me    = ic%ini_state(iVar), &
                      &                     coord =  dummycoord,        &
                      &                     n     = 1                   )
                  else
                    NodalData(1,:,1) = tem_spatial_for(             &
                      &                    me = ic%ini_state(iVar), &
                      &                     n = nquadpoints,        &
                      &                 coord = physQuadLoc         )

                    call ply_poly_project_n2m( &
                      &      me = poly_proj_list(prj_pos(iLevel)), &
                      &      dim = 3 ,                             &
                      &      nVars = 1,                            &
                      &      nodal_data = NodalData(1,:,:),        &
                      &      modal_data = modalData                )
                    ! ---> modal oversampling space
                    ! --> oversamp modal space
                    call ply_convertFromOversample(                           &
                      & modalCoeffs = modalData,                              &
                      & poly_proj   = poly_proj_list(prj_pos(iLevel)),        &
                      & nDim        = 3,                                      &
                      & state       = cube_container%statedata_list(iLevel)   &
                      &                             %state(iElem,:,iVar:iVar) )
                    ! --> modal space

                  end if
                end do
              end do !elemLoop
            end if

            deallocate(refpoints)
            deallocate(physQuadloc)
            deallocate(NodalData)
            deallocate(ModalData)

          end if

        case(atl_modg_2d_scheme_prp)
          ! get the correct amount of number of quadpoints and dofs due to
          ! the proejction
          ! refpoints are the quad volume points, where the inital condition
          ! need to be evaluated
          ! oversamo_dofs and oversamp_degree is used for the oversampling
          ! loops
          nquadpoints = poly_proj_list(prj_pos(iLevel))%body_2D%nquadpoints
          oversamp_dofs= poly_proj_list(prj_pos(iLevel))%body_2D%oversamp_dofs

          ! allocate tmpArray
          allocate(physQuadloc(nquadpoints,3))
          allocate(NodalData(1,nquadpoints,nScalars))
          allocate(ModalData(oversamp_dofs,nScalars))

          allocate(refpoints(nquadpoints,3))
          refpoints = poly_proj_list(prj_pos(iLevel))%body_2d%nodes

          ! Iterate over all the elements and project the spatial data to
          ! all the basis functions.
          do iElem = 1, n_fluids

            ! Move the quadrature points to the physical element position
            bcoord = cube_container%mesh_list(iLevel)%bary_coord(iElem,:)

            call  atl_refToPhysCoord(                                 &
              &               refpoints  = refpoints,                 &
              &               nPoints    = nquadpoints,               &
              &               baryCoord  = bCoord,                    &
              &               elemLength = cube_container             &
              &                            %mesh_list(iLevel)%length, &
              &               physPoints = physQuadLoc                )


            ! Evaluate spatial data at these points.
            do iVar=1,nScalars
              NodalData(1,:,iVar) = tem_spatial_for(             &
                &                       me = ic%ini_state(iVar), &
                &                        n = nquadpoints,        &
                &                    coord = physQuadLoc         )
            end do
            ! If the equation has primitive vars, we use them to load the
            ! initial conditions. And we need to transform them to
            ! conservative variables.
            if (equation%hasPrimitiveVariables) then
              ! Convert to conservative (inplace, overwriting the original
              ! nodal data)
              call equation%prim2cons(instate = NodalData, nElems = 1)
            end if
            ! ---> nodal values


            call ply_poly_project_n2m( &
              &  me = poly_proj_list(prj_pos(iLevel)), &
              &  dim = 2 ,                             &
              &  nVars = nScalars,                     &
              &  nodal_data = NodalData(1,:,:),        &
              &  modal_data = modalData                )
            ! ---> modal oversampling space
            call ply_convertFromOversample(                         &
              & modalCoeffs = modalData,                            &
              & poly_proj   = poly_proj_list(prj_pos(iLevel)),      &
              & nDim        = 2,                                    &
              & state       = cube_container%statedata_list(iLevel) &
              &                             %state(iElem,:,:)       )
            ! ---> modal space


          end do !elemLoop

          deallocate(refpoints)
          deallocate(physQuadloc)
          deallocate(NodalData)
          deallocate(ModalData)

        case(atl_modg_1d_scheme_prp)
          ! get the correct amount of number of quadpoints and dofs due to
          ! the proejction
          ! refpoints are the quad volume points, where the inital condition
          ! need to be evaluated
          ! oversamo_dofs and oversamp_degree is used for the oversampling
          ! loops
          nquadpoints = poly_proj_list(prj_pos(iLevel)) &
            &                         %body_1D%nquadpoints
          oversamp_dofs= poly_proj_list(prj_pos(iLevel)) &
            &                          %body_1D%oversamp_dofs

          ! allocate tmpArray
          allocate(physQuadloc(nquadpoints,3))
          allocate(NodalData(1,nquadpoints,nScalars))
          allocate(ModalData(oversamp_dofs, nScalars))

          allocate(refpoints(nquadpoints,3))
          refpoints = poly_proj_list(prj_pos(iLevel))%body_1d%nodes

          ! Iterate over all the elements and project the spatial data to
          ! all the basis functions.
          do iElem = 1, n_fluids

            ! Move the quadrature points to the physical element position
            bcoord = cube_container%mesh_list(iLevel)%bary_coord(iElem,:)

            call  atl_refToPhysCoord(                                 &
              &               refpoints  = refpoints,                 &
              &               nPoints    = nquadpoints,               &
              &               baryCoord  = bCoord,                    &
              &               elemLength = cube_container             &
              &                            %mesh_list(iLevel)%length, &
              &               physPoints = physQuadLoc                )

            ! Evaluate spatial data at these points.
            do iVar=1,nScalars
              NodalData(1,:,iVar) = tem_spatial_for(             &
                &                       me = ic%ini_state(iVar), &
                &                        n = nquadpoints,        &
                &                    coord = physQuadLoc         )
            end do
            ! If the equation has primitive vars, we use them to load the
            ! initial conditions. And we need to transform them to
            ! conservative variables.
            if (equation%hasPrimitiveVariables) then
              ! Convert to conservative (inplace, overwriting the original
              ! nodal data)
              call equation%prim2cons(instate = NodalData, nElems = 1)
            end if
            ! ---> nodal values


            call ply_poly_project_n2m( &
              &  me = poly_proj_list(prj_pos(iLevel)), &
              &  dim = 1 ,                             &
              &  nVars = nScalars,                     &
              &  nodal_data = NodalData(1,:,:),        &
              &  modal_data = modalData                )
            ! --> oversamp modal space
            call ply_convertFromOversample(                          &
               & modalCoeffs = modalData,                            &
               & poly_proj   = poly_proj_list(prj_pos(iLevel)),      &
               & nDim        = 1,                                    &
               & state       = cube_container%statedata_list(iLevel) &
               &                             %state(iElem,:,:)       )


          end do !elemLoop

          deallocate(refpoints)
          deallocate(physQuadloc)
          deallocate(NodalData)
          deallocate(ModalData)

        case default
          write(logUnit(1),*) 'ERROR in load_initial_conidtion: not able to' &
            & // ' handle initial conditions for this type of' &
            & // ' schemes, stopping ...'
          call tem_abort()
        end select

      end if anyFluids

    end do levelLoop

    write(logUnit(1),*) 'finished reading'

  end subroutine atl_load_initial_condition
  ! ------------------------------------------------------------------------ !

end module atl_initial_condition_module
