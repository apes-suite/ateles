! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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

! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2015, 2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!--------------------------------------------
!    A O S - Array of structures layout new
!-------------------------------------------
! Access to get_point value output
! Access to get_element value output

!> Routines to derive quantities from the state in the Maxwell equation system.
module atl_eqn_maxwell_derive_module
  use, intrinsic :: iso_c_binding,  only: c_f_pointer
  use env_module,                   only: rk

  use tem_varSys_module,            only: tem_varSys_type,            &
    &                                     tem_varSys_op_type
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type

  use atl_equation_module,          only: atl_equations_type

  implicit none

  private

  public :: atl_xFrom3D_getElement, atl_xFrom3D_getPoint
  public :: atl_yFrom3D_getElement, atl_yFrom3D_getPoint
  public :: atl_zFrom3D_getElement, atl_zFrom3D_getPoint
  public :: atl_eqn_maxwell_cons2prim
  public :: atl_eqn_maxwell_prim2cons

contains

  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_maxwell_prim2cons(equation, instate, outstate, material )
    ! ------------------------------------------------------------------------ !
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in) :: material(:,:)
    ! ------------------------------------------------------------------------ !

    if( size(material,1).eq.1 ) then
      if(present(outstate)) then
        ! displacement field, i.e. D from
        outstate(:,1) = instate(:,1) * material(1,1)
        outstate(:,2) = instate(:,2) * material(1,1)
        outstate(:,3) = instate(:,3) * material(1,1)

        ! magnetic field, i.e. B from H
        outstate(:,4) = instate(:,4) * material(1,2)
        outstate(:,5) = instate(:,5) * material(1,2)
        outstate(:,6) = instate(:,6) * material(1,2)
      else
        ! displacement field, i.e. D from
        instate(:,1) = instate(:,1) * material(1,1)
        instate(:,2) = instate(:,2) * material(1,1)
        instate(:,3) = instate(:,3) * material(1,1)

        ! magnetic field, i.e. B from H
        instate(:,4) = instate(:,4) * material(1,2)
        instate(:,5) = instate(:,5) * material(1,2)
        instate(:,6) = instate(:,6) * material(1,2)
      end if
    else
      if(present(outstate)) then
        ! displacement field, i.e. D from
        outstate(:,1) = instate(:,1) * material(:,1)
        outstate(:,2) = instate(:,2) * material(:,1)
        outstate(:,3) = instate(:,3) * material(:,1)

        ! magnetic field, i.e. B from H
        outstate(:,4) = instate(:,4) * material(:,2)
        outstate(:,5) = instate(:,5) * material(:,2)
        outstate(:,6) = instate(:,6) * material(:,2)
      else
        ! displacement field, i.e. D from
        instate(:,1) = instate(:,1) * material(:,1)
        instate(:,2) = instate(:,2) * material(:,1)
        instate(:,3) = instate(:,3) * material(:,1)

        ! magnetic field, i.e. B from H
        instate(:,4) = instate(:,4) * material(:,2)
        instate(:,5) = instate(:,5) * material(:,2)
        instate(:,6) = instate(:,6) * material(:,2)
      end if
    end if

  end subroutine atl_eqn_maxwell_prim2cons


  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_maxwell_cons2prim(equation, instate, outstate, material )
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! --------------------------------------------------------------------------!

    if(size(material,1).eq.1) then
      if(present(outstate)) then
        ! electric field, i.e. E from D
        outstate(:,1) = instate(:,1) / material(1,1)
        outstate(:,2) = instate(:,2) / material(1,1)
        outstate(:,3) = instate(:,3) / material(1,1)

        ! magnetizing field, i.e. H from B
        outstate(:,4) = instate(:,4) / material(1,2)
        outstate(:,5) = instate(:,5) / material(1,2)
        outstate(:,6) = instate(:,6) / material(1,2)
      else
        ! electric field, i.e. E from D
        instate(:,1) = instate(:,1) / material(1,1)
        instate(:,2) = instate(:,2) / material(1,1)
        instate(:,3) = instate(:,3) / material(1,1)

        ! magnetizing field, i.e. H from B
        instate(:,4) = instate(:,4) / material(1,2)
        instate(:,5) = instate(:,5) / material(1,2)
        instate(:,6) = instate(:,6) / material(1,2)
      end if
    else
      if(present(outstate)) then
        ! electric field, i.e. E from D
        outstate(:,1) = instate(:,1) / material(:,1)
        outstate(:,2) = instate(:,2) / material(:,1)
        outstate(:,3) = instate(:,3) / material(:,1)

        ! magnetizing field, i.e. H from B
        outstate(:,4) = instate(:,4) / material(:,2)
        outstate(:,5) = instate(:,5) / material(:,2)
        outstate(:,6) = instate(:,6) / material(:,2)
      else
        ! electric field, i.e. E from D
        instate(:,1) = instate(:,1) / material(:,1)
        instate(:,2) = instate(:,2) / material(:,1)
        instate(:,3) = instate(:,3) / material(:,1)

        ! magnetizing field, i.e. H from B
        instate(:,4) = instate(:,4) / material(:,2)
        instate(:,5) = instate(:,5) / material(:,2)
        instate(:,6) = instate(:,6) / material(:,2)
      end if
    end if
  end subroutine atl_eqn_maxwell_cons2prim


  subroutine atl_xFrom3D_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nPnts*3)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = input                                    )

    do iPoint = 1, nPnts
      res((ipoint-1)*1+1) = input((ipoint-1)*3+1)
    end do

  end subroutine atl_xFrom3D_getPoint

  subroutine atl_xFrom3D_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                     nDofs, res                          )
    ! ---------------------------------------------------------------------------
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

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nDofs*nElems*3)
    integer :: iElem, iDof
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = input                                      )

    !! We only want to have one component, thus we don't need to loop over them
    do iElem = 1, nElems
      do iDof = 1, nDofs
        res(( ielem-1)* 1* ndofs+( idof-1)* 1+1) =    &
          & input(( ielem-1)* 3* ndofs+( idof-1)* 3+1)
      end do
    end do

  end subroutine atl_xFrom3D_getElement

  subroutine atl_yFrom3D_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nPnts*3)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = input                                    )

    do iPoint = 1, nPnts
      res((ipoint-1)*1+1) = input((ipoint-1)*3+2)
    end do

  end subroutine atl_yFrom3D_getPoint

  subroutine atl_yFrom3D_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                     nDofs, res                          )
    ! ---------------------------------------------------------------------------
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

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nDofs*nElems*3)
    integer :: iElem, iDof
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = input                                      )

    !! We only want to have one component, thus we don't need to loop over them
    do iElem = 1, nElems
      do iDof = 1, nDofs
        res(( ielem-1)* 1* ndofs+( idof-1)* 1+1) =    &
          & input(( ielem-1)* 3* ndofs+( idof-1)* 3+2)
      end do
    end do

  end subroutine atl_yFrom3D_getElement

  subroutine atl_zFrom3D_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the point time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nPnts*3)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = input                                    )

    do iPoint = 1, nPnts
      res((ipoint-1)*1+1) = input((ipoint-1)*3+3)
    end do

  end subroutine atl_zFrom3D_getPoint

  subroutine atl_zFrom3D_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                     nDofs, res                          )
    ! ---------------------------------------------------------------------------
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

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nDofs*nElems*3)
    integer :: iElem, iDof
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = input                                      )

    !! We only want to have one component, thus we don't need to loop over them
    do iElem = 1, nElems
      do iDof = 1, nDofs
        res(( ielem-1)* 1* ndofs+( idof-1)* 1+1) =    &
          & input(( ielem-1)* 3* ndofs+( idof-1)* 3+3)
      end do
    end do

  end subroutine atl_zFrom3D_getElement

end module atl_eqn_maxwell_derive_module
