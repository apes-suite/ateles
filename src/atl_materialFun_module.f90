! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2019 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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

module atl_materialFun_module

  use env_module,               only: rk, labelLen

  ! Aotus modules
  use aot_out_module,           only: aot_out_type,        &
    &                                 aot_out_open_table,  &
    &                                 aot_out_close_table, &
    &                                 aot_out_val

  use tem_varSys_module,        only: tem_varSys_type

  implicit none
  private

  !> @todo PV 20150827 Rename atl_materialFun_type to atl_material_type to stay
  !!                   consistent with atl_source_type

  !> Contains settings for mode reduction. Mode reduction is used to decrease
  !! computational effort for elements where a calculation with higher modes
  !! doesn"t increase the result's resolution, because material is equal over
  !! the element's volume and faces.
  type atl_mode_reduction_type
    !> Indicates whether the mode reduction is enabled at all.
    logical :: enabled = .false.
    !> The threshold to decide whether a reduction is possible. Values of the
    !! material at matPos above threshold indicate a possible reduction.
    real(kind=rk) :: threshold = 1.0_rk
    !> The position, if present,  of the material within the materials that is
    !! used to determine whether mode reduction is possible or not.
    !! -1 if not available.
    integer :: matPos
  end type atl_mode_reduction_type


  !> Generic description of material property functions.
  type atl_materialFun_type
    !> The number of material parameters used by the equation system.
    integer :: nMat

    !> Number of scalar entries to describe the material.
    integer, allocatable :: nScalars(:)

    !> The indizes of variables in the global varSys that are used as
    !! material's, penalization's or whatever's data sources.
    !!
    !! This array is initialized while reading the present mappings from the
    !! configuration.
    integer, allocatable :: matvar_pos(:)

    !> Flag to indicate whether the material may be varying over time and
    !! thus requires updates over time.
    logical :: is_transient

    !> The names of the material parameter.
    !!
    !! This can be used together with nScalars to find the position of the
    !! different material parameter components in the material data arrays.
    !! Attention: We store the prefixed name of the material parameter, not the
    !! name of the material variable in the configuration's variable table. This
    !! name can be achieved via the input_varName stored in the variable system.
    character(labelLen), allocatable :: matParNames(:)

    !> Settings for mode reduction, i.e. reducing the modes where possible to
    !! decrease computational effort.
    type(atl_mode_reduction_type) :: mode_reduction

  end type atl_materialFun_type

  public :: atl_materialFun_type, &
    & atl_mode_reduction_type,    &
    & atl_dump_materialFun


contains

  ! ****************************************************************************
  !> Strips the prefix from material variables created by the solver when
  !! adding material to the variable system.
  function stripMatPrefix( materialName ) result( strippedName )
    ! --------------------------------------------------------------------------
    character(labelLen), intent(in) :: materialName
    character(labelLen) strippedName
    ! --------------------------------------------------------------------------

    strippedName = materialName(5:)

  end function
  ! ****************************************************************************


  ! ****************************************************************************
  !> Dump material description to lua file
  subroutine atl_dump_materialFun( conf, material, key, varSys )
    ! --------------------------------------------------------------------------
    type(aot_out_type) :: conf
    type(atl_materialFun_type), intent(in) :: material
    character(len=*), intent(in) :: key
    !> The variable system that contains the material variables
    type(tem_varSys_type), intent(in) :: varSys
    ! --------------------------------------------------------------------------
    integer :: ii, matVarPos, luaVarPos
    ! --------------------------------------------------------------------------

    call aot_out_open_table( put_conf = conf, tname = key )

    do ii = 1, size(material%matParNames)
      matVarPos = material%matvar_pos(ii)
      ! Material variables have exactly one input variable, so we can safelt
      ! access the first and only the first index here.
      luaVarPos = varSys%method%val(matVarPos)%input_varPos(1)
      call aot_out_val(                                        &
        & put_conf = conf,                                     &
        & vname    = stripMatPrefix(material%matParNames(ii)), &
        & val      = varSys%varname%val(luaVarPos)             )
    end do

    ! Close material table
    call aot_out_close_table( put_conf = conf )

  end subroutine atl_dump_materialFun
  ! ****************************************************************************


end module atl_materialFun_module
