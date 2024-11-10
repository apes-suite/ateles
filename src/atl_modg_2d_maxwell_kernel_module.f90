! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014-2016 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Parid Ndreka
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

!> Module for routines and datatypes of MOdal Discontinuous Galerkin (MODG)
!! scheme. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_2d_maxwell_kernel_module
  use env_module,                 only: rk
  use tem_aux_module,             only: tem_abort
  use tem_logging_module,         only: logUnit

  use ply_dof_module,             only: Q_space
  use ply_poly_project_module,    only: ply_poly_project_type, &
    &                                   assignment(=), &
    &                                   ply_poly_project_n2m, &
    &                                   ply_poly_project_m2n
  use ply_oversample_module,      only: ply_convert2oversample, &
    &                                   ply_convertFromoversample

  use atl_equation_module,        only: atl_equations_type
  use atl_maxwell_flux_2d_module, only: atl_maxwell_flux_2d
  use atl_modg_2d_scheme_module,  only: atl_modg_2d_scheme_type
  use atl_facedata_module,        only: atl_facedata_type
  use atl_materialPrp_module,     only: atl_material_type, &
    &                                   atl_VarMatIdx
  use atl_penalization_module,    only: atl_penalizationData_type
  use atl_scheme_module,          only: atl_scheme_type

  implicit none
  private

  public :: atl_modg_maxwell_2d_numFlux,         &
    & atl_modg_maxwell_2d_physFlux_Const,        &
    & atl_modg_maxwell_2d_physFlux_NonConst,     &
    & atl_modg_maxwell_2d_penalization_NonConst, &
    & atl_modg_maxwell_2d_penalization_Const

contains

  !> Calculate the numerical flux for Maxwell equation and MODG scheme
  subroutine atl_modg_maxwell_2d_numFlux( equation, facedata, scheme, &
    &                                     poly_proj, material         )
    ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_2d_scheme_type), intent(inout) :: scheme
    !> Data for projection method
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------
    integer :: iDir, nFaceDofs
    real(kind=rk), allocatable :: modalCoeffs(:,:,:)
    real(kind=rk), allocatable :: pntVal(:,:,:)
    real(kind=rk), allocatable :: nodalNumFlux(:,:)
    real(kind=rk), allocatable :: numFluxBuffer(:,:)
    integer :: nquadpoints, ndofs, oversamp_dofs
    ! --------------------------------------------------------------------------

    nFaceDofs = (scheme%maxPolyDegree+1)

    ! correct amount of points due to projection
    nquadpoints = poly_proj%body_1D%nquadpoints
    ndofs = poly_proj%body_1D%ndofs
    oversamp_dofs = poly_proj%body_1D%oversamp_dofs


    ! Calculate the numerical fluxes for the faces with constant material paramters
    do iDir = 1,2
      if(size(material%material_desc%computeFace(iDir,1)%leftPos) > 0) then
        call atl_maxwell_flux_2d(                                            &
          & nTotalFaces    = size(facedata%faceRep(iDir)%dat,1),             &
          & nSides         = size(material%material_desc%computeFace(iDir,1) &
          &                                             %leftPos),           &
          & nFaceDofs      = nFaceDofs,                                      &
          & faceRep        = facedata%faceRep(iDir)%dat,                     &
          & faceFlux       = facedata%faceFlux(iDir)%dat,                    &
          & leftPos        = material%material_desc%computeFace(iDir,1)      &
          &                                        %leftPos,                 &
          & rightPos       = material%material_desc%computeFace(iDir,1)      &
          &                                        %rightPos,                &
          & var            = equation%varRotation(iDir)                      &
          &                          %varTransformIndices(:),                &
          & material_left  = material%material_dat%faceMaterialData(iDir,1)  &
          &                                       %leftElemMaterialDat,      &
          & material_right = material%material_dat%faceMaterialData(iDir,1)  &
          &                                       %rightElemMaterialDat      )
      end if
    end do


    ! Calculate the numerical fluxes for the faces with non-constant material paramters
    do iDir = 1,2
      if(size(material%material_desc%computeFace(iDir,2)%leftPos) > 0) then
        ! Variable material parameters work only for Q_space polynomials
        if(.not.(scheme%basisType.eq.Q_space)) then
          call tem_abort( 'Error in modg_maxwell_numFlux: Variable material' &
            & // ' parameters work only for Q polynomial space'              )
        end if
        allocate( modalCoeffs(oversamp_dofs,equation%varSys%nScalars,2) )
        allocate( pntVal(nQuadPoints,equation%varSys%nScalars,2) )
        allocate( nodalNumFlux(nQuadPoints,equation%varSys%nScalars) )
        allocate( numFluxBuffer(oversamp_dofs,equation%varSys%nScalars) )
        call atl_maxwell_flux_2d(                                            &
          & nTotalFaces    = size(facedata%faceRep(iDir)%dat,1),             &
          & nSides         = size(material%material_desc%computeFace(iDir,2) &
          &                                             %leftPos),           &
          & nFaceDofs      = nFaceDofs,                                      &
          & faceRep        = facedata%faceRep(iDir)%dat,                     &
          & faceFlux       = facedata%faceFlux(iDir)%dat,                    &
          & leftPos        = material%material_desc%computeFace(iDir,2)      &
          &                                        %leftPos,                 &
          & rightPos       = material%material_desc%computeFace(iDir,2)      &
          &                                        %rightPos,                &
          & var            = equation%varRotation(iDir)                      &
          &                          %varTransformIndices(:),                &
          & material_left  = material%material_dat%faceMaterialData(iDir,2)  &
          &                                       %leftElemMaterialDat,      &
          & material_right = material%material_dat%faceMaterialData(iDir,2)  &
          &                                       %rightElemMaterialDat,     &
          & poly_proj      = poly_proj,                                      &
          & modalCoeffs    = modalCoeffs,                                    &
          & pntVal         = pntVal,                                         &
          & nodalNumFlux   = nodalNumFlux,                                   &
          & numFluxBuffer  = numFluxBuffer                                   )
        deallocate( modalCoeffs )
        deallocate( pntVal )
        deallocate( nodalNumFlux )
        deallocate( numFluxBuffer )
      end if
    end do


  end subroutine atl_modg_maxwell_2d_numFlux


  !> Calculate the physical flux for the MODG scheme and Maxwell equation.
  !! used for both space P and Q
  subroutine atl_modg_maxwell_2d_physFlux_const( equation, res, state, iElem,  &
                                &  iDir, penalizationData, poly_proj, material,&
                                &  nodal_data,nodal_gradData, nodal_res,       &
                                &  elemLength,  scheme_min, scheme_current     )
  ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> To store the resulting phy flux in modal form
    real(kind=rk), intent(inout)     :: res(:,:)
    !> The state of the equation
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current Element
    integer, intent(in) :: iElem
    !> The current Direction
    integer, intent(in) :: iDir
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    real(kind=rk), intent(inout)     :: nodal_res(:,:)
    !> Length of the element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! --------------------------------------------------------------------------!
    ! Rotation indices for physical flux calculation in y direction
    integer :: rot(7)
    integer :: nDofs, nScalars
    real(kind=rk) :: inv_mu, inv_epsi
    ! --------------------------------------------------------------------------!

    ! get the rotation for the physical flux calculation in y
    rot(1:7) = equation%varRotation(iDir)%varTransformIndices(1:7)

    nDofs = size(state,1)
    nScalars = equation%varSys%nScalars

    ! Calculate the physical flux for all elements with constant material
    ! parameters.
    inv_mu   = 1.0_rk &
      &         / material%material_dat%elemMaterialData(1) &
      &                                %materialDat(iElem,1,1)
    inv_epsi = 1.0_rk &
      &         / material%material_dat%elemMaterialData(1) &
      &                                %materialDat(iElem,1,2)

    call compute_physFlux_2d(                    &
      &                   nDofs     = nDofs,     &
      &                   nScalars  = nScalars,  &
      &                   state_der = res,       &
      &                   state     = state,     &
      &                   rot       = rot,       &
      &                   inv_mu    = inv_mu,    &
      &                   inv_epsi  = inv_epsi   )


  end subroutine atl_modg_maxwell_2d_physFlux_Const

  !> Calculate the physical flux for the MODG scheme and Maxwell equation.
  !! used for both space P and Q
  subroutine atl_modg_maxwell_2d_physFlux_NonConst( equation, res, state, &
    & iElem, iDir, penalizationData, poly_proj, material, nodal_data,     &
    & nodal_gradData, nodal_res, elemLength,  scheme_min, scheme_current  )
  ! --------------------------------------------------------------------------
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> To store the resulting phy flux in modal form
    real(kind=rk), intent(inout)     :: res(:,:)
    !> The state of the equation
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current Element
    integer, intent(in) :: iElem
    !> The current Direction
    integer, intent(in) :: iDir
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    real(kind=rk), intent(inout)     :: nodal_res(:,:)
    !> Length of the element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! --------------------------------------------------------------------------!
    ! Rotation indices for physical flux calculation in y direction
    integer :: rot(7)
    integer :: nDofs, nScalars
    ! --------------------------------------------------------------------------!

    ! get the rotation for the physical flux calculation in y
    rot(1:7) = equation%varRotation(iDir)%varTransformIndices(1:7)

    !nTotalElems = size(statedata%state,1)
    nDofs = size(state,1)
    nScalars = equation%varSys%nScalars

    ! Calculate the physical flux for all elements with constant material
    ! parameters.
    ! We do it for the fluxes in all three spatial directions.
    ! This is achieved by using a properly defined variable rotation.
    ! Variable material parameters work only for Q_space polynomials
    if(.not.(scheme_current%modg_2d%basisType.eq.Q_space)) then
      write(logUnit(1),*) 'Error in modg_maxwell_physFlux: Variable '
      write(logUnit(1),*) 'material parameters work only for '
      write(logUnit(1),*) 'Q polynomial space, stopping ...'
      call tem_abort()
    end if

    call compute_physFlux_nonConst_2d(                                       &
      &   nDofs            = nDofs,                                          &
      &   nScalars         = nScalars,                                       &
      &   state_der        = res,                                            &
      &   state            = state,                                          &
      &   rot              = rot,                                            &
      &   nElems           = material%material_desc%computeElems(2)%nElems,  &
      &   elems            = material%material_desc%computeElems(2)          &
      &                              %totElemIndices ,                       &
      &   material         = material%material_dat%elemMaterialData(2)       &
      &                              %materialDat,                           &
      &   poly_proj        = poly_proj,                                      &
      &   modalCoeffs      = scheme_min%temp_over,                           &
      &   iElem            = iElem ,                                         &
      &   nodalPhysFlux    = scheme_min%temp_nodal                           )


  end subroutine atl_modg_maxwell_2d_physFlux_NonConst


  !> Compute the physical flux in x direction.
  !!
  !! For other directions a properly defined variable permutation can be used.
  !! This routine covers only constant material parameters.
  subroutine compute_physFlux_2d( nDofs, nScalars, state_der, state, rot, &
      &                           inv_mu, inv_epsi                        )
    ! ---------------------------------------------------------------------------
    !> dimensions
    integer, intent(in) :: nDofs, nScalars
    !> Array to store the fluxes in.
    real(kind=rk), intent(inout) :: state_der(:,:)
    !> State to compute the fluxes from.
    real(kind=rk), intent(in) :: state(nDofs,nScalars)
    !> Rotationing to index the variables.
    integer, intent(in) :: rot(7)
    real(kind=rk), intent(in) :: inv_mu, inv_epsi
    ! ---------------------------------------------------------------------------
    integer :: iDoF
    ! ---------------------------------------------------------------------------

    do iDoF=1,nDoFs
        ! D_x
        state_der(iDoF,rot(1)) = 0.0_rk
        ! PML
        state_der(iDoF,rot(4:7)) = 0.0_rk
        ! D_y
        state_der(iDoF,rot(2)) = state(iDoF,rot(3)) * inv_mu
        ! B_z
        state_der(iDoF,rot(3)) = state(iDoF,rot(2)) * inv_epsi
    end do

  end subroutine compute_physFlux_2d


  !> Function for physical flux of the Maxwell equations in terms of D and B.
  function atl_physFluxMaxwell_2d(state, material) result(flux)
    ! ---------------------------------------------------------------------------
    !> State to compute the fluxes from (D,B).
    real(kind=rk), intent(in) :: state(7)
    !> Material parameters (mu, epsilon) the flux calculation
    real(kind=rk), intent(in) :: material(3)
    !> The resulting flux in x direction
    real(kind=rk) :: flux(7)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: inv_mu, inv_epsi
    ! ---------------------------------------------------------------------------

    ! Get magnetic permeability
    inv_mu = 1.0_rk / material(1)
    ! Get electric permitivity
    inv_epsi = 1.0_rk / material(2)

    ! D_x
    flux(1) = 0.0_rk

    ! PML
    flux(4:7) = 0.0_rk

    ! D_y
    flux(2) = state(3) * inv_mu

    ! B_z
    flux(3) = state(2) * inv_epsi

  end function atl_physFluxMaxwell_2d



  !> Compute the physical flux in x direction.
  !!
  !! For other directions a properly defined variable permutation can be used.
  !! This routine covers non-constant material parameters.
  subroutine compute_physFlux_nonConst_2d( nDofs, nScalars, iElem, state_der, &
    &                                      state, rot, nElems, elems,         &
    &                                      material, poly_proj, modalCoeffs,  &
    &                                      nodalPhysFlux )
    ! ---------------------------------------------------------------------------
    !> dimensions
    integer, intent(in) ::  nDofs, nScalars
    !> Array to store the fluxes in.
    real(kind=rk), intent(inout) :: state_der(:,:)
    !> State to compute the fluxes from.
    real(kind=rk), intent(in) :: state(nDofs,nScalars)
    !> Rotationing to index the variables.
    !integer, intent(in) :: rotY(7)
    integer, intent(in) :: rot(7)
    !> Number of elements.
    integer, intent(in) :: nElems
    integer, intent(in) :: iElem
    !> Element positions in the total state vector
    integer, intent(in) :: elems(nElems)
    !> Material parameters (mu, epsilon) for all elements
    real(kind=rk), intent(in) :: material(nElems,nDofs,3)
    !> Data for projection method
    type(ply_poly_project_type) :: poly_proj
    !> Working array for modal coefficients of the current element in the loop.
    real(kind=rk), intent(inout) :: modalCoeffs(poly_proj%body_2D%oversamp_dofs,size(state,2),1)
    !> Working array for nodal representation of the physical flux along the 3 spatial directions.
    real(kind=rk), intent(inout) :: nodalPhysFlux(poly_proj%body_2D%nquadpoints,size(state,2),2)
    ! ---------------------------------------------------------------------------
    integer :: iPoint, glob_elem_ind
    ! ---------------------------------------------------------------------------
    glob_elem_ind = elems(iElem)


    ! get the modal coefficients of the current cell (for all variables
    ! of the Maxwell equation, therefore we use ":" for the third index).
    ! ATTENTION: have to be duplicated as the FPT is modifying the input vector.

    call ply_convert2oversample(state       = state(:,:),        &
      &                         poly_proj   = poly_proj,         &
      &                         nDim        = 2,                 &
      &                         modalCoeffs = modalCoeffs(:,:,1) )

    ! Now, we transform the modal representation of this element to nodal
    ! space by making use of fast polynomial transformations (FPT)
    call ply_poly_project_m2n(me = poly_proj,         &
      &                       dim = 2 ,               &
      &                       nVars = 3,              &
      !&                       nVars = nScalars,              &
      &                       nodal_data= nodalPhysFlux(:,:,2),   &
      &                       modal_data= modalCoeffs(:,:,1) )
    ! --> oversampling nodal space

    ! ... the PML does not involve physical fluxes
    !pointVal(:,4:7) = 0.0_rk
    nodalPhysFlux(:,4:7,2) = 0.0_rk

    ! Calculate the physical flux point by point within this cell - x direction
    do iPoint = 1, poly_proj%body_2D%nquadpoints
      !nodalPhysFlux(iPoint,:) = atl_physFluxMaxwell_2d( pointVal(iPoint,:),      &
      !                                             & material(iElem,iPoint,:) )
      nodalPhysFlux(iPoint,rot,1) = atl_physFluxMaxwell_2d(                &
        &                                nodalPhysFlux(iPoint,rot,2),      &
                                         & material(iElem,iPoint,:)        )
    end do
    ! Transform the nodal physical flux back to modal space
    call ply_poly_project_n2m(me = poly_proj,                   &
      &                       dim = 2,                          &
      &                       nVars = 3,                        &
      &                       nodal_data= nodalPhysFlux(:,:,1), &
      &                       modal_data= modalCoeffs(:,:,1)    )
    ! ... the PML does not involve physical fluxes

    modalCoeffs(:,4:7,1) = 0.0_rk


    call ply_convertFromOversample(modalCoeffs = modalCoeffs(:,:,1), &
     &                             poly_proj= poly_proj,             &
     &                             nDim = 2,                         &
     &                             state = state_der(:,:)            )
    ! --> oversampl modal space

  end subroutine compute_physFlux_nonConst_2d



  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_maxwell_2d_penalization_NonConst(equation, poly_proj, &
    &         nodal_data, scheme_min, penalizationData, iElem, material     )
    ! --------------------------------------------------------------------------!
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The current Element
    integer, intent(in) :: iElem
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------!
    integer :: iPoint,  glob_elem_ind
!PV!    real(kind=rk) :: mat(material%material_desc%computeElems(2)   &
!PV!     &                              %nElems, poly_proj%body_2D%nquadpoints,3)
    real(kind=rk), allocatable :: mat(:,:,:)
    !real(kind=rk) :: temp(poly_proj%body_2D%nquadpoints,7)
    ! --------------------------------------------------------------------------!
    !> @todo PV 20150820 Get the correct penalization data here
    mat = material%material_dat%elemMaterialData(atl_VarMatIdx)%materialDat
    glob_elem_ind = material%material_desc               &
      &                     %computeElems(atl_VarMatIdx) &
      &                     %totElemIndices(iElem)


!    temp = nodal_data
!    temp(:,4:7) = 0.0_rk


    do iPoint = 1,poly_proj%body_2D%nquadpoints
      scheme_min%temp_nodal(iPoint,1:2,1) = (-1.0_rk) * mat(iElem,iPoint,3) &
                  & * nodal_data(iPoint,1:2) / mat(iElem,iPoint,2)
    end do

    call ply_poly_project_n2m( me         = poly_proj,                     &
      &                        dim        = 2,                             &
      &                        nVars      = 2,                             &
      &                        nodal_data = scheme_min%temp_nodal(:,:2,1), &
      &                        modal_data = scheme_min%temp_over(:,:2,1)   )

    ! ... lossy material has no impact on B_z and the remaining PML variables.
    scheme_min%temp_over(:,3:7,1) = 0.0_rk

    ! --> oversampl modal space
    call ply_convertFromOversample(                                          &
      & modalCoeffs = scheme_min%temp_over(:,:2,1),                          &
      & poly_proj   = poly_proj,                                             &
      & nDim        = 2,                                                     &
      & state       = penalizationdata%penalization_data(glob_elem_ind,:,:2) )


  end  subroutine atl_modg_maxwell_2d_penalization_NonConst

  ! Calculate the penalization terms (for density, momentum, energy)
  ! The penalization terms are calculated in the sammer manner as
  ! the physical fluxes, i.e. in a nodal-modal approach
  subroutine atl_modg_maxwell_2d_penalization_Const(equation, poly_proj,     &
    &         nodal_data, scheme_min, penalizationData, iElem, material     )
    ! --------------------------------------------------------------------------!
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> Poly project
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The state in nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    !> The scheme information
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> The penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The current Element
    integer, intent(in) :: iElem
    !> Material description for the faces on the current level
    type(atl_material_type), intent(inout) :: material
    ! --------------------------------------------------------------------------!


    write(logUnit(7),*) " Penalization term for material " &
      &                 //"with constant material parameter."
    write(logUnit(7),*) " Not yet implemented "

  end  subroutine atl_modg_maxwell_2d_penalization_Const
end module atl_modg_2d_maxwell_kernel_module

