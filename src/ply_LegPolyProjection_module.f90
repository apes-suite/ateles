! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
!
! Parts of this file were written by Daniel Fleischer, Daniel Petró, Peter Vitt
! and Raphael Haupt for University of Siegen.
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
!> Module for projection of Q Legendre Polynomials from parent cell
!! to child cells.
!!
module ply_LegPolyProjection_module

  ! include treelm modules
  use env_module,          only: rk
  use tem_param_module,    only: PI
  use tem_aux_module,      only: tem_abort
  use treelmesh_module,    only: treelmesh_type
  use tem_param_module,    only: childPosition

  implicit none

  private

  public :: ply_gauleg

  ! ************************************************************************ !
  !> Parameter to specify Legendre polynomials as the degrees of freedoms
  !! of the elements. The multidimensional polynomias are build as
  !! Q-polynomials.
  !! The projection is a L2-Projection onto the ansatz function of the finer
  !! elements.
  integer, parameter :: ply_QLegendrePoly_prp = 1
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Datatype storing the coefficients arising for the projection
  !! of solutions on a parent cell to its children during the subsampling
  !! routines.
  type ply_ProjCoeff_type
    !> Array holding all the projection coefficients for the projection
    !! of a degree of freedom on the parent element (first index) to
    !! a degree of freedom on the child element (second index) for
    !! a given child element (third index).
    !! Therefore the dimension of this array is (nDofs, nDofs, 8).
    real(kind=rk), allocatable :: projCoeff(:,:,:)
  end type ply_ProjCoeff_type
  ! ************************************************************************ !


  ! ************************************************************************ !
  type ply_subsample_type
    !> Is subsampling active
    logical :: isActive = .false.

    !> The current sampling lvl.
    integer :: sampling_lvl

    !> The type of projection we use to subsample the elemental data.
    integer :: projectionType = ply_QLegendrePoly_prp

    !> Maximal Level down to which subsampling should be done.
    integer :: caplevel = 20

    !> Minimal subsampling depth:
    integer :: minsub = 0

    !> Maximal subsampling depth:
    integer :: maxsub = 0

    !> Maximum allowed oscillation of the solution.
    !! For adaptive subsampling only.
    real(kind=rk) :: eps_osci

    !> Factor for the reduction of the degrees of freedom in one subsampling
    !! step (per spatial direction).
    real(kind=rk) :: dofReducFactor

    !> Indicator for limitation of total memory consumption
    logical :: adaptiveDofReduction

    !> Absolute upper bound level to refine to.
    integer :: AbsUpperBoundLevel
  end type ply_subsample_type
  ! ************************************************************************ !

  ! ************************************************************************ !
  type ply_array_type
    real(kind=rk),allocatable :: dat(:)
  end type ply_array_type
  ! ************************************************************************ !

  public :: ply_subsample_type
  public :: ply_QPolyProjection
  public :: ply_array_type

contains

  ! ************************************************************************ !
  !> Subsampling by L2-Projection of the Q-Tensorproduct Legendre polynomials.
  subroutine ply_QPolyProjection( subsamp, dofReduction, tree, meshData,   &
    &                             varDofs, ndims, varcomps, refine_tree,   &
    &                             new_refine_tree, newMeshData, newVarDofs )
    ! -------------------------------------------------------------------- !
    !> Parameters for the subsampling.
    type(ply_subsample_type), intent(in) :: subsamp

    !> Factor for reducion of degrees of freedom.
    real(kind=rk), intent(in) :: dofReduction(:)

    !> The tree the data is written for.
    type(treelmesh_type), intent(in) :: tree

    !> The data to sub-sample.
    type(ply_array_type), intent(in) :: meshData(:)

    !> The number of degrees of freedom for each scalar variable.
    integer, intent(in) :: varDofs(:)

    !> Number of dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    !> Number of components
    integer, intent(in) :: varcomps(:)

    !> Logical array that marks elements for refinement from the last sampling
    ! level.
    logical, intent(in) :: refine_tree(:)

    !> Logical array that marks elements for refinment.
    logical, intent(in) :: new_refine_tree(:)

    !> The subsampled data.
    type(ply_array_type), allocatable, intent(out) :: newMeshData(:)

    !> The number of dofs for the subsampled dofs.
    integer, allocatable, intent(out) :: newVarDofs(:)
    ! -------------------------------------------------------------------- !
    integer :: nChilds, nVars, nDofs, nComponents
    integer :: iVar
    type(ply_ProjCoeff_type) :: projection
    type(ply_ProjCoeff_type) :: projection_oneDof
    ! Working tree and working data
    real(kind=rk), allocatable :: workDat(:)
    real(kind=rk), allocatable :: newWorkDat(:)
    integer :: nChildDofs, oneDof
    ! -------------------------------------------------------------------- !
    if (subsamp%projectionType.ne.ply_QLegendrePoly_prp) then
      call tem_abort( 'ERROR in ply_QPolyProjection: subsampling is ' &
        & // 'only implemented for Q-Legendre-Polynomials'            )
    end if

    nVars = size(varDofs)
    allocate(newVarDofs(nVars))
    allocate(newMeshData(nVars))

    varLoop: do iVar=1,nVars
      nDofs = varDofs(iVar)
      nComponents = varcomps(iVar)

      allocate(workDat(size(meshData(iVar)%dat)))
      workDat = meshData(iVar)%dat

      ! All elements in refine_tree marked with .TRUE. will get
      ! a projection to the next sampling level while all other elements
      ! will get their integral mean value ( projection till there is only
      ! one dof left )

      ! Reduce the number of dofs per direction in each subsample
      nChildDofs = (ceiling(nint(nDofs**(1.0_rk/real(ndims, kind=rk))) &
        &                      * dofReduction(iVar)))**nDims

      if (nChildDofs > nDofs) then
        nChildDofs = nDofs
      elseif (nChildDofs < 1) then
        nChildDofs = 1
      end if

      if (subsamp%sampling_lvl == subsamp%maxsub) then
        nChildDofs = 1
      end if

      ! Set the correct number of Childs for the projection to reduced Dofs.
      nChilds = 2**ndims

      ! init the projection coefficients for the current number of child dofs.
      call ply_initQLegProjCoeff( doftype    = subsamp%projectionType, &
        &                         nDofs      = nDofs,                  &
        &                         ndims      = ndims,                  &
        &                         nChilds    = nChilds,                &
        &                         nChildDofs = nChildDofs,             &
        &                         projection = projection              )

      ! init the projection coefficients for the reduction to polynomial
      ! degree of 0.
      oneDof = 1
      ! Set the correct number of childs for the projection to one dof.
      nChilds = 1
      call ply_initQLegProjCoeff( doftype    = subsamp%projectionType, &
        &                         ndofs      = nDofs,                  &
        &                         ndims      = ndims,                  &
        &                         nChilds    = nChilds,                &
        &                         nChildDofs = oneDof,                 &
        &                         projection = projection_oneDof       )

      ! ... subsample the data
      call ply_subsampleData( tree              = tree,              &
        &                     meshData          = workDat,           &
        &                     nDofs             = nDofs,             &
        &                     nChildDofs        = nChildDofs,        &
        &                     nComponents       = nComponents,       &
        &                     projection        = projection,        &
        &                     projection_oneDof = projection_oneDof, &
        &                     refine_tree       = refine_tree,       &
        &                     new_refine_tree   = new_refine_tree,   &
        &                     ndims             = ndims,             &
        &                     subsamp           = subsamp,           &
        &                     newMeshData       = newWorkDat         )

      allocate(newMeshData(iVar)%dat(size(newWorkDat)))
      newMeshData(iVar)%dat = newWorkDat

      newVarDofs(iVar) = nChildDofs

      deallocate(workDat)
      deallocate(projection%projCoeff)
      deallocate(projection_oneDof%projCoeff)

    end do varLoop

  end subroutine ply_QPolyProjection
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to initialize the projection coefficients for a usage in the
  !! subsampling routine to project degrees of freedoms of a parent cell
  !! to the degrees of freedoms of a child cell if the degrees of
  !! freedoms are Q-Legendre polynomials.
  subroutine ply_initQLegProjCoeff( doftype, nDofs, ndims, nChilds, &
    &                               nChildDofs, projection          )
    ! -------------------------------------------------------------------- !
    !> The type of degrees of freedom we have in our cells.
    integer, intent(in) :: dofType

    !> The number of degrees of freedom for the parent cells.
    integer, intent(in) :: nDofs

    !> The  number of dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    !> The number of child cells.
    integer, intent(in) :: nChilds

    !> The number of degrees of freedom for the child cells.
    integer, intent(in) :: nChildDofs

    !> The subsampling coefficients that will be initialized by this routine.
    type(ply_ProjCoeff_type), intent(out) :: projection
    ! -------------------------------------------------------------------- !
    integer :: iParentDof, iChildDof, iChild
    integer :: xShift, yShift, zShift
    integer :: xParentAnsFunc, yParentAnsFunc, zParentAnsFunc
    integer :: xChildAnsFunc, yChildAnsFunc, zChildAnsFunc
    ! First index is Legendre polynomial on the parent element, second index
    ! is the Legendre polynomial on the child element, third index is left
    ! or right projection.
    real(kind=rk), allocatable :: projCoeffOneDim(:,:,:)
    real(kind=rk) :: dimexp
    ! -------------------------------------------------------------------- !
    select case(dofType)
    case(ply_QLegendrePoly_prp)
      allocate(projection%projCoeff(nDofs, nChildDofs, nChilds))
      projection%projCoeff = 0.0_rk

      ! Create projection of one-dimensional Legendre polynomials for the
      ! reference interval [-1,+1]  (child-element) and the double length ref
      ! element [-1,+3] (parent element).
      dimexp = 1.0_rk/real(ndims, kind=rk)
      projCoeffOneDim = ply_QLegOneDimCoeff( nint(nDofs**dimexp),     &
        &                                    nint(nChildDofs**dimexp) )

      ! Loop over the children of this element
      childLoop: do iChild = 1, nChilds

        ! get the right index for the x,y,z shift of the current child with
        ! respect to the parent cell.
        if(childPosition(iChild,1).eq.-1) then
          xShift = 1
        else
          xShift = 2
        end if
        if(childPosition(iChild,2).eq.-1) then
          yShift = 1
        else
          yShift = 2
        end if
        if(childPosition(iChild,3).eq.-1) then
          zShift = 1
        else
          zShift = 2
        end if
        ! Loop over the parent dofs and calculate the projection coefficients
        ! for each of the child dofs
        parentDofLoop: do iParentDof = 1, nDofs
          ! convert the number of the parent dof to ansatz function numbers in
          ! the spatial direction.
          call ply_dofToQPoly( dof      = iParentDof,     &
            &                  nDofs    = nDofs,          &
            &                  ndims    = ndims,          &
            &                  xAnsFunc = xParentAnsFunc, &
            &                  yAnsFunc = yParentAnsFunc, &
            &                  zAnsFunc = zParentAnsFunc  )

          childDofLoop: do iChildDof = 1, nChildDofs
            ! convert the number of the child dof to ansatz function numbers in
            ! the spatial direction.
            call ply_dofToQPoly( dof      = iChildDof,     &
              &                  nDofs    = nChildDofs,    &
              &                  ndims    = ndims,         &
              &                  xAnsFunc = xChildAnsFunc, &
              &                  yAnsFunc = yChildAnsFunc, &
              &                  zAnsFunc = zChildAnsFunc  )

            ! reuse the one-dimensional projection coefficients to build up the
            ! 3D
            projection%projCoeff(iParentDof, iChildDof, iChild)            &
              & = projCoeffOneDim(xParentAnsFunc, xChildAnsFunc, xShift)   &
              &   * projCoeffOneDim(yParentAnsFunc, yChildAnsFunc, yShift) &
              &   * projCoeffOneDim(zParentAnsFunc, zChildAnsFunc, zShift)

          end do childDofLoop
        end do parentDofLoop
      end do childLoop

    case default
      call tem_abort( 'ERROR in ply_initProjCoeff: initialization of '      &
        & // 'projection coefficients for subsampling is implemented only ' &
        & // 'for Q-Legendre polynomials'                                   )
    end select
    deallocate(projCoeffOneDim)
  end subroutine ply_initQLegProjCoeff
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to create one-dimensional projection coefficient for a coarse
  !! element to a fine element.
  function ply_QLegOneDimCoeff( nDofsOneDim, nChildDofsOneDim ) &
    & result(projCoeffOneDim)
    ! -------------------------------------------------------------------- !
    !> The number of dofs in one dimension.
    integer , intent(in) :: nDofsOneDim

    !> The number of dofs in one dimension for the children.
    integer , intent(in) :: nChildDofsOneDim

    !> Projected one-dimensional coefficients.
    !!
    !! First index is Legendre polynomial on the parent element, second index
    !! is the Legendre polynomial on the child element, third index is left
    !! or right projection.
    real(kind=rk), allocatable :: projCoeffOneDim(:,:,:)
    ! -------------------------------------------------------------------- !
    integer :: nIntP, iParentFunc, iChildFunc
    real(kind=rk), allocatable :: points(:), weights(:),              &
      &                           pointsLeft(:), pointsRight(:),      &
      &                           parentFuncVal(:,:), childFuncVal(:,:)
    ! -------------------------------------------------------------------- !
    allocate(projCoeffOneDim(nDofsOneDim, nChildDofsOneDim,2))

    ! Create the gauss legendre quadrature points for the reference element
    ! [-1,+1]
    nIntP = ceiling(( max(nDofsOneDim,nChildDofsOneDim) )/2.0)**2
    call ply_gauleg(-1.0_rk, +1.0_rk, points, weights, nIntP)

    ! Now, project ansatz function of the parent for the left child.
    ! We apply a Gaussian quadrature and take case of composition of
    ! the mapping from reference to physical element (for child and
    ! parent).
    ! So, we evaluate the parent Legendre polynomial at shifted quadrature
    ! points ( shifted by 0.5 x - 0.5 ) and evaluate the child Legendre
    ! quadrature points at the original Gauss-Legendre points.
    allocate( pointsLeft(size(points)) )
    pointsLeft(:) = 0.5_rk * points(:) - 0.5_rk
    parentFuncVal = ply_legVal( pointsLeft, nIntP, nDofsOneDim-1 )
    childFuncVal = ply_legVal( points, nIntP, nChildDofsOneDim-1 )
    do iParentFunc = 1, nDofsOneDim
      do iChildFunc = 1, nChildDofsOneDim
        ! ... calculate the integral by the quadrature rule.
        projCoeffOneDim(iParentFunc, iChildFunc,1)           &
          & = sum( weights(:) * parentFuncVal(iParentFunc,:) &
          &        * childFuncVal(iChildFunc,:) )
        ! ... and normalize by the norm of the child function.
        projCoeffOneDim(iParentFunc, iChildFunc,1)       &
          & = projCoeffOneDim(iParentFunc, iChildFunc,1) &
          &   * (1.0_rk / ply_QLegSqNorm(iChildFunc))
      end do
    end do
    deallocate(parentFuncVal)
    deallocate(childFuncVal)
    deallocate(pointsLeft)

    ! now, project ansatz function of the parent for the right child. We apply
    ! the same procedure as for the left child, but with a different shift of
    ! the quadrautre point for the parent function.
    allocate( pointsRight(size(points)) )
    pointsRight(:) = 0.5_rk * points(:) + 0.5_rk
    parentFuncVal = ply_legVal( pointsRight, nIntP, nDofsOneDim-1 )
    childFuncVal = ply_legVal( points, nIntP, nChildDofsOneDim-1 )
    do iParentFunc = 1, nDofsOneDim
      do iChildFunc = 1, nChildDofsOneDim
        ! ... calculate the integral by the quadrature rule.
        projCoeffOneDim(iParentFunc, iChildFunc,2)           &
          & = sum( weights(:) * parentFuncVal(iParentFunc,:) &
          &        * childFuncVal(iChildFunc,:) )
        ! ... and normalize by the norm of the child function.
        projCoeffOneDim(iParentFunc, iChildFunc,2)       &
          & = projCoeffOneDim(iParentFunc, iChildFunc,2) &
          &   * (1.0_rk / ply_QLegSqNorm(iChildFunc))
      end do
    end do
    deallocate(parentFuncVal)
    deallocate(childFuncVal)
    deallocate(pointsRight)

  end function ply_QLegOneDimCoeff
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function to calculate the squared L2-Norm of a given Legendre polynomial
  !! on the reference element [-1,+1].
  function ply_QLegSqNorm( polyIndex ) result(sqNorm)
    ! -------------------------------------------------------------------- !
    !> The Legendre polynomial index to calculate the squared norm for.
    !! The first polynomial has index 1.
    integer, intent(in) :: polyIndex

    !> The squared L2 Norm of the Legendre polynomial.
    real(kind=rk) :: sqNorm
    ! -------------------------------------------------------------------- !

    sqNorm = 2.0_rk / ( 2.0_rk * polyIndex  - 1.0_rk)

  end function ply_QLegSqNorm
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Evaluate a given set of Legendre polynomials a given set of 1D points.
  !!
  function ply_legVal( points, nPoints, maxPolyDegree ) result( val )
    ! -------------------------------------------------------------------- !
    !> A given set of 1D points.
    real(kind=rk), intent(in) :: points(:)

    !> The number of points to evaluate the polynomials at.
    integer, intent(in) :: nPoints

    !> The maximal polynomial degree to evaluate for.
    integer, intent(in) :: maxPolyDegree

    !> Function values for for all Legendre polynomials up to degree
    !! maxPolyDegree at all given points.
    !! Therefore the dimension of this array is (maxPolyDegree+1, nPoints)
    real(kind=rk), allocatable :: val(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iAns
    ! -------------------------------------------------------------------- !

    allocate(val(maxPolyDegree+1, nPoints))

    ! the first Legendre polynomial is constant.
    val(1,:) = 1.0_rk
    if(maxPolyDegree > 0) then
      ! the second Legendre ploynomial is identity.
      val(2,:) = points(:)
      ! apply the recursion
      do iAns = 3, maxPolyDegree+1
        val(iAns,:) = ( (2.0_rk*(iAns-1.0_rk)-1.0_rk) &
          &               * points(:)*val(iAns-1,:)   &
          &             - ((iAns-1)-1)                &
          &               * val(iAns-2,:) )           &
          &           / (iAns-1)
      end do
    end if

  end function ply_legVal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to convert linearized dof index to ansatz function number for
  !! Q-Polynomials.
  subroutine ply_dofToQPoly( dof, nDofs, ndims, xAnsFunc, yAnsFunc, zAnsFunc )
    ! -------------------------------------------------------------------- !
    !> The linearized degree of freedom index.
    integer, intent(in) :: dof

    !> The number of dofs for all directions.
    integer, intent(in) :: nDofs

    !> The number of Dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    !> The ansatz function number in x direction.
    integer, intent(out) :: xAnsFunc

    !> The ansatz function number in y direction.
    integer, intent(out) :: yAnsFunc

    !> The ansatz function number in z direction.
    integer, intent(out) :: zAnsFunc
    ! -------------------------------------------------------------------- !
    integer :: nDofsPerDir
    ! -------------------------------------------------------------------- !

    ! The number of dofs in each direction
    nDofsPerDir = nint(nDofs**(1.0_rk/real(ndims,kind=rk)))

    ! now, we compute the ansatz function numbers
    ! Works for polynomials of all dimensionality, as dof will never reach
    ! sufficiently high values for lower dimensions.
    zAnsFunc = (dof-1) / (nDofsPerDir**2)  + 1
    yAnsFunc = (dof-1-(zAnsFunc-1)*(nDofsPerDir**2)) / nDofsPerDir + 1
    xAnsFunc = dof - (zAnsFunc-1)*(nDofsPerDir**2) - (yAnsFunc-1)*nDofsPerDir

  end subroutine ply_dofToQPoly
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Routine to subsample mesh information for one refinement level.
  subroutine ply_subsampleData( tree, meshData, nDofs, nChildDofs,          &
    &                           nComponents, projection, projection_oneDof, &
    &                           refine_tree, new_refine_tree, ndims,        &
    &                           subsamp, newMeshData                        )
    ! -------------------------------------------------------------------- !
    !> The tree the data is written for.
    type(treelmesh_type), intent(in) :: tree

    !> The data to sub-sample.
    real(kind=rk), intent(in) :: meshData(:)

    !> The number of degrees of freedom for each scalar variable.
    integer, intent(in) :: nDofs

    !> The number of degrees of freedom per scalar variable on the child
    !! elements.
    integer, intent(in) :: nChildDofs

    !> Number of components
    integer, intent(in) :: nComponents

    !> Projection coefficients for the given data.
    type(ply_ProjCoeff_type), intent(in) :: projection

    !> Projection coeffiecients for the the reduction to polynomial
    !! degree of 0.
    type(ply_ProjCoeff_type), intent(in), optional :: projection_oneDof

    !> Logical array that marks all elements for refinement from the last
    !! sampling lvl.
    logical, intent(in), optional :: refine_tree(:)

    !> Logical array that marks all elements for refinement
    logical, intent(in) :: new_refine_tree(:)

    !> The number of dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    !> Parameters for the subsampling.
    type(ply_subsample_type), intent(in) :: subsamp

    !> The subsampled data.
    real(kind=rk), allocatable, intent(out) :: newMeshData(:)
    ! -------------------------------------------------------------------- !
    integer :: nElems, nChilds, nParentElems, nElemsToRefine, nElemsNotToRefine
    integer :: iElem, iParentElem, iChild
    integer :: lowElemIndex, upElemIndex, lowChildIndex, upChildIndex
    integer :: oneDof, noChilds, childpos
    real(kind=rk), allocatable :: childData(:)
    ! -------------------------------------------------------------------- !
    nChilds = 2**ndims
    nElems = tree%nElems
    nElemsToRefine = count(new_refine_tree)
    nElemsNotToRefine = nElems - nElemsToRefine
    nParentElems = size(refine_tree)

    ! Now, we set the correct data for the newMeshData.
    allocate(newMeshData((nElemsToRefine * nChilds * nChildDofs &
      &                  + nElemsNotToRefine) * nComponents))
    newMeshData = 0.0_rk
    upChildIndex = 0
    upElemIndex = 0
    childPos = 0

    if (subsamp%sampling_lvl > 1) then

      elementLoop: do iParentElem=1,nParentElems
        ! Check if the parent cell was already refined...
        if (refine_tree(iParentElem)) then
          ! Parent cell was already refined so it contains data with nDofs
          childLoop: do iChild=1,nChilds
            childPos = childPos + 1
            ! Check if the child elems will be refined...
            if (new_refine_tree(childPos)) then

              ! Child cell will be refined.
              !
              ! Need to project current elem data to new childs with reduced
              ! dofs.
              allocate(childData(nChildDofs*nChilds*nComponents))
              ! Create lower and upper indices for all data of iElem in
              ! meshData.
              lowElemIndex = upElemIndex + 1
              upElemIndex = (lowElemIndex-1) + nDofs * nComponents

              ! Project these dofs from the coarse element to the
              ! finer elements.
              call ply_projDataToChild(                              &
                &  parentData  = meshData(lowElemIndex:upElemIndex), &
                &  nParentDofs = nDofs,                              &
                &  nChildDofs  = nChildDofs,                         &
                &  nComponents = nComponents,                        &
                &  nChilds     = nChilds,                            &
                &  projection  = projection,                         &
                &  childData   = childData                           )

              ! Iterate over all childDofs and set the data corectly in
              ! newMeshData
              lowChildIndex = upChildIndex + 1
              upChildIndex = (lowChildIndex-1) + nChilds * nChildDofs &
                &                                        * nComponents

              newMeshData(lowChildIndex:upChildIndex) = childData
              deallocate(childData)
            else
              ! Child cell won't be refined.
              !
              ! Need projecton from current dofs to 1 dof.
              allocate(childData(nComponents))
              ! Create lower and upper indices for all data of iElem in
              ! meshData.
              lowElemIndex = upElemIndex + 1
              upElemIndex = (lowElemIndex-1) + nDofs * nComponents

              ! Projection from nDofs to oneDof (integral mean valuea).
              oneDof = 1
              noChilds = 1
              call ply_projDataToChild(                              &
                &  parentData  = meshData(lowElemIndex:upElemIndex), &
                &  nParentDofs = nDofs,                              &
                &  nChildDofs  = oneDof,                             &
                &  nComponents = nComponents,                        &
                &  nChilds     = noChilds,                           &
                &  projection  = projection_oneDof,                  &
                &  childData   = childData                           )

              ! Iterate over all childDofs and set the data corectly in
              ! newMeshData
              lowChildIndex = upChildIndex + 1
              upChildIndex = (lowChildIndex-1) + nComponents

              newMeshData(lowChildIndex:upChildIndex) = childData
              deallocate(childData)
            end if
          end do childLoop

        else
          ! Parent cell wasn't refined so it contains data with only one dof
          ! left.
          ! Simple copying.
          allocate(childData(nComponents))

          lowElemIndex = upElemIndex + 1
          upElemIndex = (lowElemIndex-1) + nComponents

          childData = meshData(lowElemIndex:upElemIndex)

          lowChildIndex = upChildIndex + 1
          upChildIndex = (lowChildIndex-1) + nComponents

          newMeshData(lowChildIndex:upChildIndex) = childData
          childpos = childpos + 1
          deallocate(childData)

        end if

      end do elementLoop

    else

      elemLoop: do iElem=1,nElems
        if (new_refine_tree(iElem)) then
          allocate(childData(nChildDofs*nChilds*nComponents))
          ! Create lower and upper indices for all data of iElem in meshData.
          lowElemIndex = upElemIndex + 1
          upElemIndex = (lowElemIndex-1) + nDofs * nComponents

          ! Project these dofs from the coarse element to the
          ! finer elements.
          call ply_projDataToChild(                              &
            &  parentData  = meshData(lowElemIndex:upElemIndex), &
            &  nParentDofs = nDofs,                              &
            &  nChildDofs  = nChildDofs,                         &
            &  nComponents = nComponents,                        &
            &  nChilds     = nChilds,                            &
            &  projection  = projection,                         &
            &  childData   = childData                           )

          ! Iterate over all childDofs and set the data corectly in newMeshData
          lowChildIndex = upChildIndex + 1
          upChildIndex = (lowChildIndex-1) + nChilds * nChildDofs * nComponents

          newMeshData(lowChildIndex:upChildIndex) = childData
          deallocate(childData)
        else
          allocate(childData(nComponents))
          ! Create lower and upper indices for all data of iElem in meshData.
          lowElemIndex = upElemIndex + 1
          upElemIndex = (lowElemIndex-1) + nDofs * nComponents

          ! Projection from nDofs to oneDof (integral mean valuea).
          oneDof = 1
          noChilds = 1
          call ply_projDataToChild(                              &
            &  parentData  = meshData(lowElemIndex:upElemIndex), &
            &  nParentDofs = nDofs,                              &
            &  nChildDofs  = oneDof,                             &
            &  nComponents = nComponents,                        &
            &  nChilds     = noChilds,                           &
            &  projection  = projection_oneDof,                  &
            &  childData   = childData                           )

          ! Iterate over all childDofs and set the data corectly in newMeshData
          lowChildIndex = upChildIndex + 1
          upChildIndex = (lowChildIndex-1) + nComponents

          newMeshData(lowChildIndex:upChildIndex) = childData
          deallocate(childData)
        end if
      end do elemLoop
    end if

  end subroutine ply_subsampleData
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to project elemental data from a parent cell to one of
  !! its children.
  subroutine ply_projDataToChild( parentData, nParentDofs, nChildDofs,        &
    &                             nComponents, nChilds, projection, childData )
    ! -------------------------------------------------------------------- !
    !> Linearized data for a single variable (can have multiple components)
    !! and a single degree of freedom of the parent cell.
    real(kind=rk), intent(in) :: parentData(:)

    !> The number of dofs of the parent element.
    integer, intent(in) :: nParentDofs

    !> The total number of dofs for the child cells.
    integer, intent(in) :: nChildDofs

    !> The number of componentns of the given variable.
    integer, intent(in) :: nComponents

    !> The number of children.
    integer, intent(in) :: nChilds

    !> The information about the projection coefficients for the parent
    !! dofs to the child dofs.
    type(ply_ProjCoeff_type), intent(in) :: projection

    !> The created childData.
    real(kind=rk), intent(out) :: childData(:)
    ! -------------------------------------------------------------------- !
    integer :: iChildDof, iComp, iChild, iParentDof
    integer :: childDof_pos, parentDof_pos
    real(kind=rk) :: projCoeff
    ! -------------------------------------------------------------------- !
    childData(:) = 0.0_rk

    childLoop: do iChild = 1, nChilds
      parentDofLoop: do iParentDof = 1, nParentDofs
        childDofLoop: do iChildDof = 1, nChildDofs
          ! Get the projection coefficient for iChild, parentDof and iChildDof.
          projCoeff = projection%projCoeff( iParentDof, iChildDof, iChild )
          compLoop: do iComp = 1, nComponents
            ! Position of the child dof in childData.
            childDof_pos = iComp + (iChildDof - 1) * nComponents &
              &          + (iChild - 1) * nChildDofs * nComponents
            ! Position of the parent dof in parentData.
            parentDof_pos = iComp + (iParentDof - 1) * nComponents

            childData( childDof_pos ) = childData( childDof_pos ) &
              &                       + projCoeff * parentData( parentDof_pos )

          end do compLoop
        end do childDofLoop
      end do parentDofLoop
    end do childLoop

  end subroutine ply_projDataToChild
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> subroutine to create gauss points and weights for one-dimensional
  !! integration on the interval [x1,x2].
  subroutine ply_gauleg( x1, x2, x, w, nIntP )
    ! -------------------------------------------------------------------- !
    !> The coordinates of the gauss points on the interval [-1,1].
    !! The array has the length nIntP.
    real(kind=rk), allocatable, intent(inout) :: x(:)

    !> The quadrature weights. The array has the length nIntP.
    real(kind=rk), allocatable, intent(inout) :: w(:)

    !> The number of integration points.
    integer, intent(in) :: nIntP

    !> lower limit of integration interval
    real(kind=rk), intent(in) :: x1

    !> upper limit of integration interval
    real(kind=rk), intent(in) :: x2
    ! -------------------------------------------------------------------- !
    ! some work variables
    real(kind=rk) :: z1,z,xm,xl,pp,p3,p2,p1;
    ! the relative precision of the points
    real(kind=rk) :: EPS
    integer :: m, i, j
    ! -------------------------------------------------------------------- !

    allocate(x(nIntP), w(nIntP))

    EPS= 1.0 / (10.0**(PRECISION(1.0_rk)-2) )
    m = (nIntP+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);

    do i = 1, m

      z=cos(PI*((i-1)+0.75_rk)/(nIntP+0.5_rk));

      loopToExit : do
        p1=1.0_rk; p2=0.0_rk;
        do j=0 , nIntP-1
          p3=p2;
          p2=p1;
          p1=((2.0_rk*j+1.0_rk)*z*p2-j*p3)/(j+1.0_rk);
        end do
        pp=nIntP*(z*p1-p2)/(z*z-1.0_rk);
        z1=z;
        z=z1-p1/pp;
        if ( abs(z-z1) < EPS ) then
          exit loopToExit
        end if
      end do loopToExit

      x(i)=xm-xl*z;
      x(nIntP-i+1)=xm+xl*z;
      w(i)=2.0_rk*xl/((1.0_rk-z*z)*pp*pp);
      w(nIntp-i+1)=w(i);

    end do

  end subroutine ply_gauleg
  ! ************************************************************************ !

end module ply_LegPolyProjection_module
! **************************************************************************** !
