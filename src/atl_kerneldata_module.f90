! Copyright (c) 2011 Metin Cakircali <m.cakircali@grs-sim.de>
! Copyright (c) 2011-2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2014, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Laura Didinger <l.didinger@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2019 Robin Weihe <robin.weihe@student.uni-siegen.de>
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

! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
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
!
! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of
!   modal coefficients in.


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
!> author: Jens Zudrop
!! module containing all the informations and subroutines about the state
!! of a solver.
module atl_kerneldata_module
  use env_module,      only: rk
  use tem_time_module, only: tem_time_type
  use ply_dof_module,  only: P_Space, Q_space

  implicit none

  private


  !> summary: Data type to describe the state of an equation in a solver.
  !!
  !! This data type shares some properties with cube_container_type type like
  !! the number of elements on this refinement level of a cell type.
  type atl_statedata_type
    !> Local time of this kernel
    type(tem_time_type) :: local_time
    !> Local output time of this kernel
    real(kind=rk) :: local_output_time

    !> State collects the states of the variables of our simulation.
    !! So the dimension of this array is as follows:
    !! The first dimension is the number of elements. If we assume that
    !! our atl_kerneldata_type variable is located inside
    !! cube_container_type%state_list(index) then the number of
    !! elements is given by
    !! cube_container_type%mesh_list(index)%descriptor%nElems_fluid
    !! + cube_container_type%mesh_list(index)%descriptor%nElems_ghostFromCoarser
    !! + cube_container_type%mesh_list(index)%descriptor%nElems_ghostFromFiner
    !! + cube_container_type%mesh_list(index)%descriptor%nElems_halo .
    !! The second dimension is the number of degree of freedoms. This
    !! information is stored in side atl_kerneldata_type%scheme%nDoFs
    !! and no where else, since p-adaptivity might change the number of
    !! degrees of freedoms for other parts of the mesh.
    !! The third dimension is the number of variables we have to simulate
    !! for our equation(s). It is stored inside Equations_type%nScalars.
    real(kind=rk), allocatable :: state(:,:,:)

  end type atl_statedata_type


  !> summary: Data type to describe the kernelstate of a specific kernel.
  !!
  !! This data type shares some properties with cube_container_type type like
  !! the number of elements on this refinement level of a cell type.
  type atl_kerneldata_type

    !> The total number of cells (including only fluid, ghost, halo and
    !! boundary cells).
    integer :: nTotal

    !> The number of scalar variables of the current equation.
    integer :: nVars

    !> Maximal polynomial degree of the data in this kerneldata.
    integer :: maxPolyDegree

    !> Chosen tensor kind of the polynomial representation in this
    !! kerneldata.
    integer :: poly_space

    !> Number of dimensions of the polynomial in this kerneldata.
    integer :: nDims

    !> The number of degrees of freedom per scalar variable of your equation.
    integer :: nDofs

    !> The number of derived quantities the kernel will use in the future.
    integer :: nDervQuant

    !> array of derived states. Could be anything like derivatives
    !! face values, etc. The only thing that is important is that
    !! the kernel has to handle the data consistently.
    !! The first dimension is the number of elemnts (including fluid,
    !! ghost, halo and boundary cells).
    !! The second dimension is the number of derived quantities (e.g.
    !! the fave value and derivatives at the faces). The exact meaning
    !! of this dimension is specified by the kernel.
    !! The third dimension is the number of informations per derived
    !! quantity per cell (e.g. the number of faces times the quadrature
    !! points). The exact meaning of this dimension is specified by
    !! the kernel.
    !! The fourth dimension is the number of varibales of the equation.
    real(kind=rk), allocatable :: state_der(:,:,:)

    !> Flag to indicate, whether to compute
    !! the deviation.
    logical :: need_deviation = .false.

    !> Flag to indicate, whether to compute maximal estimates for derivatives.
    logical :: need_maxgrad = .false.

    !> Maximal deviation bound of the polynomial in state.
    real(kind=rk), allocatable :: deviation(:,:)

    !> Limit for maximal size of derivative in state.
    real(kind=rk), allocatable :: maxgrad(:,:)
  end type atl_kerneldata_type

  public :: atl_statedata_type
  public :: atl_init_statedata
  public :: atl_kerneldata_type, atl_init_kerneldata
  public :: atl_kerneldata_update_maxdev
  public :: atl_kerneldata_update_maxgrad
  public :: atl_kerneldata_update_estimates


contains


  !> Initialize the statedata.
  subroutine atl_init_statedata(statedata, nTotal, nDofs, nVars, time)
    type(atl_statedata_type), intent(inout) :: statedata

    !> The total number of cells (including only fluid, ghost, halo and
    !! boundary cells).
    integer, intent(in) :: nTotal

    !> The number of scalar variables of the current equation.
    integer, intent(in) :: nVars

    !> The number of degrees of freedom per scalar variable of your equation.
    integer, intent(in) :: nDofs

    !> current time
    type(tem_time_type), intent(in) :: time
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    !now allocate enough memory for variable's array
    !...and set it to default values
    allocate(statedata%state( nTotal, &
      &                       nDoFs,  &
      &                       nVars)  )
    statedata%local_time = time
    statedata%state = 0.0_rk

  end subroutine atl_init_statedata


  !> summary: init routine for the kerneldata type.
  !!
  !! Subroutine to init the atl_kerneldata_type. This routine should be
  !! used in the init routine of the kernel, since the dimensions of the
  !! fields in the atl_kerneldata_type depend on the combination of scheme
  !! and equation.
  subroutine atl_init_kerneldata( kerneldata, statedata, nTotal, nVars, nDofs, &
    &                             nDervQuant, time,                            &
    &                             maxPolyDegree, nDims, poly_space,            &
    &                             need_deviation, need_maxgrad                 )
    ! -------------------------------------------------------------------------
    !> The data type to initialize.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The data type to initialize.

    type(atl_statedata_type), intent(inout) :: statedata

    !> The total number of cells (including only fluid, ghost, halo and
    !! boundary cells).
    integer, intent(in) :: nTotal

    !> The number of scalar variables of the current equation.
    integer, intent(in) :: nVars

    !> The number of degrees of freedom per scalar variable of your equation.
    integer, intent(in) :: nDofs

    !> The number of derived quantities the kernel will use in the future.
    integer, intent(in) :: nDervQuant

    !> current time
    type(tem_time_type), intent(in) :: time

    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree

    !> Number of dimensions for the polynomials in the kerneldata.
    integer, intent(in) :: nDims

    !> Kind of tensor product used to represent multidimensional polynomials
    !! in this kerneldata.
    integer, intent(in) :: poly_space

    !> Should the maximal deviation be computed for the state?
    logical, intent(in) :: need_deviation

    !> Should the maximal gradient be computed for the state?
    logical, intent(in) :: need_maxgrad
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    call atl_init_statedata(statedata = statedata, &
      &                     nTotal    = nTotal,    &
      &                     nDofs     = nDofs,     &
      &                     nVars     = nVars,     &
      &                     time      = time       )

    ! the derived data
    kerneldata%need_deviation = need_deviation
    kerneldata%need_maxgrad   = need_maxgrad
    allocate(kerneldata%deviation(nTotal, nVars))
    kerneldata%deviation = 0.0_rk
    allocate(kerneldata%maxgrad(nTotal, nVars))
    kerneldata%maxgrad = 0.0_rk

    allocate(kerneldata%state_der( nTotal, &
      &                            nDervQuant+modulo(maxPolyDegree,2), &
      &                            nVars ) )
    kerneldata%nTotal           = nTotal
    kerneldata%nVars            = nVars
    kerneldata%nDofs            = nDofs
    kerneldata%nDims            = nDims
    kerneldata%poly_space       = poly_space
    kerneldata%maxPolyDegree    = maxPolyDegree
    kerneldata%nDervQuant       = nDervQuant
    kerneldata%state_der(:,:,:) = 0.0_rk

  end subroutine atl_init_kerneldata


  !> Find the maximal deviation for the polynomials representing the state
  !! in each element.
  !!
  !! This additional information is useful for estimates during the
  !! simulation.
  subroutine atl_kerneldata_update_maxdev(statedata, kerneldata)
    !> The statedata to update the bounds for.
    !!
    !! deviation of kerneldata will be recomputed according to
    !! the statedata.
    type(atl_statedata_type), intent(in) :: statedata

    !> The kerneldata storing the deviation, this will be updated
    !! according to the data provided in statedata.
    type(atl_kerneldata_type), intent(inout) :: kerneldata

    integer :: nElems, iElem
    integer :: nVars, iVar

    if (kerneldata%need_deviation) then
      nVars  = size(statedata%state, 3)
      nElems = size(statedata%state, 1)

      do iVar=1,nVars
        do iElem=1,nElems
          kerneldata%deviation(iElem,iVar) &
            &  = sum(abs(statedata%state(iElem,2:,iVar)))
        end do
      end do
    end if

  end subroutine atl_kerneldata_update_maxdev


  !> Find the maximal gradient estimation for the polynomials representing
  !! the state in each element.
  !!
  !! This additional information is useful for estimates during the
  !! simulation.
  subroutine atl_kerneldata_update_maxgrad(statedata, kerneldata)
    !> The statedata to update the bounds for.
    !!
    !! maxgrad of kerneldata will be recomputed according to
    !! the statedata.
    type(atl_statedata_type), intent(in) :: statedata

    !> The kerneldata storing the gradients, this will be updated
    !! according to the data provided in statedata.
    type(atl_kerneldata_type), intent(inout) :: kerneldata

    real(kind=rk) :: maxgrad
    integer :: nElems, iElem
    integer :: nVars, iVar, iDof
    integer :: i,j,k
    integer :: n
    integer :: nDofs

    if (kerneldata%need_maxgrad) then

      nVars  = size(statedata%state, 3)
      nDofs  = kerneldata%nDofs
      nElems = size(statedata%state, 1)
      n = kerneldata%maxPolyDegree + 1

      select case(kerneldata%nDims)
      case (1)
        do iVar=1,nVars
          do iElem=1,nElems
            maxgrad = 0.0_rk
            do iDof=2,nDofs
              maxgrad = maxgrad &
                &     + statedata%state(iElem,iDof,iVar)**2 &
                &       * 0.25_rk * (iDof*(iDof-1))**2
            end do
          end do
        end do

      case (2)
        if (kerneldata%poly_space /= P_space) then
          ! If the data is not in a P_space Tensor, use Q-Tensor representation
          do iVar=1,nVars
            do iElem=1,nElems
              maxgrad = 0.0_rk
              do iDof=2,nDofs
                j = iDof / n
                i = iDof - j*n
                maxgrad = maxgrad &
                  &     + statedata%state(iElem,iDof,iVar)**2 &
                  &       * 0.25_rk * (  (i*(i+1))**2 &
                  &                    + (j*(j+1))**2 )
              end do
              kerneldata%maxgrad(iElem, iVar) = maxgrad
            end do
          end do
        else
          ! P Polynomial representation
          do iVar=1,nVars
            do iElem=1,nElems
              maxgrad = 0.0_rk
              i = 1
              j = 1
              do iDof=2,nDofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.
  if (i .ne. 1) then
    ! next item
    i = i - 1
    j = j + 1
  else
    ! next layer
    i = j + 1
    j = 1
  end if
                maxgrad = maxgrad &
                  &     + statedata%state(iElem,iDof,iVar)**2 &
                  &       * 0.25_rk * (  (i*(i-1))**2 &
                  &                    + (j*(j-1))**2 )
              end do
              kerneldata%maxgrad(iElem, iVar) = maxgrad
            end do
          end do
        end if

      case (3)
        if (kerneldata%poly_space /= P_space) then
          ! If the data is not in a P_space Tensor, use Q-Tensor representation
          do iVar=1,nVars
            do iElem=1,nElems
              maxgrad = 0.0_rk
              do iDof=2,nDofs
                k = iDof / n**2
                j = (iDof - k*n**2) / n
                i = iDof - k*n**2 - j*n
                maxgrad = maxgrad &
                  &     + statedata%state(iElem,iDof,iVar)**2 &
                  &       * 0.25_rk * (  (i*(i+1))**2 &
                  &                    + (j*(j+1))**2 &
                  &                    + (k*(k+1))**2 )
              end do
              kerneldata%maxgrad(iElem, iVar) = maxgrad
            end do
          end do
        else
          ! P Polynomial representation
          do iVar=1,nVars
            do iElem=1,nElems
              maxgrad = 0.0_rk
              i = 1
              j = 1
              k = 1
              do iDof=2,nDofs
  ! - ansatz indices are arranged in layers. within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. within a block, ansfuncz is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. each item represents one particular
  !   combination of ansfuncx, -y, and -z degrees.

  if (i .ne. 1) then
    ! next item
    i = i - 1
    j = j + 1
  elseif (j .ne. 1) then
    ! next block
    i = j - 1
    j = 1
    k = k + 1
  else
    ! next layer
    i = k + 1
    j = 1
    k = 1
  end if
                maxgrad = maxgrad &
                  &     + statedata%state(iElem,iDof,iVar)**2 &
                  &       * 0.25_rk * (  (i*(i-1))**2 &
                  &                    + (j*(j-1))**2 &
                  &                    + (k*(k-1))**2 )
              end do
              kerneldata%maxgrad(iElem, iVar) = maxgrad
            end do
          end do
        end if
      end select

    end if

  end subroutine atl_kerneldata_update_maxgrad



  subroutine atl_kerneldata_update_estimates(statedata, kerneldata)
    !> The statedata to update the estimates for.
    !!
    !! maxgrad of kerneldata will be recomputed according to
    !! the statedata.
    type(atl_statedata_type), intent(in) :: statedata

    !> The kerneldata storing the estimates, this will be updated
    !! according to the data provided in statedata.
    type(atl_kerneldata_type), intent(inout) :: kerneldata

    call atl_kerneldata_update_maxdev( &
      &    statedata  = statedata,     &
      &    kerneldata = kerneldata     )

    call atl_kerneldata_update_maxgrad( &
      &    statedata  = statedata,      &
      &    kerneldata = kerneldata      )

  end subroutine atl_kerneldata_update_estimates

end module atl_kerneldata_module
