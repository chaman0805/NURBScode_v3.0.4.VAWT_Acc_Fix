!------------------------------------------------------------------------
!    Module for storing arrays and allocation routines
!------------------------------------------------------------------------     
module types_beam

  implicit none

  ! Declare type mesh
  type :: mesh_beam
    ! Order in U for each patch
    integer, allocatable :: P(:)

    ! patch type
    ! 0-cable (only tension or compression); 1-rod..
    integer, allocatable :: PTYPE(:)

    ! Knot vectors in u for each element
    real(8), allocatable :: U_KNOT(:,:)

    ! Size of the knot vector for each element (e.g. NUK=P+MCP+1)
    integer, allocatable :: NUK(:)

    ! Control Net
    real(8), allocatable :: B_NET_D_alphaf(:,:), B_NET_Dt_alphaf(:,:), B_NET_DDt_alpham(:,:), &
                            B_NET(:,:), B_NET_U(:,:), B_NET_D(:,:), B_NET_Dt(:,:), B_NET_DDt(:,:), &
                            B_NET_D_old(:,:), B_NET_Dt_old(:,:), B_NET_DDt_old(:,:)

    ! Boundary condition indicator for global nodes and edges
    ! respectively
    integer, allocatable :: IBC(:,:)

    ! array to store force vectors on the wind turbine blade
    real(8), allocatable :: FORCE(:,:)

    ! Global connectivity array
    integer, allocatable :: IEN(:,:), INN(:)

    ! number of shape functions for every element
    integer, allocatable :: NSHL(:)

    ! Bezier extraction operator
    real(8), allocatable :: Ext(:,:,:)

    ! Array for closest points
    integer, allocatable :: CLE(:,:)
    real(8), allocatable :: CLP(:,:,:)

    integer :: NGAUSS, NNODE, NNODE_LOC, NEL, maxNSHL

  end type mesh_beam


  ! Declare type mesh for multi-patch
  type :: mesh_mp_beam
    ! degree in U for each patch
    integer, allocatable :: P(:)

    ! number of control points in U for each patch
    ! (no need for T-spline)
    integer, allocatable :: MCP(:)

    ! Total number of control points and elements for each patch
    integer, allocatable :: NNODE(:), NEL(:)

    ! patch type
    ! 0-cable (only tension or compression); 1-rod...
    integer, allocatable :: PTYPE(:)

    ! Knot vectors in u for each patch
    real(8), allocatable :: U_KNOT(:,:)

    ! Control Net
    real(8), allocatable :: B_NET(:,:,:)

    ! Boundary condition indicator for global nodes and edges
    ! respectively
    integer, allocatable :: IBC(:,:,:)

    ! array to store force vectors on the wind turbine blade
    real(8), allocatable :: FORCE(:,:,:)

    ! Mapping between patches and global reduced numbering
    ! e.g., MAP(1,2) = 3 means 1st patch, 2nd node points to 
    !   global reduced node number 3
    integer, allocatable :: MAP(:,:)
  end type mesh_mp_beam

  ! Declare type beam
  type :: beam
    ! number of patches for cable (C)
    ! number of patches for Rod (R)
    ! number of patches
    integer :: NPC, NPR, NPT

    ! material matrix for beam
    integer :: NMat
    real(8), allocatable :: Ec(:), Dens(:), Rad(:), C_dp(:)
    real(8) :: f_fact, G_fact(3)
  end type beam

end module types_beam
