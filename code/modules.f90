!------------------------------------------------------------------------
! 
!------------------------------------------------------------------------   
module class_def
  
  type NURBSpatch
  
    integer :: P, Q, R
    integer :: MCP, NCP, OCP
    real(8), allocatable :: U_KNOT(:), V_KNOT(:), W_KNOT(:)
    
  end type NURBSpatch

  type bnd_class

    integer :: FACE_ID

    integer :: NFACE  
  
    integer, allocatable :: FACE_IEN(:,:)
    
    integer, allocatable :: F2E(:)
    integer, allocatable :: FACE_OR(:)
  
    integer :: NNODE
    integer, allocatable :: BNODES(:)
  
    ! mapping between partitioned (local) boundary node/element
    ! and unpartitioned boundary (shell) node/element
    integer, allocatable :: L2SNODE(:), L2SELEM(:)

  end type bnd_class
  
end module class_def



!------------------------------------------------------------------------
! Module for defining shell types and variables
!------------------------------------------------------------------------     
module defs_shell

  implicit none

  ! Declare type mesh
  type :: mesh
    ! degree in U and V for each patch
    integer, allocatable :: P(:), Q(:)

    ! patch type
    ! 1-blade; 0-bending strips; 2-shear web; ...
    integer, allocatable :: PTYPE(:)

    ! Knot vectors in u and v directions for each element
    real(8), allocatable :: U_KNOT(:,:), V_KNOT(:,:)

    ! Size of the knot vector for each elements (e.g. NUK=P+MCP+1)
    integer, allocatable :: NUK(:), NVK(:)

    ! Control Net
    ! B_NET is reference config, B_NET_D is current config.
    ! For the pre-bend case, reference config. could change, so 
    ! we created B_NET_U for undeformed, original config.
    real(8), allocatable :: B_NET(:,:), B_NET_U(:,:), B_NET_D(:,:)

    ! Boundary condition indicator for global nodes and edges
    ! respectively
    integer, allocatable :: IBC(:,:)

    ! array to store force vectors on the wind turbine blade
    real(8), allocatable :: FORCE(:,:)

    ! Global connectivity array
    integer, allocatable :: IEN(:,:), INN(:,:)

    ! number of shape functions for every element
    integer, allocatable :: NSHL(:)

    ! Bezier extraction operator
    real(8), allocatable :: Ext(:,:,:)

    ! Array for closest points (and the corresponding element)
    integer, allocatable :: CLE(:,:)
    real(8), allocatable :: CLP(:,:,:)

    ! Array for element list
    integer, allocatable :: Elm_Close(:,:,:)
    integer :: NEL_Close, NEL_Close_t

    ! Array for element gauss point location and for buiding the 
    ! element list based on radial location
    real(8) :: Elm_Size
    real(8), allocatable :: Elm_Loc(:,:,:)
    integer, allocatable :: RAD_ELM_LIST(:,:,:), RAD_ELM_NUM(:,:)

    integer :: NGAUSS = 0, NNODE = 0, NEL = 0, maxNSHL = 0, NNODE_LOC = 0

    ! Store node number of the tip
    integer :: TipLoc, TipLocTr

    ! array for Solution vectors
    real(8), allocatable :: dsh(:,:), dshold(:,:), &
                            ush(:,:), ushold(:,:), &
                            ash(:,:), ashold(:,:)

    ! Surface ID and Boundary number
    integer :: FaceID, iBound
  end type mesh


  ! Declare type mesh for multi-patch
  type :: mesh_mp
    ! degree in U and V for each patch
    integer, allocatable :: P(:), Q(:)

    ! number of control points in U and V for each patch
    ! (no need for T-spline)
    integer, allocatable :: MCP(:), NCP(:)

    ! Total number of control points and elements for each patch
    integer, allocatable :: NNODE(:), NEL(:)

    ! patch type
    ! 1-blade surface; 0-bending strips; 2-shear web; ...
    integer, allocatable :: PTYPE(:)

    ! Knot vectors in u and v directions for each patch
    real(8), allocatable :: U_KNOT(:,:), V_KNOT(:,:)

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
  end type mesh_mp


  ! Declare type shell (for wind turbine blade)
  type :: shell_bld

    type(mesh_mp) :: mNRB
    type(mesh)    :: NRB

    type(mesh)    :: TSP, BEZ

    type(mesh)    :: FEM

    ! number of patches for Blade Surface (S). Matches the fluid mesh
    ! number of patches for blade structure (B). May include shear webs
    ! number of total patches including bending strips (T)
    integer :: NPS, NPB, NPT

    integer :: M_Flag, T_Flag
    real(8) :: RHSGtol, G_fact(3), RHSGNorm

    ! row, col, and total of nonzero entries for sparse structure
    integer, allocatable :: row(:), col(:)
    integer :: icnt

    ! The right hand side load vector G and left hand stiffness matrix K
    real(8), allocatable :: RHSG(:,:), LHSK(:,:), &
                            RHSG_EXT(:,:), RHSG_GRA(:,:)

!    ! Solution vectors
!    real(8), allocatable :: yg(:,:), dg(:,:), tg(:), &
!                            mg(:), dl(:,:)

    ! material matrix for composite
    integer :: NMat
!    real(8), allocatable :: matA(:,:,:), matB(:,:,:), matD(:,:,:)
    real(8), allocatable :: matA(:,:,:,:), matB(:,:,:,:), matD(:,:,:,:), Density(:,:), Thickness(:,:)

    ! number of newton iterations for shell
    integer, allocatable :: Nnewt(:)

    ! Torque computed on shell mesh
    real(8) :: Tq1, Tq2

    ! Blade rotation. 0 degree is the straight-up position
    real(8) :: BldRot

    integer :: bmap
  end type shell_bld


  ! Declare type shell (for non-matching boundaries)
  type :: shell_nmb

    type(mesh), allocatable :: FEM(:)

  end type shell_nmb
end module defs_shell



!------------------------------------------------------------------------
!     Module for storing arrays and allocation routines       
!------------------------------------------------------------------------  
module aAdjKeep
  
  use class_def
  use defs_shell

  implicit none
  save
   
  ! Mesh
  real(8), allocatable :: xg(:,:), wg(:)

  integer, allocatable :: IEN(:,:), EPID(:), EIJK(:,:), NodeID(:)
  
  type(bnd_class),  allocatable :: bound(:) 
  type(NURBSpatch), allocatable :: patch(:)
    
  ! Contraint flags
  integer, allocatable :: IPER(:)
  integer, allocatable :: IBC(:,:)
  
  ! Type flags
  integer, allocatable :: EL_TYP(:), D_FLAG(:), P_FLAG(:)
  
  ! Spars Struc
  integer, allocatable :: row(:), col(:)

  ! The right hand side load vector G and left hand stiffness matrix K
  real(8), allocatable :: rhsGu(:,:), rhsGm(:,:), rhsGp(:), rhsGls(:)
  real(8), allocatable :: lhsK11(:,:), lhsK12(:,:), lhsK22(:,:), &
                          lhsG(:,:),   lhsD1(:,:),  lhsD2(:,:),  &
                          lhsM(:),     lhsLS(:),    LHSPi(:,:),  &
                          LHSPti(:,:), LHSKi(:,:),  lhsMi(:),    &
                          LHSLSi(:),   LHSmass(:),               &
                          LHSlsu(:,:), LHSPls(:), LHSUls(:,:)

  ! Solution vectors
  real(8), allocatable :: dg(:,:),   dgold(:,:),   &
                          ug(:,:),   ugold(:,:),   &
                          ugm(:,:),  ugmold(:,:),  &
                          acg(:,:),  acgold(:,:),  &
                          acgm(:,:), acgmold(:,:), &
                          pg(:),     pgold(:),     &
                          phig(:),   phigold(:),   &
                          rphig(:),  rphigold(:)

  real(8), allocatable :: uavg(:,:), pavg(:)

  ! Rigid body
  real(8) :: vbn0(3),  vbn1(3) 
  real(8) :: dbn0(3),  dbn1(3) 
  real(8) :: wbn0(3),  wbn1(3) 
  real(8) :: Rn0(3,3), Rn1(3,3) 
      
  ! First P-K Stress
  real(8), allocatable :: FPKS(:,:,:)

  ! Array for Prism
  integer, allocatable :: ELMNSHL(:), ELMNGAUSS(:)

  ! global information for individual blades
  type(mesh) :: blade(3)

end module aAdjKeep

!----------------------------------------------------------------------
! Module for common variables  
!---------------------------------------------------------------------- 
module commonvars

  implicit none
  save

  ! Mesh 
  integer :: NSD, NNODE, NELEM, NBOUND, NPATCH, NSHLB, icnt, &
             NGAUSSTET, NGAUSSPRI, NGAUSSB, NBlade, maxNSHL
  logical :: iga 
  real(8) :: DetJ, DetJb, DetJinv, hglob

  ! Time step 
  real(8) :: Delt, Dtgl, rhoinf, beti, gami, alfi, almi, &
             mgam, ogam, lambda, time, &
             conv_time, mono_time, move_time, shel_time
  integer :: Nstep, ifq, ifq_sh, ifq_tq, mono_iter

  ! Navier-Stokes solver
  real(8) :: mua, rhoa, muw, rhow
  real(8) :: NS_kdc_w, NS_kdc_a
    
  real(8) :: NS_GMRES_tol, NS_NL_Utol, NS_NL_Ptol
  integer :: NS_GMRES_itermin,NS_GMRES_itermax, NS_NL_itermax, &
             NS_hess_flag
  
  ! LevelSet Convection solver  
  real(8) :: LSC_kdc
  real(8) :: LSC_GMRES_tol, LSC_NL_tol
  integer :: LSC_GMRES_itermin,LSC_GMRES_itermax, LSC_NL_itermax, &
             LSC_pred_step

  ! Rigid body solver
  real(8) :: RB_NL_Ftol,RB_NL_Mtol
  
  ! Mesh solver  
  real(8) :: Mesh_GMRES_tol,Mesh_NL_tol
  integer :: Mesh_GMRES_itermin,Mesh_GMRES_itermax

  ! Level Set Redistance solver  
  real(8) :: LSRD_kdc, LSRD_penalty_fac
  real(8) :: LSRD_pseudoDTGL_fac  
  real(8) :: LSRD_GMRES_tol, LSRD_NL_tol
  integer :: LSRD_GMRES_itermin,LSRD_GMRES_itermax, LSRD_NL_itermax
  
  ! Level Set massfix     
  real(8) :: Mass_init
  real(8) :: Mass_NL_tol
  integer :: Mass_NL_itermax    
  
  ! Rigid body         
  real(8) :: massb, Ibhat(3,3), xcg(3), Fb(3),Mb(3)
  integer :: VBC(3), MBC(3)
  
  ! Domain and Hull
  real(8) :: domain_top,  domain_bottom, &
             domain_left, domain_right, &
             domain_front,domain_back, &
             water_level,water_depth,air_height,hull_length 
    
  ! Wave generating wall    
  real(8) :: wave_periode, wave_length, wave_amp, wave_angle 
     
  ! Gravity     
  real(8) :: gravity, gravvec(3)

  ! Interface
  real(8) :: mp_eps, mu, rho, dmudphi, drhodphi 
    
  ! Setup   
  real(8) :: Froude,Uin   
  integer :: BCtype(99)  
  
  integer, allocatable :: BCugType(:,:)
  real(8), allocatable :: BCugValu(:,:)

  ! Flags
  logical :: move, mono, conv, shel
    
  ! rotation
  real(8) :: theta, thetd, thedd, thetaOld, thetdOld, theddOld, &
             maxthetd, torque1, torque2, torque3, torque4, torque_SH, &
             moment1, moment2

  logical :: solshell, nonmatch

  ! Rotation

!  real(8) :: x0_c(3), xt_c(3) !x0_1(3), x0_2(3), x0_3(3), xt_1(3), xt_2(3), xt_3(3), 
!A  real(8) :: e0(3,3), et(3,3), Rt(3,3), Identity(3,3) temp_vec1(3), temp_vec2(3), temp1, temp2

  real(8) :: Rmat(3,3), RmatOld(3,3), Rtang(3,3), RtangOld(3,3), xrot(3), xrotOld(3), xh_E(3), &
             Omega(3), OmegaOld(3), Omegatang(3), OmegatangOld(3), normal(3), Identity(3,3), &
             Rfix(3,3), norm(2,3)

end module commonvars



!----------------------------------------------------------------------
! Module for common parameters (e.g. pi)
!---------------------------------------------------------------------- 
module commonpars
  implicit none
  save

  real(8), parameter :: pi = 3.14159265358979323846264338328d0
end module commonpars



!----------------------------------------------------------------------
!     Module for mpi
!----------------------------------------------------------------------
module mpi
  implicit none
  save

  include "mpif.h"

  integer, parameter :: mpi_master = 0
  integer, parameter :: maxtask    = 50
  integer, parameter :: maxseg     = 15000
      
  integer :: numnodes, myid, mpi_err
  integer :: status(MPI_STATUS_SIZE)
  integer :: lstrideu, lstridep, numtask
  integer :: lfrontu, maxfrontu
  integer :: lfrontp, maxfrontp
  integer :: nlworku, nlworkp
  integer :: itag, iacc, iother, numseg, isgbeg, itkbeg, isgend
  integer :: sevsegtypeu(maxtask,15)
  integer :: sevsegtypep(maxtask,15)
  logical :: ismaster
            
  integer, allocatable :: ilworku(:), ilworkp(:)
end module mpi



!------------------------------------------------------------------------
!    Module for beam
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


!------------------------------------------------------------------------
!    Module for structure types
!------------------------------------------------------------------------     
module types_structure

  implicit none

  ! Declare type shell
  type :: structure

    integer :: S_Flag, Nnewt, Nincr, Rincr, Nstep, C_Flag, M_Flag, &
               T_Flag, Sol_Flag, FF_Flag, Nincr_step
    real(8) :: RHSGtol, DGtol, RHSGNorm, DGNorm, &
               RHSFact, Delta_t

    ! row, col, and total of nonzero entries for sparse structure
    integer, allocatable :: row(:), col(:)
    integer :: icnt

    integer, allocatable :: ndiag(:)

    ! The right hand side load vector G and left hand stiffness matrix K
    real(8), allocatable :: RHSG(:,:), LHSK(:,:), LHSKV(:), &
                            RHSG_EXT(:,:), RHSG_GRA(:,:)

    ! Solution vectors
    real(8), allocatable :: yg(:,:), dg(:,:), tg(:), &
                            mg(:), dl(:,:), ygv(:), ddu(:,:)

  end type structure

end module types_structure
