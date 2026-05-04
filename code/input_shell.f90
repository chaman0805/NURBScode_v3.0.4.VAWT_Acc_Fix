!=======================================================================
! Subroutine to call all the necessary shell input subroutines
!=======================================================================
subroutine input_shell_nmb(NSD, NM, Center_Rot)
  use aAdjKeep
  use mpi
  use defs_shell
  implicit none

  type(shell_nmb), intent(out) :: NM
  integer,         intent(in)  :: NSD

  real(8), intent(in) :: Center_Rot(3)

  integer :: so, sf, i, ni, nj, nk, ii, j, k, ip, ier, &
             sumNDZ, sumNEL, &
             loop2, loop3, lnod, eloc, eglob, &
             p, q, mcp, ncp, nshl, nel, nnode, &
             maxP, maxQ, maxMCP, maxNCP, maxNNODE, &
             NEL_CLOSE

  character(len=30) :: fname, cname

  allocate(NM%FEM(2))

  !*************************************
  NM%FEM(1)%iBound = 4    ! inner
  NM%FEM(2)%iBound = 5    ! outer
  !*************************************

  NM%FEM(1)%FaceID = bound(NM%FEM(1)%iBound)%Face_ID  ! 11
  NM%FEM(2)%FaceID = bound(NM%FEM(2)%iBound)%Face_ID  ! 12

  ! get the fem shell mesh
  call shell_input_fem(NSD, NM%FEM(1))
  call shell_input_fem(NSD, NM%FEM(2))

  NM%FEM(:)%NGAUSS = 1    ! in 2D (total)

  !--------------------------------
  ! Find the elememt size
  !--------------------------------
  if (ismaster) write(*,*) "NM: Compute element size" 
  call fem_find_elm_size(NM%FEM(1), NSD)
  call fem_find_elm_size(NM%FEM(2), NSD)

  !------------------------------------------------------
  ! Compute the elememt gauss points physical location
  !------------------------------------------------------
  if (ismaster) write(*,*) "NM: Compute the elememt gauss points physical location" 
  do i = 1, 2
    allocate(NM%FEM(i)%Elm_Loc(NM%FEM(i)%NEL,NM%FEM(i)%NGAUSS,4))
    NM%FEM(i)%Elm_Loc = 0.0d0

    call fem_find_elm_loc(NM%FEM(i), NSD)
  end do

  !--------------------------------------------------
  ! Build element list based on the radial location
  !--------------------------------------------------
  allocate(NM%FEM(1)%RAD_ELM_LIST(NM%FEM(1)%NEL,NM%FEM(1)%NGAUSS,5000), &
           NM%FEM(1)%RAD_ELM_NUM( NM%FEM(1)%NEL,NM%FEM(1)%NGAUSS))
  allocate(NM%FEM(2)%RAD_ELM_LIST(NM%FEM(2)%NEL,NM%FEM(2)%NGAUSS,5000), &
           NM%FEM(2)%RAD_ELM_NUM( NM%FEM(2)%NEL,NM%FEM(2)%NGAUSS))
!!!  ! for some reason, this initialization takes forever...
!!!  NM%FEM(1)%RAD_ELM_LIST = -1; NM%FEM(1)%RAD_ELM_NUM = -1
!!!  NM%FEM(2)%RAD_ELM_LIST = -1; NM%FEM(2)%RAD_ELM_NUM = -1

  if (ismaster) write(*,*) "NM: Build element list based on the radial location" 

  call f2f_find_elm_rad(NM%FEM(1), NM%FEM(2), nsd, Center_Rot)
  call f2f_find_elm_rad(NM%FEM(2), NM%FEM(1), nsd, Center_Rot)

  deallocate(NM%FEM(1)%Elm_Loc, NM%FEM(2)%Elm_Loc)

  if (ismaster) then
    write(*,*) "** SHELL ***************************************"
    write(*,*) "total FEM NNODE = ", NM%FEM(:)%NNODE
    write(*,*) "total FEM NEL   = ", NM%FEM(:)%NEL
    write(*,*) "FEM%NGAUSS      = ", NM%FEM(:)%NGAUSS
    write(*,*) "Element Size    = ", NM%FEM(:)%Elm_Size
    write(*,*) "Radial Elem Num = ", maxval(NM%FEM(1)%RAD_ELM_NUM), &
                                     maxval(NM%FEM(2)%RAD_ELM_NUM)
    write(*,*) "************************************************"
  end if

!!!  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!!!  stop

  !--------------------------------------------------------------------
  ! Array for close element list (use around 50 elements for the list)
  !--------------------------------------------------------------------
  NM%FEM(1)%NEL_Close = NM%FEM(2)%NEL*0.003d0
  NM%FEM(2)%NEL_Close = NM%FEM(1)%NEL*0.003d0

  do i = 1, 2
    if (ismaster) then
      write(*,*) "f2f: Number of closest elements:", NM%FEM(i)%NEL_Close
    end if
    allocate(NM%FEM(i)%Elm_Close(NM%FEM(i)%NEL,NM%FEM(i)%NGAUSS,NM%FEM(i)%NEL_Close))
    NM%FEM(i)%Elm_Close = -1
  end do

  !-----------------------------------------------------------
  ! array for closest points (and corresponding element list
  !-----------------------------------------------------------
  do i = 1, 2
    allocate(NM%FEM(i)%CLE(NM%FEM(i)%NEL,NM%FEM(i)%NGAUSS),   &
             NM%FEM(i)%CLP(NM%FEM(i)%NEL,NM%FEM(i)%NGAUSS,2))
    NM%FEM(i)%CLE = 0; NM%FEM(i)%CLP = 0.0d0
  end do

  !----------------------------------------
  ! Read in the closest points
  !----------------------------------------
!!$  call read_close_point_FEM_FEM(NM%FEM(1), NM%FEM(2), NSD)
!!$  call read_close_point_FEM_FEM(NM%FEM(2), NM%FEM(1), NSD)
end subroutine input_shell_nmb




!========================================================================
! Main routine to call all the subroutines   
! Find closest points between FEM and T-Spline
! f2f means for a Gauss point on f1, we want to find the closest point
! on f2
!========================================================================
subroutine read_close_point_FEM_FEM(FEM1, FEM2, NSD)
  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(inout) :: FEM1, FEM2
  integer,    intent(in)    :: NSD

  integer :: i, j, k, ier, nf1, nf2, nf3, itmp1, itmp2

  integer, allocatable :: CLOSE_ELM(:,:,:)
  integer :: NEL_CLOSE

  real(8) :: clock1, clock2


  character(len=30) :: fname1, fname2, cname, ch1, ch2, fmat

  nf1 = 11; nf2 = 21; nf3 = 31

  write(cname,'(I8)') FEM1%FaceID
  fname2 = 'f2f_close_point.'//trim(adjustl(cname))
  write(cname,'(I8)') FEM2%FaceID
  fname2 = trim(adjustl(fname2))//'.'//trim(adjustl(cname))//'.dat'

  open(nf1, file=fname2, status='old', iostat=ier)

  if (ier == 0) then

    if (ismaster) then
      write(*,*) "f2f: Reading the closest-points list", FEM1%FaceID, &
                                                         FEM2%FaceID
    end if

    read(nf1,*) itmp1, itmp2
    if (itmp1/=FEM1%NEL .or. itmp2/=FEM1%NGAUSS) then
      write(*,*) "ERROR: FEM1%NEL or FEM1%NGAUSS does not match the"
      write(*,*) "       numbers in f2f_close_point.dat"
      stop
    end if
    
    do i = 1, FEM1%NEL
      do j = 1, FEM1%NGAUSS
        read(nf1,*) FEM1%CLE(i,j), FEM1%CLP(i,j,:)
      end do
    end do

    ! check if the closest point is between -1 and 1
    if (maxval(abs(FEM1%CLP)) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

  ! if the file does not exist
  else
    write(*,*) "f2f: Closest-points list does not exist. Build it!!!"
    stop
  end if

  close(nf1)
end subroutine read_close_point_FEM_FEM




!=======================================================================
! Subroutine to call all the necessary shell input subroutines
!=======================================================================
subroutine input_shell_blade(NSD, SH, BEA, mNRB_BEA, NRB_BEA)
  use mpi
  use defs_shell
  use types_beam
  use commonpars
  implicit none

  type(shell_bld), intent(out) :: SH
  type(beam), intent(out) :: BEA
  type(mesh_mp_beam), intent(out) :: mNRB_BEA
  type(mesh_beam), intent(out)    :: NRB_BEA

  integer, intent(in) :: NSD

  integer :: so, sf, i, ni, nj, nk, ii, j, k, ip, ier, &
             iel, ct, igauss, jgauss, &
             sumNDZ, sumNEL, nf,  &
             loop2, loop3, lnod, eloc, eglob, &
             p, q, mcp, ncp, nshl, nel, nnode, &
             maxP, maxQ, maxMCP, maxNCP, maxNNODE

  integer :: inewt_SH, iincr_SH, istep, &
             maxP_BEA, maxMCP_BEA, maxNCP_BEA, maxNNODE_BEA, &
             maxP_SHL, maxQ_SHL, maxMCP_SHL, maxNCP_SHL, maxNNODE_SHL, &
             sumNDZ_BEA, sumNEL_BEA, sumNDZ_SHL, sumNEL_SHL, &
             Nincr_step1, Nincr_step2
  real(8) :: RhsNrm0, dis0, dis, sin_theta, cos_theta, theta, length, area, radius, update_fact
  real(8) :: rho_inf, alpha_m, alpha_f, gamma_t, beta_t
  real(8) :: Step1_time, Step2_time, zero_dis
  character(len=30) :: fname41, fname51, fname61, fname, cname

  nf = 21

  ! tolerance of converging RHS
!A  SH%RHSGtol = 1.0d-3

  SH%NMat = 7 !Maximum PTYPE: 2 - Shear Web
  fname = "input.dat"
  open(nf, file=fname, status='old')
  ! read in the number of patches for NUBRS
  read(nf,*) SH%NPT
  read(nf,*)
  read(nf,*)
  read(nf,*)
  read(nf,*)
  read(nf,*)
  read(nf,*)
  read(nf,*)
  read(nf,*)
  read(nf,*) SH%G_fact
  close(nf)

  open(nf, file="input/input_beam.dat", status='old')
  ! read in the total number of patches
  read(nf,*) BEA%NPT
  close(nf)

  write(*,*) BEA%NPT, SH%NPT

  ! Ec longitudinal stiffness for beam
  ! matM(a,3,3): a - material ID (different layups)
  open(nf, file="input/input_beam.mat.dat", status='old')
  read(nf,*) BEA%NMat
  allocate(BEA%Ec(BEA%NMat), BEA%Dens(BEA%NMat), BEA%Rad(BEA%NMat), BEA%C_dp(BEA%NMat) )
  do i = 1, BEA%NMat
    read(nf,*) BEA%Ec(i)
    read(nf,*) BEA%Dens(i)
    read(nf,*) BEA%Rad(i)
    read(nf,*) BEA%C_dp(i)
  end do
  close(nf)


  ! Include gravity load in each direction
  ! (1.0, 0.0, 0.0) means gravity in +x direction
!A  SH%G_fact(1:3) = (/0.0d0, -1.0d0, 0.0d0/)
  BEA%G_fact = SH%G_fact

!  ! Flag for the method (1:NUBRS, 2:T-Spline)
!  read(nf,*) SH%M_Flag

  ! Flag for the thickness (1:Michael, else:Old)
  SH%T_Flag = 0

  !----------------------------------------------------------
  ! use first three processors to solve the blades
  ! since the blades are rotated from the one mike
  ! gave us, need to keep the information of the
  ! rotation for other purpose (thickness, BC, etc..) 
  !----------------------------------------------------------
  if (myid == 0) then
    SH%BldRot = 0.0d0
  else if (myid == 1) then
    SH%BldRot = 120.0d0
  else if (myid == 2) then
    SH%BldRot = 240.0d0
  end if
  ! convert degree to radians
  SH%BldRot = SH%BldRot*pi/180.0d0

  !---------------------------
  ! Preprocessing and Setup
  !---------------------------
  ! A, B and D homogenized material matrix for composite
  ! matM(a,3,3): a - material ID (different layups)
!  open(nf, file="input.sh.mat.dat", status='old')
!  read(nf,*) SH%NMat
!  allocate(SH%matA(SH%NMat,3,3), SH%matB(SH%NMat,3,3), &
!           SH%matD(SH%NMat,3,3))
!  do i = 1, SH%NMat
!    call composite(nf, i, SH%matA(i,:,:), SH%matB(i,:,:), &
!                          SH%matD(i,:,:))
!  end do
!  close(nf)

  if(ismaster) write(*,*) 'get BEA input'

  call beam_input_nrb_mp(BEA%NPT, NSD, maxP_BEA, maxMCP_BEA, &
                           maxNNODE_BEA, NRB_BEA%maxNSHL, mNRB_BEA)

  ! here we use ptype to indicate the material type
  ! need to check if they matches
  if (maxval(mNRB_BEA%PTYPE) > BEA%NMat) then
    write(*,*) "ERROR: Undifined PTYPE!!!", 'Here beam'
    stop
  end if 

  if(ismaster) write(*,*) 'get NRB input'

  call shell_input_nrb_mp(SH%NPT, NSD, maxP, maxQ, maxMCP, &
                          maxNCP, maxNNODE, SH%NRB%maxNSHL, SH%mNRB)

!!$  ! here we use ptype to indicate the material type
!!$  ! need to check if they matches
!A  if (maxval(SH%mNRB%PTYPE) > SH%NMat) then
!A    write(*,*) "ERROR: Undifined PTYPE!!!", 'Here shell'
!A    stop
!A  end if

  ! remove duplcate nodes
  call reduce_node(NSD, BEA%NPT, mNRB_BEA%NNODE, mNRB_BEA%NEL, &
                   mNRB_BEA%B_NET, maxNNODE_BEA, mNRB_BEA%MAP, &
                   NRB_BEA%NNODE, NRB_BEA%NNODE_LOC, sumNDZ_BEA, sumNEL_BEA, &
                   SH%NPT, SH%mNRB%NNODE, SH%mNRB%NEL, &
                   SH%mNRB%B_NET, maxNNODE, SH%mNRB%MAP, &
                   SH%NRB%NNODE, SH%NRB%NNODE_LOC, sumNDZ, sumNEL)
   
  ! number of total elements (non-overlaping)
  SH%NRB%NEL = sum(SH%mNRB%NEL)
  NRB_BEA%NEL = sum(mNRB_BEA%NEL)

  !================================================================

  ! remove patch structure for NURBS
  call beam_input_nrb(mNRB_BEA, NRB_BEA, BEA%NPT, NSD, maxP_BEA, maxMCP_BEA)

  call shell_input_nrb(SH%mNRB, SH%NRB, SH%NPT, NSD, maxP, maxQ, &
                       maxMCP, maxNCP, SH%BldRot)

!!$  ! here we use ptype to indicate the material type
!!$  ! need to check if they matches
  if (maxval(SH%NRB%PTYPE) > SH%NMat) then
    write(*,*) "ERROR: Undifined PTYPE!!!"
    stop
  end if

  ! get the fem shell mesh
  call shell_input_fem_blade(NSD, SH%FEM, SH%BldRot)
  if (maxval(SH%FEM%PTYPE) > SH%NMat) then
    write(*,*) "ERROR: Undifined PTYPE!!!"
    stop
  end if

  if(ismaster) write(*,*) 'combine NURB for beam and shell, and set up BCs'

  call combine_beam_shell_nrb(NSD, NRB_BEA, SH%NRB)

  ! get the t-spline mesh and bezier extraction
!!!  call shell_input_tsp(NSD, SH%TSP, SH%BEZ, SH%BldRot)

!  if (maxval(SH%TSP%PTYPE) > SH%NMat) then
 !   write(*,*) "ERROR: Undifined PTYPE!!!"
  !  stop
  !end if

!!!  call shell_genSparStruc(SH%TSP, NSD, SH)
!A  call shell_genSparStruc(SH%NRB, NSD, SH) !A instead we generate spar matrix for STRU

    NRB_BEA%NGAUSS = maxP_BEA+1
    SH%NRB%NGAUSS = 3    ! in 1D
!!  SH%TSP%NGAUSS = 3    ! in 1D (reduced integration)
    SH%FEM%NGAUSS = 3    ! in 2D (total)
 
!    SH%NRB%B_NET(:,1:3) = SH%NRB%B_NET(:,1:3)*100.0d0/61.0d0



  if (ismaster) then
    write(*,*) "** SHELL ***************************************"
    write(*,*) "total   NNODE = ", sum(SH%mNRB%NNODE)
    write(*,*) "reduced NNODE = ", SH%NRB%NNODE
    write(*,*) "total   NEL   = ", SH%NRB%NEL
    write(*,*)
!!!    write(*,*) "total T-Spline NNODE = ", SH%TSP%NNODE
!!!    write(*,*) "total T-Spline NEL   = ", SH%TSP%NEL
!!!    write(*,*)
    write(*,*) "total FEM NNODE = ", SH%FEM%NNODE
    write(*,*) "total FEM NEL   = ", SH%FEM%NEL
!!$    write(*,*) "NRB%NGAUSS, TSP%NGAUSS=", NRB%NGAUSS, TSP%NGAUSS
    write(*,*)
    write(*,*) "FEM%NGAUSS, NRB%NGAUSS=", SH%FEM%NGAUSS, SH%NRB%NGAUSS
    write(*,*) "************************************************"
  end if

  ! A, B and D homogenized material matrix for composite
  ! matM(a,3,3): a - Gauss point of each NRB element
  allocate(SH%matA(SH%NRB%NEL,SH%NRB%NGAUSS**2,3,3), SH%matB(SH%NRB%NEL,SH%NRB%NGAUSS**2,3,3), &
           SH%matD(SH%NRB%NEL,SH%NRB%NGAUSS**2,3,3), SH%Density(SH%NRB%NEL,SH%NRB%NGAUSS**2), &
           SH%Thickness(SH%NRB%NEL,SH%NRB%NGAUSS**2))

  ! Read Shell material properties:
  open(nf, file = 'Shell_Properties.dat', status = 'old')
  do iel = 1, SH%NRB%NEL
    ct = 0
    do igauss = 1, SH%NRB%NGAUSS
      do jgauss = 1, SH%NRB%NGAUSS
        ct = ct + 1
        read(nf,*) (SH%matA(iel,ct,1:3,j), j = 1, 3), (SH%matB(iel,ct,1:3,j), j = 1, 3), &
                   (SH%matD(iel,ct,1:3,j), j = 1, 3), SH%Density(iel,ct), SH%Thickness(iel,ct)
      end do
    end do
  end do
  close(nf)
 
  ! array for closest points and element list
  allocate(SH%NRB%CLE(SH%NRB%NEL,SH%NRB%NGAUSS**2),   &
           SH%NRB%CLP(SH%NRB%NEL,SH%NRB%NGAUSS**2,2), &
           SH%FEM%CLE(SH%FEM%NEL,SH%FEM%NGAUSS),   &
           SH%FEM%CLP(SH%FEM%NEL,SH%FEM%NGAUSS,2))
  SH%NRB%CLE = 0; SH%NRB%CLP = 0.0d0
  SH%FEM%CLE = 0; SH%FEM%CLP = 0.0d0
!!!  allocate(SH%FEM%CLE(SH%FEM%NEL,SH%FEM%NGAUSS),   &
!!!           SH%FEM%CLP(SH%FEM%NEL,SH%FEM%NGAUSS,2), &
!!!           SH%TSP%CLE(SH%TSP%NEL,SH%TSP%NGAUSS**2),   &
!!!           SH%TSP%CLP(SH%TSP%NEL,SH%TSP%NGAUSS**2,2))
!!!  SH%FEM%CLE = 0; SH%FEM%CLP = 0.0d0
!!!  SH%TSP%CLE = 0; SH%TSP%CLP = 0.0d0

  !----------------------------------------
  ! Find or Read in the closest points
  !----------------------------------------
!!!  call find_close_point_NRB_TSP(NRB, TSP, BEZ, NSD)
!!!  call find_close_point_FEM_TSP(SH%FEM, SH%TSP, SH%BEZ, NSD)
  call find_close_point_nrb_fem(SH%NRB, SH%FEM)
end subroutine input_shell_blade
