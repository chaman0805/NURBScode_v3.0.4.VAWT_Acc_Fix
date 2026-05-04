!-----------------------------------------------------------------
! subroutine to compute A and D matrix of symmetric composite
! which will be used to compute mabrain and bending stiffness
! Dm = A*thi, Db = D*thi^3
! thickness is controled outside. Each ply is about 0.1 mm
!-----------------------------------------------------------------  
subroutine composite(mf, iMat, Amat, Bmat, Dmat)
  use commonpars
  implicit none
  integer, intent(in)  :: iMat, mf
  real(8), intent(out) :: Amat(3,3), Bmat(3,3), Dmat(3,3)
  real(8), allocatable :: theta(:), Qbar(:,:,:)
  integer :: i, j, k, iply, NPly, M_Flag
  real(8) :: E1, E2, G12, nu12
  real(8) :: RmatS(3,3), RmatE(3,3), Qmat(3,3), rtmp

  ! read in number of plies
  read(mf,*) NPly

  ! isotropic or orthotropic material
  ! 1-isotropic; 2-orthotropic
  read(mf,*) M_Flag

  if (M_Flag == 1) then
    ! e.g. isotropic case
    ! E1 = 19.0d9; nu12 = 0.29d0
    read(mf,*) E1, nu12
    E2  = E1
    G12 = E1/(2.0d0*(1.0d0+nu12))

  else if (M_Flag == 2) then
    ! e.g. orthotropic case
    ! material data of E-glass/epoxy from Daniel & Ishai (1994)
    ! E1   = 39.0d9    ! fiber direction
    ! E2   =  8.6d9    ! usually epoxy (matrix) effect
    ! G12  =  3.8d9
    ! nu12 = 0.28d0
    read(mf,*) E1, E2, G12, nu12
  
  else
    write(*,*) "ERROR: Wrong Material Flag"
    stop
  end if

  ! read in angle (from x, ccw) of each ply
  ! ply numbers are counted from the bottom!!!!!!!
  ! initialized as zero degree plies
  ! [+-45/0/90_2/0_3]_s
  allocate(theta(Nply), Qbar(3,3,Nply))
  theta = 0.0d0; Qbar = 0.0d0

  read(mf,*) (theta(i), i = 1, Nply)
  theta = theta/180.0d0*pi

  !---------------------------------------------------------------

  ! get the local stiffness matrix
  ! this is also Qbar with theta = 0
  call composite_local_stiff(E1, E2, G12, nu12, Qmat)

  do iply = 1, Nply
    ! get the rotation matrix for stress and strain (global to local)
    ! and their inverse (local to global)
    call composite_rotation_mat( theta(iply), RmatS, RmatE)

    ! compute global stiffness and compliance matrices
    Qbar(:,:,iply) = matmul(transpose(RmatE),matmul(Qmat,RmatE))
  end do

  ! compute A and D matrix (without thickness)
  ! eventually, Dm = A*thi, Db = D*thi^3
  Amat = 0.0d0; Bmat = 0.0d0; Dmat = 0.0d0
  do iply = 1, Nply
    Amat = Amat + Qbar(:,:,iply)

    ! k - n/2 - 1/2
    rtmp = real(iply,8) - 0.5d0*real(Nply,8) - 0.5d0

    Bmat = Bmat + Qbar(:,:,iply)*rtmp

    Dmat = Dmat + Qbar(:,:,iply)*(1.0d0/12.0d0 + rtmp**2)
  end do

  Amat = Amat/dble(Nply)
  Bmat = Bmat/dble(Nply)**2
  Dmat = Dmat/dble(Nply)**3

!  write(*,*) '*** Amat ***************'
!  do i = 1, 3
!    write(*,"(3F8.3)") (Amat(i,j)/1.0d9, j = 1, 3)
!  end do
!  write(*,*) '*** Bmat ***'
!  do i = 1, 3
!    write(*,"(3F8.3)") (Bmat(i,j)/1.0d9, j = 1, 3)
!  end do
!  write(*,*) '*** Dmat ***'
!  do i = 1, 3
!    write(*,"(3F8.3)") (Dmat(i,j)/1.0d9, j = 1, 3)
!  end do

  deallocate(theta, Qbar)

end subroutine composite


!=======================================================================

!----------------------------------------------------------------------
! subroutine to compute local (lamina) stiffness [Q] and compliance [S] 
! matrices relative to the material principal axes (i.e., x1-x2). 
! (based on p.39)
!----------------------------------------------------------------------
subroutine composite_local_stiff(E1, E2, G12, nu12, Qmat)
  implicit none
  real(8), intent( in) :: E1, E2, G12, nu12
  real(8), intent(out) :: Qmat(3,3)
  real(8) :: tmp
  
  ! local stiffness matrix [Q] = [S]^-1
  Qmat = 0.0d0

  tmp = 1.0d0 - nu12*(E2/E1*nu12)
  Qmat(1,1) = E1/tmp
  Qmat(2,2) = E2/tmp
  Qmat(3,3) = G12
  Qmat(1,2) = nu12*E2/tmp
  Qmat(2,1) = Qmat(1,2)
end subroutine composite_local_stiff

!----------------------------------------------------------------------
! subroutine to compute reduced transformation matrix [T_sigma] and 
! [T_epsilon] (p.40)
!----------------------------------------------------------------------
subroutine composite_rotation_mat(the, Tsig, Teps)
  implicit none
  real(8), intent( in) :: the
  real(8), intent(out) :: Tsig(3,3), Teps(3,3)

  ! Rotation matrix for stress (global to local coord.) T_sigma on p.40
  ! For local to global, need Tinv = Tten^-1 = Tten(-the)
  ! Therefore, we can use the same routine, but with -theta,
  ! or just use Transpose(T_epsilon). They should be the same.
  Tsig = 0.0d0

  Tsig(1,1) = cos(the)**2
  Tsig(1,2) = sin(the)**2
  Tsig(1,3) = 2.0d0*sin(the)*cos(the)

  Tsig(2,1) = sin(the)**2
  Tsig(2,2) = cos(the)**2
  Tsig(2,3) = -2.0d0*sin(the)*cos(the)

  Tsig(3,1) = -sin(the)*cos(the)
  Tsig(3,2) =  sin(the)*cos(the)
  Tsig(3,3) = cos(the)**2 - sin(the)**2

  ! Transformation matrix for engineering strain (global to local)
  Teps = 0.0d0

  Teps(1,1) = cos(the)**2
  Teps(1,2) = sin(the)**2
  Teps(1,3) = sin(the)*cos(the)

  Teps(2,1) = sin(the)**2
  Teps(2,2) = cos(the)**2
  Teps(2,3) = -sin(the)*cos(the)

  Teps(3,1) = -2.0d0*sin(the)*cos(the)
  Teps(3,2) =  2.0d0*sin(the)*cos(the)
  Teps(3,3) = cos(the)**2 - sin(the)**2
end subroutine composite_rotation_mat
