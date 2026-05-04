!====================================================================
! program to read in fem mesh (linear triangles so far)
!====================================================================
subroutine shell_input_fem(NSD, FEM)
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: FEM
  integer,    intent(in)    :: NSD
  integer :: i, j, k, ier, mf, P, Q, iel, itmp
  character(len=30) :: fname, cname

  !-----------------------------------------------------
  ! read in FEM format
  !-----------------------------------------------------
  mf = 12
  write(cname,'(I8)') FEM%FaceID
  fname = 'bmesh.'//trim(adjustl(cname))//'.dat'
  open(mf, file=fname, status='old')
     
  read(mf,*) itmp, FEM%maxNSHL
  read(mf,*) FEM%NNODE, FEM%NEL

  allocate(FEM%IEN(FEM%NEL,FEM%maxNSHL), FEM%NSHL(FEM%NEL), &
           FEM%PTYPE(FEM%NEL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: FEM%IEN'
  FEM%IEN   = 0
  FEM%NSHL  = FEM%maxNSHL
  FEM%PTYPE = 1

  allocate(FEM%B_NET(  FEM%NNODE,NSD+1), &
           FEM%B_NET_U(FEM%NNODE,NSD+1), &
           FEM%B_NET_D(FEM%NNODE,NSD+1), &
           FEM%FORCE(  FEM%NNODE,NSD),   &
           FEM%IBC(    FEM%NNODE,NSD), stat=ier)
  if (ier /= 0) stop 'Allocation Error: FEM%B_NET'
  FEM%B_NET = 0.0d0; FEM%B_NET_U = 0.0d0; FEM%B_NET_D = 0.0d0 
  FEM%FORCE = 0.0d0
  FEM%IBC   = 0

  do i = 1, FEM%NNODE
    read(mf,*) (FEM%B_NET(i,j), j = 1, NSD)
  end do
  FEM%B_NET(:,4) = 1.0d0
  FEM%B_NET_U    = FEM%B_NET
  FEM%B_NET_D    = FEM%B_NET

  do iel = 1, FEM%NEL
    read(mf,*) (FEM%IEN(iel,j), j = 1, FEM%NSHL(iel))
  end do

  close(mf)

!!&  ! use the coordinates to setup IBC
!!&  do i = 1, FEM%NNODE
!!&    if (abs(FEM%B_NET(i,2)-2.0d0) < 1.0d-12) then
!!&      FEM%IBC(i,:) = 1
!!&    end if
!!&  end do
end subroutine shell_input_fem



!====================================================================
! program to read in fem mesh (linear triangles so far)
!====================================================================
subroutine shell_input_fem_blade(NSD, FEM, deg)

  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(out) :: FEM
  integer,    intent(in)  :: NSD
  real(8),    intent(in)  :: deg

  integer :: i, j, k, ier, mf, P, Q, iel, itmp
  real(8) :: rtmp, xtmp, utmp
  character(len=30) :: fname, cname

  !-----------------------------------------------------
  ! read in FEM format
  !-----------------------------------------------------
  mf = 12

  write(cname,'(I8)') myid+21
  fname = 'bmesh.'//trim(adjustl(cname))//'.dat'
  open(mf, file=fname, status='old')
     
  read(mf,*) itmp, FEM%maxNSHL
  read(mf,*) FEM%NNODE, FEM%NEL

  allocate(FEM%IEN(FEM%NEL,FEM%maxNSHL), FEM%NSHL(FEM%NEL), &
           FEM%PTYPE(FEM%NEL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: FEM%IEN'
  FEM%IEN   = 0
  FEM%NSHL  = FEM%maxNSHL
  FEM%PTYPE = 1

  allocate(FEM%B_NET(  FEM%NNODE,NSD+1), &
           FEM%B_NET_U(FEM%NNODE,NSD+1), &
           FEM%B_NET_D(FEM%NNODE,NSD+1), &
           FEM%FORCE(  FEM%NNODE,NSD),   &
           FEM%IBC(    FEM%NNODE,NSD), stat=ier)
  if (ier /= 0) stop 'Allocation Error: FEM%B_NET'
  FEM%B_NET = 0.0d0; FEM%B_NET_U = 0.0d0; FEM%B_NET_D = 0.0d0 
  FEM%FORCE = 0.0d0
  FEM%IBC   = 0

  do i = 1, FEM%NNODE
    read(mf,*) (FEM%B_NET(i,j), j = 1, NSD)
    FEM%B_NET(i,1:3) = FEM%B_NET(i,1:3)
  end do
  FEM%B_NET(:,4) = 1.0d0
  FEM%B_NET_U    = FEM%B_NET
  FEM%B_NET_D    = FEM%B_NET

  do iel = 1, FEM%NEL
    read(mf,*) (FEM%IEN(iel,j), j = 1, FEM%NSHL(iel))
  end do

  close(mf)

  ! use the coordinates to setup IBC
!A  do i = 1, FEM%NNODE
 !A   rtmp = FEM%B_NET_U(i,3)!sin(-deg)*FEM%B_NET_U(i,1) + cos(-deg)*FEM%B_NET_U(i,2)
  !A  if (abs(rtmp-10.0d0) < 1.0d-3) then
   !A   FEM%IBC(i,:) = 1
!A    end if
!A  end do

    if(myid == 0) then
      FEM%TipLoc = 400
    else if(myid == 1) then
      FEM%TipLoc = 2139
    else 
      FEM%TipLoc = 6295
    end if

    if(myid == 0) then
      FEM%TipLocTr = 368
    else if(myid == 1) then
      FEM%TipLocTr = 2133
    else 
      FEM%TipLocTr = 6289
    end if

!  write(*,*) "Leading edge Node Number (FEM)=", FEM%TipLoc, &
!             FEM%B_NET_U(FEM%TipLoc,1:2), myid+1

!  write(*,*) "Tailing edge Node Number (FEM)=", FEM%TipLocTr, &
!             FEM%B_NET_U(FEM%TipLocTr,1:2), myid+1
end subroutine shell_input_fem_blade
