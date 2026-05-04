!======================================================================
! program to read in t-spline mesh from M.A. Scott
!======================================================================
subroutine shell_input_tsp(NSD, TSP, BEZ, deg)

  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(out) :: TSP, BEZ
  integer,    intent(in)  :: NSD
  real(8),    intent(in)  :: deg

  integer :: i, j, k, ier, mf, P, Q, iel
  character(len=30) :: fname, cname

  !----------------------------------------------------------
  ! read in Mike Scott's extraction format
  !----------------------------------------------------------
  mf = 12
  write(cname,'(I8)') myid + 21
  fname = 'tmesh.'//trim(adjustl(cname))//'.dat'
  open(mf, file=fname, status='old')

  ! The degree of the T-spline basis functions in each 
  ! parametric direction. (not used in the code)
  read(mf,*) P, Q

  ! The number of global T-spline basis functions or control
  ! points in the T-mesh.
  read(mf,*) TSP%NNODE

  ! The number of bezier elements which constitute the T-spline.
  read(mf,*) TSP%NEL

  ! These two entries are added by me for allocating arrays
  ! The maximum number of global T-spline basis functions which 
  ! are non-zero over eacg element
  read(mf,*) TSP%maxNSHL
  ! The maximum number of Bernstein basis functions
  ! For cubic (p=3) functions this number is always 16
  read(mf,*) BEZ%maxNSHL


  allocate(BEZ%P(TSP%NEL), BEZ%Q(TSP%NEL), BEZ%NSHL(TSP%NEL), &
           TSP%NSHL(TSP%NEL), TSP%IEN(TSP%NEL,TSP%maxNSHL),&
           TSP%PTYPE(TSP%NEL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: TSP%IEN'
  BEZ%P = 0; BEZ%Q = 0; BEZ%NSHL = 0
  TSP%NSHL = 0; TSP%IEN = 0; TSP%PTYPE = 1

  allocate(TSP%B_NET(  TSP%NNODE,NSD+1),&
           TSP%B_NET_U(TSP%NNODE,NSD+1),&
           TSP%B_NET_D(TSP%NNODE,NSD+1),&
           TSP%FORCE(  TSP%NNODE,NSD),&
           TSP%IBC(    TSP%NNODE,NSD), stat=ier)
  if (ier /= 0) stop 'Allocation Error: TSP%B_NET'
  TSP%B_NET = 0.0d0; TSP%B_NET_U = 0.0d0; TSP%B_NET_D = 0.0d0 
  TSP%FORCE = 0.0d0
  TSP%IBC   = 0

  allocate(BEZ%Ext(TSP%NEL,TSP%maxNSHL,BEZ%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: BEZ%Ext'
  BEZ%Ext = 0.0d0

  ! The control point (x, y, z, w) associated with each global 
  ! T-spline basis function in the T-mesh. The ith T-spline  denoted
  ! control point is by the token "gi". NOTE: These control points
  ! are NOT in homogeneous form i.e in homogeneous form the control
  ! points have the form (xw, yw, zw, w).
  do i = 1, TSP%NNODE
    read(mf,*) (TSP%B_NET(i,j), j = 1, NSD+1)
  end do
  TSP%B_NET_U = TSP%B_NET
  TSP%B_NET_D = TSP%B_NET

  ! Each bezier element in the T-spline is now enumerated. We only 
  ! annotate the first element as all others follow 
  ! this format exactly.
  do iel = 1, TSP%NEL

    ! The first two integers specify the degree of the bezier 
    ! element in s followed by the degree in t. The third number 
    ! denotes the number of global T-spline basis functions which 
    ! are non-zero over this element. Note that this number often
    ! varies from element to element in a T-spline. 
    ! This is why we need to find maxNSHL_TS in advanced for 
    ! allocating arrays.
    read(mf,*) BEZ%P(iel), BEZ%Q(iel), TSP%NSHL(iel)

    BEZ%NSHL(iel) = (BEZ%P(iel)+1)*(BEZ%Q(iel)+1)

    ! check TSP%maxNSHL and BEZ%maxNSHL
    if (TSP%NSHL(iel) > TSP%maxNSHL) then
      write(*,*) "ERROR: TSP%maxNSHL is wrong!!! Should be:", &
                 TSP%NSHL(iel), TSP%maxNSHL
      stop
    end if

    if (BEZ%NSHL(iel) > BEZ%maxNSHL) then
      write(*,*) "ERROR: BEZ%maxNSHL is wrong!!! Should be:", &
                 BEZ%NSHL(iel), BEZ%maxNSHL
      stop
    end if

    ! For each global T-spline basis function which is non-zero
    ! over this element the global index of the basis function 
    ! is listed.
    read(mf,*) (TSP%IEN(iel,j), j = 1, TSP%NSHL(iel))

    ! The complete extraction operator is now listed for the 
    ! bezier element. Each row corresponds to the decomposition
    ! of a global T-spline basis function into the berstein
    ! basis defined over the element. 
 
    ! The indexing of the bernstein basis functions over the bezier
    ! element where the degrees in s and t are 3 proceeds as 
    ! diagrammed below with B_k = B_[i,j] 
    ! where k = (p + 1) * (i - 1) + j and i,j ranging from 1,..,4.

    ! t
    ! ^ B_13 = B_[4,1] B_14 = B_[4,2] B_15 = B_[4,3] B_16 = B_[4,4]
    ! |
    ! |  B_9 = B_[3,1] B_10 = B_[3,2] B_11 = B_[3,3] B_12 = B_[3,4]
    ! |
    ! |  B_5 = B_[2,1]  B_6 = B_[2,2]  B_7 = B_[2,3]  B_8 = B_[2,4]
    ! |
    ! |  B_1 = B_[1,1]  B_2 = B_[1,2]  B_3 = B_[1,3]  B_4 = B_[1,4]
    ! ---------------> s

    !  In this case, for the first row of the operator we have that 
    !  N_1 = 0.25 * B_1 and for the second row, 
    !  N_2 = 0.25 * B_1 + 0.5 * B_2 + B_3, etc.
    do j = 1, TSP%NSHL(iel)
      read(mf,*) (BEZ%Ext(iel,j,k), k = 1, BEZ%NSHL(iel))
    end do

  end do

  ! read IBC
  do i = 1, TSP%NNODE
    read(mf,*) (TSP%IBC(i,j), j = 1, NSD)
  end do

  close(mf)

  ! find the node number for the tip
  ! rotate the y coord of the blade, then find the location of 
  ! maximum y coord.
  TSP%TipLoc = maxloc((sin(-deg)*TSP%B_NET_U(:,1) + &
                       cos(-deg)*TSP%B_NET_U(:,2)), dim=1)
  write(*,*) "Tip Node Number (TSP)=", TSP%TipLoc, &
             TSP%B_NET_U(TSP%TipLoc,2), myid+1

end subroutine shell_input_tsp
