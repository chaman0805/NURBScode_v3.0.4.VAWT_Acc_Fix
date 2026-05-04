subroutine shell_input_nrb(mNRB, NRB, NPATCH, NSD, maxP, maxQ, maxMCP, &
                           maxNCP, deg)

  use mpi
  use defs_shell
  implicit none

  type(mesh_mp), intent(in)    :: mNRB
  type(mesh),    intent(inout) :: NRB

  integer, intent(in) :: NPATCH, NSD, maxP, maxQ, maxMCP, maxNCP
  real(8), intent(in) :: deg

  integer :: ier, i, ip, p, q, mcp, ncp, nnode, nel, nshl, eloc, eglob
  real(8) :: rtmp, xtmp, utmp

  ! IEN matches element number a local node number with
  ! patch node number
  integer, allocatable :: IEN_SH(:,:)

  ! INN relate global node number to the (i,j) "NURBS coordinates"
  integer, allocatable :: INN_SH(:,:)


  allocate(NRB%P(NRB%NEL), NRB%Q(NRB%NEL), NRB%NSHL(NRB%NEL), &
           NRB%U_KNOT(NRB%NEL,maxMCP+maxP+1), &
           NRB%V_KNOT(NRB%NEL,maxNCP+maxQ+1), & 
           NRB%NUK(NRB%NEL), NRB%NVK(NRB%NEL),  &
           NRB%IEN(NRB%NEL,NRB%maxNSHL), NRB%INN(NRB%NEL,2), &
           NRB%PTYPE(NRB%NEL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: NRB%IEN'
  NRB%P = 0; NRB%Q = 0; NRB%NSHL = 0
  NRB%NUK = 0; NRB%NVK = 0;
  NRB%U_KNOT = 0.0d0; NRB%V_KNOT = 0.0d0
  NRB%IEN = 0

  eglob = 0  
  do ip = NPATCH, 1, -1

    p     = mNRB%P(ip)
    q     = mNRB%Q(ip)
    mcp   = mNRB%MCP(ip)
    ncp   = mNRB%NCP(ip)
    nnode = mNRB%NNODE(ip)  ! number of local nodes
    nel   = mNRB%NEL(ip)    ! number of local elements
    nshl  = (p+1)*(q+1)     ! number of local shape functions

    allocate(INN_SH(nnode,2), IEN_SH(nel,nshl), stat=ier)
    if (ier /= 0) stop 'Allocation Error: INN_SH'
    IEN_SH = 0
    INN_SH = 0
  
    ! generate IEN and Coordinates
    call genIEN_INN_shell(p, q, nshl, nnode, nel, mcp, ncp, &
                          INN_SH, IEN_SH)

    do eloc = 1, nel
      !  global element number
      eglob = eglob + 1     
      ! NRB%IEN and IEN_SH can have different numbers of shape functions
      ! so it is necessary to indicate 1:nshl
      NRB%IEN(eglob,1:nshl) = mNRB%MAP(ip,IEN_SH(eloc,1:nshl))
      NRB%INN(eglob,:)  = INN_SH(IEN_SH(eloc,1),:) 

      ! build the global elements data 
      NRB%P(eglob)    = p
      NRB%Q(eglob)    = q
      NRB%NSHL(eglob) = nshl

      NRB%NUK(eglob) = p+mcp+1
      NRB%NVK(eglob) = q+ncp+1

      ! every element has a PTYPE, which will be used for
      ! indicating the material type. e.g. if ptype = i,
      ! this element uses the ith material, and ptype = 0
      ! is reserved for the bending strips
      NRB%PTYPE(eglob) = mNRB%PTYPE(ip)

      NRB%U_KNOT(eglob,:) = mNRB%U_KNOT(ip,:)
      NRB%V_KNOT(eglob,:) = mNRB%V_KNOT(ip,:)
    end do

    deallocate(IEN_SH)
    deallocate(INN_SH)    
  end do ! End loop over patches

  allocate(NRB%B_NET(  NRB%NNODE,NSD+1), &
           NRB%B_NET_U(NRB%NNODE,NSD+1), &
           NRB%B_NET_D(NRB%NNODE,NSD+1), &
           NRB%FORCE(  NRB%NNODE,NSD), &
           NRB%IBC(    NRB%NNODE,NSD))
  NRB%B_NET   = 0.0d0    ! reference config
  NRB%B_NET_U = 0.0d0    ! Undeformed (used in pre-bend)
  NRB%B_NET_D = 0.0d0    ! current config (deformed)
  NRB%FORCE   = 0.0d0
  NRB%IBC     = 0

  ! build the reduced node information
  do ip = 1, NPATCH
    ! do not loop throught bending strips
    if (mNRB%PTYPE(ip) /= 0) then
      do i = 1, mNRB%NNODE(ip)
        NRB%B_NET(mNRB%MAP(ip,i),:) = mNRB%B_NET(ip,i,:)
        NRB%IBC(  mNRB%MAP(ip,i),:) = mNRB%IBC(  ip,i,:)
      end do
    end if
  end do

  NRB%B_NET_U = NRB%B_NET
  NRB%B_NET_D = NRB%B_NET

  ! use the coordinates to setup IBC
!  do i = 1, NRB%NNODE
!    rtmp = NRB%B_NET_U(i,3)!sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2)
!    if (abs(rtmp-10.0d0) < 1.0d-3) then
!      NRB%IBC(i,:) = 1
!    end if
!  end do

  ! find the node number for the tip


!   NRB%TipLoc = 640

!   NRB%TipLocTr = 73

!    rtmp = sin(-deg)*NRB%B_NET_U(1,1) + cos(-deg)*NRB%B_NET_U(1,2)!Y location
 !   xtmp = cos(-deg)*NRB%B_NET_U(1,1) - sin(-deg)*NRB%B_NET_U(1,2) !X coordinate
  !  NRB%TipLoc = 1
   ! NRB%TipLocTr = 1
    !do i = 1, NRB%NNODE
 !!     if(abs(sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2) - rtmp) > 1.0d-4) then
  !    if((sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2)) > rtmp) then
   !     NRB%TipLoc = i
    !    NRB%TipLocTr = i
     !   rtmp = sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2)
      !  xtmp = cos(-deg)*NRB%B_NET_U(i,1) - sin(-deg)*NRB%B_NET_U(i,2)
       ! utmp = xtmp
  !!    else if(abs(sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2) -  rtmp) <= 1.0d-4)  then
  !    else if((sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2)) ==  rtmp) then
  !!      if((cos(-deg)*NRB%B_NET_U(i,1) - sin(-deg)*NRB%B_NET_U(i,2) - xtmp) <= 1.0d-4) then
    !    if((cos(-deg)*NRB%B_NET_U(i,1) - sin(-deg)*NRB%B_NET_U(i,2)) <= xtmp) then
     !     NRB%TipLoc = i
      !    rtmp = sin(-deg)*NRB%B_NET_U(i,1) + cos(-deg)*NRB%B_NET_U(i,2)
       !   xtmp = cos(-deg)*NRB%B_NET_U(i,1) - sin(-deg)*NRB%B_NET_U(i,2)
!        end if
 !       if(NRB%B_NET_U(i,1) > utmp) then
  !        NRB%TipLocTr = i
   !       rtmp = NRB%B_NET_U(i,2)
    !      utmp = NRB%B_NET_U(i,1)
     !   end if
      !end if
 !   end do

!  write(*,*) "Leading edge Node Number (NRB)=", NRB%TipLoc, &
!             NRB%B_NET_U(NRB%TipLoc,1:2), myid+1

!  write(*,*) "Trailing edge Node Number (NRB)=", NRB%TipLocTr, &
!             NRB%B_NET_U(NRB%TipLocTr,1:2), myid+1
end subroutine shell_input_nrb
