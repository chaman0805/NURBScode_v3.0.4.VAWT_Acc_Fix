subroutine beam_input_nrb(mNRB, NRB, NPATCH, NSD, maxP, maxMCP, Nincr)

  use types_beam
  implicit none

  type(mesh_mp_beam), intent(in)    :: mNRB
  type(mesh_beam),    intent(inout) :: NRB

  integer, intent(in) :: NPATCH, NSD, maxP, maxMCP, Nincr

  integer :: ier, i, ip, p, mcp, nnode, nel, nshl, eloc, eglob

  ! IEN matches element number a local node number with
  ! patch node number
  integer, allocatable :: IEN_SH(:,:)

  ! INN relate global node number to the (i) "NURBS coordinates"
  integer, allocatable :: INN_SH(:)


  allocate(NRB%P(NRB%NEL), NRB%NSHL(NRB%NEL), &
           NRB%U_KNOT(NRB%NEL,maxMCP+maxP+1), &
           NRB%NUK(NRB%NEL),  &
           NRB%IEN(NRB%NEL,NRB%maxNSHL), NRB%INN(NRB%NEL), &
           NRB%PTYPE(NRB%NEL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: NRB%IEN'
  NRB%P = 0;  NRB%NSHL = 0
  NRB%NUK = 0
  NRB%U_KNOT = 0.0d0
  NRB%IEN = 0

  eglob = 0  
  do ip = NPATCH, 1, -1
    
    !write(*,*) "ip = ", ip
    p     = mNRB%P(ip)
    mcp   = mNRB%MCP(ip)
    nnode = mNRB%NNODE(ip)  ! number of local nodes
    nel   = mNRB%NEL(ip)    ! number of local elements
    nshl  = (p+1)*1     ! number of local shape functions

    allocate(INN_SH(nnode), IEN_SH(nel,nshl), stat=ier)
    if (ier /= 0) stop 'Allocation Error: INN_SH'
    IEN_SH = 0
    INN_SH = 0
  
    ! generate IEN and Coordinates
    call genIEN_INN_beam(p, nshl, nnode, nel, mcp, &
                          INN_SH, IEN_SH)

    do eloc = 1, nel
      !  global element number
      ! write(*,*) "eloc = ", eloc
      eglob = eglob + 1 
      ! write(*,*) "eglob = ", eglob
      ! NRB%IEN and IEN_SH can have different numbers of shape functions
      ! so it is necessary to indicate 1:nshl
      NRB%IEN(eglob,1:nshl) = mNRB%MAP(ip,IEN_SH(eloc,1:nshl))
      NRB%INN(eglob)  = INN_SH(IEN_SH(eloc,1)) 
      ! write(*,*) "IEN_SH = ", IEN_SH(eloc,:)
      ! write(*,*) "INN = ", NRB%INN(eglob)

      ! build the global elements data 
      NRB%P(eglob)    = p
      NRB%NSHL(eglob) = nshl

      NRB%NUK(eglob) = p+mcp+1


      ! every element has a PTYPE, which will be used for
      ! indicating the material type. e.g. if ptype = i,
      ! this element uses the ith material, and ptype = 0
      ! is reserved for the bending strips
      NRB%PTYPE(eglob) = mNRB%PTYPE(ip)

      NRB%U_KNOT(eglob,:) = mNRB%U_KNOT(ip,:)

    end do

    deallocate(IEN_SH)
    deallocate(INN_SH)    
  end do ! End loop over patches

  allocate(NRB%B_NET_D_alphaf(NRB%NNODE,NSD+1), &
           NRB%B_NET_Dt_alphaf(NRB%NNODE,NSD), &
           NRB%B_NET_DDt_alpham(NRB%NNODE,NSD), &
           NRB%B_NET(  NRB%NNODE,NSD+1), &
           NRB%B_NET_U(NRB%NNODE,NSD+1), &
           NRB%B_NET_D(NRB%NNODE,NSD+1), &
           NRB%B_NET_Dt(NRB%NNODE,NSD), &
           NRB%B_NET_DDt(NRB%NNODE,NSD), &
           NRB%FORCE(  NRB%NNODE,NSD), &
           NRB%IBC(    NRB%NNODE,NSD), &
           NRB%B_NET_D_old(NRB%NNODE,NSD+1), &
           NRB%B_NET_Dt_old(NRB%NNODE,NSD), &
           NRB%B_NET_DDt_old(NRB%NNODE,NSD))
  NRB%B_NET_D_alphaf = 0.0d0
  NRB%B_NET_Dt_alphaf = 0.0d0
  NRB%B_NET_DDt_alpham = 0.0d0
  NRB%B_NET   = 0.0d0    ! reference config
  NRB%B_NET_U = 0.0d0    ! Undeformed (used in pre-bend)
  NRB%B_NET_D = 0.0d0    ! current config (deformed)
  NRB%B_NET_Dt = 0.0d0
  NRB%B_NET_DDt = 0.0d0
  NRB%FORCE   = 0.0d0
  NRB%IBC     = 0
  NRB%B_NET_D_old = 0.0d0    
  NRB%B_NET_Dt_old = 0.0d0
  NRB%B_NET_DDt_old = 0.0d0


  ! build the reduced node information
  do ip = 1, NPATCH
    ! do not loop throught bending strips
    if (mNRB%PTYPE(ip) /= 2 .and. mNRB%PTYPE(ip) /= 4) then
      do i = 1, mNRB%NNODE(ip)
        NRB%B_NET(mNRB%MAP(ip,i),:) = mNRB%B_NET(ip,i,:)
        NRB%IBC(  mNRB%MAP(ip,i),:) = mNRB%IBC(ip,i,:)
      end do

      ! For the force, only loop through the blade surface
      ! DO NOT loop through bending strips and SHEAR WEBS
      
      do i = 1, mNRB%NNODE(ip)
        NRB%FORCE(mNRB%MAP(ip,i),:) = mNRB%FORCE(ip,i,:)
        !write (*,*) 'NRB%force=', NRB%FORCE(mNRB%MAP(ip,i),:)
      end do
    end if 
  end do


  NRB%B_NET_U = NRB%B_NET
  !NRB%B_NET_D(1,:,:) = NRB%B_NET(:,1:nsd+1)

  ! use the coordinates to setup IBC_SH because one also
  ! need to take the bending strips into consideration
  !do i = 1, NRB%NNODE
  !  if ((abs(NRB%B_NET(i,1)-0.0d0)<=1.0d-12 .and. abs(NRB%B_NET(i,2)-0.0d0)<=1.0d-12) .and. abs(NRB%B_NET(i,3)-0.0d0)<=1.0d-12) then
  !    NRB%IBC(i,:) = 1
    !elseif (abs(NRB%B_NET(i,1)-0.0d0)<=1.0d-12) then
    !  NRB%IBC(i,1) = 1
    !elseif (abs(NRB%B_NET(i,2)-0.0d0)<=1.0d-12) then
    !  NRB%IBC(i,2) = 1
    !else
    !  NRB%IBC(i,1:2) = 1
  !  end if
  !end do
  


end subroutine beam_input_nrb
