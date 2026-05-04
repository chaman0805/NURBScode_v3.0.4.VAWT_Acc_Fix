!--------------------------------------------------------------------
! program to extract the beam mesh
!--------------------------------------------------------------------
subroutine beam_input_nrb_mp(NP, NSD, maxP, maxMCP, &
                              maxNNODE, maxNSHL, mNRB)

  use types_beam
  implicit none

  type(mesh_mp_beam), intent(out) :: mNRB

  integer, intent(in)  :: NP, NSD
  integer, intent(out) :: maxP, maxMCP, maxNNODE, maxNSHL

  integer :: i, j, k, l, mf, ip, ier, tmp, ct

  character(len=30) :: fname, cname

  allocate(mNRB%P(NP), &
           mNRB%MCP(NP),  &
           mNRB%NNODE(NP), mNRB%NEL(NP), &
           mNRB%PTYPE(NP))

  mNRB%P     = 0
  mNRB%MCP   = 0
  mNRB%NNODE = 0
  mNRB%NEL   = 0
  mNRB%PTYPE = -1

  ! first loop through all the patches to find the max number of
  ! parameter. This will be used for allocating other arrays.
  do ip = 1, NP

    mf = 11
  
    ! Read in preliminary information
    write(cname,'(I8)') ip
    fname = 'input/smesh_beam.'//trim(adjustl(cname))//'.dat'
  
    open(mf, file=fname, status='old')
    ! number of spatial dimensions. Usually NSD = 3
    read(mf,*)   
    ! degree of curves in u direction
    read(mf,*)  mNRB%P(ip)    
    ! number of control points in u direction
    read(mf,*)  mNRB%MCP(ip)

    mNRB%NNODE(ip) = mNRB%MCP(ip)*1
    mNRB%NEL(ip)   = (mNRB%MCP(ip)-mNRB%P(ip))*1


    close(mf)
  end do ! end loop patches

  ! these maximum values
  maxP     = maxval(mNRB%P)
  maxMCP   = maxval(mNRB%MCP)
  maxNNODE = maxval(mNRB%NNODE)
  maxNSHL  = (maxP+1)*1


  ! Allocate arrays for knot vectors and for control net
  allocate(mNRB%U_KNOT(NP,maxMCP+maxP+1))
  mNRB%U_KNOT = 0.0d0


  allocate(mNRB%B_NET(NP,maxNNODE,NSD+1))
  allocate(mNRB%IBC(NP,maxNNODE,NSD))
  allocate(mNRB%FORCE(NP,maxNNODE,NSD))
  mNRB%B_NET  = 0.0d0
  mNRB%IBC    = 0
  mNRB%FORCE  = 0.0d0
 
  ! now loop through the patches again to read in the rest of 
  ! the information
  do ip = 1, NP

    mf = 12  
    ! Read in preliminary information
    write(cname,'(I8)') ip

    ! write (*,*) ip

    fname = 'input/smesh_beam.'//trim(adjustl(cname))//'.dat'
    open(mf, file=fname, status='old')
    read(mf,'(//)')    ! skip "3" lines

    ! Read knot vectors and control net
    read(mf,*) (mNRB%U_KNOT(ip,i), i = 1, mNRB%MCP(ip)+mNRB%P(ip)+1)
    
    ! write (*,*) (mNRB%U_KNOT(ip,:))  

    ct = 0
   
    do i = 1, mNRB%MCP(ip)
      ct = ct + 1
      read(mf,*) (mNRB%B_NET(ip,ct,l), l = 1, NSD+1)
!!!   B_NET_SH(ip,i,j,2) = B_NET_SH(ip,i,j,2) + 2.0d0
    end do

    ! read in the patch type
    ! 0-cable (only tension or compression); 1-rod...
    read(mf,*) mNRB%PTYPE(ip)

  
    close(mf) 

    ! Read in the force
    ! Notice: the force on the bending strip is zero
    ! When building the reduced node system, be careful not to
    ! overwrite the force on the blade with force on the strip
    write(cname,'(I8)') ip
    fname = 'input/sforce_beam.'//trim(adjustl(cname))//'.dat' 
    open(mf, file=fname, status='old', iostat=ier)

    ! if the force files exist
    if (ier == 0) then
      ! if the patch has fluid force
      if (mNRB%PTYPE(ip) == 1) then
        do i = 1, mNRB%NNODE(ip)
          read(mf,*) (mNRB%FORCE(ip,i,j), j = 1, NSD)
          !write (*,*) 'force=', mNRB%FORCE(ip,i,:)
        end do
      end if
    end if

    close(mf)

  end do ! end loop patches


  ! remove duplcate nodes
  ! the reduced node count that counts only master nodes 
  ! in file "input_mp"
  allocate(mNRB%MAP(NP,maxNNODE), stat=ier)
  if (ier /= 0) stop 'Allocation Error: mNRB%MAP'
  mNRB%MAP = -1
end subroutine beam_input_nrb_mp
