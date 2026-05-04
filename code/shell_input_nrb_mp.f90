!--------------------------------------------------------------------
! program to extract the shell mesh
!--------------------------------------------------------------------
subroutine shell_input_nrb_mp(NP, NSD, maxP, maxQ, maxMCP, maxNCP, &
                              maxNNODE, maxNSHL, mNRB)

  use mpi
  use defs_shell
  implicit none

  type(mesh_mp), intent(out) :: mNRB

  integer, intent(in)  :: NP, NSD
  integer, intent(out) :: maxP, maxQ, maxMCP, maxNCP, maxNNODE, maxNSHL

  integer :: i, j, k, l, mf, ip, ier, tmp, ct

  character(len=30) :: fname, cname, mname

  allocate(mNRB%P(NP),     mNRB%Q(NP), &
           mNRB%MCP(NP),   mNRB%NCP(NP), &
           mNRB%NNODE(NP), mNRB%NEL(NP), &
           mNRB%PTYPE(NP))

  mNRB%P     = 0
  mNRB%Q     = 0
  mNRB%MCP   = 0
  mNRB%NCP   = 0
  mNRB%NNODE = 0
  mNRB%NEL   = 0
  mNRB%PTYPE = -1

  ! first loop through all the patches to find the max number of
  ! parameter. This will be used for allocating other arrays.
  do ip = 1, NP

    mf = 11
  
    ! Read in preliminary information
    write(cname,'(I8)') ip
    write(mname,'(I8)') myid+21
    fname = 'smesh.'//trim(adjustl(mname))//'.'//trim(adjustl(cname))//'.dat'  

    open(mf, file=fname, status='old')
    ! number of spatial dimensions. Usually NSD = 3
    read(mf,*)   
    ! degree of curves in u and v direction
    read(mf,*)  mNRB%P(ip), mNRB%Q(ip)     
    ! number of control points in each direction
    read(mf,*)  mNRB%MCP(ip), mNRB%NCP(ip)

    mNRB%NNODE(ip) = mNRB%MCP(ip)*mNRB%NCP(ip)
    mNRB%NEL(ip)   = (mNRB%MCP(ip)-mNRB%P(ip))*(mNRB%NCP(ip)-mNRB%Q(ip))

    close(mf)
  end do ! end loop patches

  ! these maximum values incluid bending stripes
  maxP     = maxval(mNRB%P)
  maxQ     = maxval(mNRB%Q)
  maxMCP   = maxval(mNRB%MCP)
  maxNCP   = maxval(mNRB%NCP)
  maxNNODE = maxval(mNRB%NNODE)
  maxNSHL  = (maxP+1)*(maxQ+1)

  ! Allocate arrays for knot vectors and for control net
  allocate(mNRB%U_KNOT(NP,maxMCP+maxP+1))
  allocate(mNRB%V_KNOT(NP,maxNCP+maxQ+1))
  mNRB%U_KNOT = 0.0d0
  mNRB%V_KNOT = 0.0d0

  allocate(mNRB%B_NET(NP,maxNNODE,NSD+1))
  allocate(mNRB%IBC(  NP,maxNNODE,NSD))
  mNRB%B_NET  = 0.0d0
  mNRB%IBC    = 0

  ! now loop through the patches again to read in the rest of 
  ! the information
  do ip = 1, NP

    mf = 12  
    ! Read in preliminary information
    write(cname,'(I8)') ip
    write(mname,'(I8)') myid+21
    fname = 'smesh.'//trim(adjustl(mname))//'.'//trim(adjustl(cname))//'.dat'
    open(mf, file=fname, status='old')
    read(mf,'(//)')    ! skip "3" lines

    ! Read knot vectors and control net
    read(mf,*) (mNRB%U_KNOT(ip,i), i = 1, mNRB%MCP(ip)+mNRB%P(ip)+1)
    read(mf,*) (mNRB%V_KNOT(ip,i), i = 1, mNRB%NCP(ip)+mNRB%Q(ip)+1)

    ct = 0
    do j = 1, mNRB%NCP(ip)
      do i = 1, mNRB%MCP(ip)
        ct = ct + 1
        read(mf,*) (mNRB%B_NET(ip,ct,l), l = 1, NSD+1)
        mNRB%B_NET(ip,ct,1:3) = mNRB%B_NET(ip,ct,1:3)
!!!    B_NET_SH(ip,i,j,2) = B_NET_SH(ip,i,j,2) + 2.0d0
      end do
    end do

    ! read in the patch type
    ! 1-blade; 0-bending strips; 2-shear web; ...
    read(mf,*) mNRB%PTYPE(ip) 
    close(mf) 
  end do ! end loop patches


  ! remove duplcate nodes
  ! the reduced node count that counts only master nodes 
  ! in file "input_mp"
  allocate(mNRB%MAP(NP,maxNNODE), stat=ier)
  if (ier /= 0) stop 'Allocation Error: mNRB%MAP'
  mNRB%MAP = -1
end subroutine shell_input_nrb_mp
