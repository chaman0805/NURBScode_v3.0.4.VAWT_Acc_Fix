!======================================================================
! program to read in all patches, remove duplicat nodes and creat
! mapping
! NP1 : total number of patches
! NP0 : number of patches without bending strips
!======================================================================
subroutine reduce_node(NSD, NP1_BEA, NNODZ_BEA, NEL_BEA, B_NET_BEA, maxNNODZ_BEA, mMap_BEA, &
                       redNNODZ_BEA, redNNODZ_LOC_BEA, sumNNODZ_BEA, sumNEL_BEA, &
                        NP1_SHL, NNODZ_SHL, NEL_SHL, B_NET_SHL, maxNNODZ_SHL, mMap_SHL, &
                       redNNODZ_SHL, redNNODZ_LOC_SHL, sumNNODZ_SHL, sumNEL_SHL)

  
  implicit none

  integer, intent( in) :: NSD, NP1_BEA, NNODZ_BEA(NP1_BEA), NEL_BEA(NP1_BEA), maxNNODZ_BEA
  integer, intent(out) :: sumNNODZ_BEA, sumNEL_BEA, redNNODZ_BEA, redNNODZ_LOC_BEA, mMap_BEA(NP1_BEA,maxNNODZ_BEA)
  real(8), intent( in) :: B_NET_BEA(NP1_BEA,maxNNODZ_BEA,NSD+1)
  integer, intent( in) :: NP1_SHL, NNODZ_SHL(NP1_SHL), NEL_SHL(NP1_SHL), maxNNODZ_SHL
  integer, intent(out) :: sumNNODZ_SHL, sumNEL_SHL, redNNODZ_SHL, redNNODZ_LOC_SHL, mMap_SHL(NP1_SHL,maxNNODZ_SHL)
  real(8), intent( in) :: B_NET_SHL(NP1_SHL,maxNNODZ_SHL,NSD+1)

  real(8), parameter :: tol = 1.0d-6

  integer :: i, j, k, l, ip, ier
  integer :: DFlag
  real(8) :: dist
  
  real(8), allocatable :: xg(:,:), xl(:,:)
  
  character(len=30) :: fname
  
  ! determine maximum local NNODZ (for allocating arrays)
  ! also compute the total number of NNODZ (including duplicate nodes)
  sumNNODZ_BEA = sum(NNODZ_BEA)
  sumNNODZ_SHL = sum(NNODZ_SHL)

  ! since elements will not overlap, the total number of elements should
  ! simply be sum(NEL)
  sumNEL_BEA = sum(NEL_BEA)
  sumNEL_SHL = sum(NEL_SHL)
  !--------------------------------------------------------------------
  ! Now, loop through all the patches, read in mesh files, compare
  ! coordinates for overlaping nodes. Remove overlaping nodes and 
  ! creat the mapping.
  ! Note: Loop through NP to 1. i.e., Higher numbering of patches is 
  !   the master to the lower numbering. It has to be this way, 
  !   since in our algorithm, master owns the node and all information
  !   on that node. Slave are zeros on overlapping nodes.
  !--------------------------------------------------------------------

  allocate(xg(sumNNODZ_BEA+sumNNODZ_SHL,NSD+1), stat=ier)
  if(ier /= 0) stop 'Allocation Error: xg'
  xg = 0.0d0

  redNNODZ_BEA = 0
  do ip = NP1_BEA, 1, -1
  
!!!    write(*,*) '====== Patch', ip, beam '======'
  
    allocate(xl(NNODZ_BEA(ip),NSD+1), stat=ier)
    if (ier /= 0) stop 'Allocation Error: U_KNOT'
    xl = 0.0d0
  
    do i = 1, NNODZ_BEA(ip)
      xl(i,:) = B_NET_BEA(ip,i,:)
    end do
  
    if ((i-1)/=NNODZ_BEA(ip)) then
      write(*,*) "issue: NPC*MCP & NNODZ don't match"
    end if
  
    !--------------------------------------------------------------
    ! Now, check for duplicate nodes, or assign to global node
    ! NOTE: always loop backwards so as to maintain master slave relationship,
    !	   in which the master has a higher number than the slave 
    !--------------------------------------------------------------

    do i = 1, NNODZ_BEA(ip)

      DFlag = 0
      do j = 1, redNNODZ_BEA
        ! this should be fine even for the first node
        ! since xg(1:3) could be zero, but xg(4) should not
        dist = sqrt(sum((xg(j,1:3)-xl(i,1:3))**2))
        ! duplicate node found
        if (dist < tol) then
          DFlag = 1
          exit
        end if

      end do ! end loop global non-duplicat nodes (redNNODZ)
  
      ! New node found
      if (DFlag == 0) then
        redNNODZ_BEA = redNNODZ_BEA + 1
        xg(redNNODZ_BEA,:)  = xl(i,:)
        mMap_BEA(ip,i) = redNNODZ_BEA
      ! allocate duplicate node
      else if (DFlag == 1) then
        mMap_BEA(ip,i) = j
      else
        stop 'wrong DFlag'
      end if
    end do ! end loop local nodes (NNODZ(ip))
  
    deallocate(xl)
  end do  ! end loop patches
  redNNODZ_LOC_BEA=redNNODZ_BEA
  redNNODZ_SHL=redNNODZ_BEA

  redNNODZ_LOC_SHL=0
  do ip = NP1_SHL, 1, -1
  
!!!    write(*,*) '====== Patch', ip, shell '======'

  
    allocate(xl(NNODZ_SHL(ip),NSD+1), stat=ier)
    if (ier /= 0) stop 'Allocation Error: U_KNOT'
    xl = 0.0d0
 
    do i = 1, NNODZ_SHL(ip)
      xl(i,:) = B_NET_SHL(ip,i,:)
    end do

  
    if ((i-1)/=NNODZ_SHL(ip)) then
      write(*,*) "issue: NPC*MCP & NNODZ don't match"
    end if
  
    !--------------------------------------------------------------
    ! Now, check for duplicate nodes, or assign to global node
    ! NOTE: always loop backwards so as to maintain master slave relationship,
    !	   in which the master has a higher number than the slave 
    !--------------------------------------------------------------
    do i = 1, NNODZ_SHL(ip)

      DFlag = 0
      do j = 1, redNNODZ_SHL
        ! this should be fine even for the first node
        ! since xg(1:3) could be zero, but xg(4) should not
        dist = sqrt(sum((xg(j,1:3)-xl(i,1:3))**2))
        ! duplicate node found
        if (dist < tol) then
          DFlag = 1
          exit
        end if
      end do ! end loop global non-duplicat nodes (redNNODZ)
  
      ! New node found
      if (DFlag == 0) then
        redNNODZ_LOC_SHL = redNNODZ_LOC_SHL + 1
        redNNODZ_SHL = redNNODZ_SHL + 1
        xg(redNNODZ_SHL,:)  = xl(i,:)
        mMap_SHL(ip,i) = redNNODZ_SHL
      ! allocate duplicate node
      else if (DFlag == 1) then
        mMap_SHL(ip,i) = j
      else
        stop 'wrong DFlag'
      end if
    end do ! end loop local nodes (NNODZ(ip))
  
    deallocate(xl)
  end do  ! end loop patches

  redNNODZ_BEA = redNNODZ_SHL

  deallocate(xg)
end subroutine reduce_node
