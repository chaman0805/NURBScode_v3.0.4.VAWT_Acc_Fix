subroutine genSparStruc_STRU(NEL_BEA, NNODE_BEA, mNSHL_BEA, IEN_BEA, NSHL_BEA, &
                        NEL_SHL, NNODE_SHL, mNSHL_SHL, IEN_SHL, NSHL_SHL, &
                              colm, rowp, icnt)
  implicit none
  integer, intent(in)  :: NEL_BEA, NNODE_BEA, mNSHL_BEA, IEN_BEA(NEL_BEA,mNSHL_BEA), NSHL_BEA(NEL_BEA), &
                          NEL_SHL, NNODE_SHL, mNSHL_SHL, IEN_SHL(NEL_SHL,mNSHL_SHL), NSHL_SHL(NEL_SHL)
  integer, intent(out) :: colm(NNODE_SHL+1), rowp(NNODE_SHL*50*max(mNSHL_BEA, mNSHL_SHL)), icnt
  integer, allocatable :: row_fill_list(:,:)
  integer :: tmpr(NNODE_SHL), adjcnt(NNODE_SHL), mloc(1)
  integer :: i, j, k, imin, ibig, ncol

  allocate(row_fill_list(NNODE_SHL,50*max(mNSHL_BEA, mNSHL_SHL)))

  row_fill_list = 0 
  adjcnt = 0

  !write (*,*) 'icnt =', icnt
  ! compute sparse matrix data structures     
  call Asadj_STRU(NEL_BEA, NNODE_BEA, mNSHL_BEA, IEN_BEA, NSHL_BEA, &
                  NEL_SHL, NNODE_SHL, mNSHL_SHL, IEN_SHL, NSHL_SHL, &
                  row_fill_list, adjcnt)

  ! build the colm array     
  colm(1) = 1
     
  do i = 1, NNODE_SHL
    colm(i+1) = colm(i) + adjcnt(i)
  end do

  ! sort the rowp into increasing order     
  ibig = 10*NNODE_SHL
  ! icnt: total nonzero entries
  icnt = 0
  do i = 1, NNODE_SHL
    ncol = adjcnt(i)
    tmpr(1:ncol) = row_fill_list(i,1:ncol)
    do j = 1, ncol
      icnt = icnt + 1
      imin = minval(tmpr(1:ncol))
      mloc = minloc(tmpr(1:ncol))
      rowp(icnt)    = imin
      tmpr(mloc(1)) = ibig
    end do
  end do

  deallocate(row_fill_list)


end subroutine genSparStruc_STRU




subroutine Asadj_STRU(NEL_BEA, NNODE_BEA, mNSHL_BEA, IEN_BEA, NSHL_BEA, &
                  NEL_SHL, NNODE_SHL, mNSHL_SHL, IEN_SHL, NSHL_SHL, &
                  row_fill_list, adjcnt)
  implicit none
  integer, intent(in)  :: NEL_BEA, NNODE_BEA, mNSHL_BEA, IEN_BEA(NEL_BEA,mNSHL_BEA), NSHL_BEA(NEL_BEA), &
                          NEL_SHL, NNODE_SHL, mNSHL_SHL, IEN_SHL(NEL_SHL,mNSHL_SHL), NSHL_SHL(NEL_SHL)
  integer, intent(out) :: row_fill_list(NNODE_SHL, 50*max(mNSHL_BEA, mNSHL_SHL)), adjcnt(NNODE_SHL)
  integer :: ndlist(max(mNSHL_BEA, mNSHL_SHL))
  integer :: k, i, j, l, ibroke, knd, jnd, jlngth, ni, nj, nk

  do i = 1, NEL_BEA

    ! gen list of global "nodes" for this element
    do j = 1, nshl_BEA(i)
      ndlist(j) = IEN_BEA(i,j)  
    end do
      
    do j = 1, nshl_BEA(i)
      ! jnd is the global "node" we are working on
      jnd = ndlist(j)   
      jlngth = adjcnt(jnd)  ! current length of j's list
      do k = 1, nshl_BEA(i)
        knd = ndlist(k)
        ibroke = 0
      
        ! row_fill_list is, for each node, the
        ! list of nodes that I have already
        ! detected interaction with
        do l= 1,jlngth    
          if (row_fill_list(jnd,l) == knd) then
            ibroke = 1
            exit
          end if
        end do
      
        ! to get here k was not in  j's list so add it     
        if (ibroke == 0) then
          jlngth = jlngth + 1 ! lengthen list
          if (jlngth > 50*nshl_BEA(i)) then
            write(*,*) 'increase overflow factor in genadj'
            stop
          end if        
          row_fill_list(jnd,jlngth)=knd ! add unique entry to list
        end if
      end do    ! finished checking all the k's for this j

      ! update the counter
      adjcnt(jnd)=jlngth  

      end do    ! done with j's
      
!*    endif
    
  end do    ! done with elements in this block


  do i = 1, NEL_SHL

    ! gen list of global "nodes" for this element
    do j = 1, nshl_SHL(i)
      ndlist(j) = IEN_SHL(i,j)  
    end do
      
    do j = 1, nshl_SHL(i)
      ! jnd is the global "node" we are working on
      jnd = ndlist(j)   
      jlngth = adjcnt(jnd)  ! current length of j's list
      do k = 1, nshl_SHL(i)
        knd = ndlist(k)
        ibroke = 0
      
        ! row_fill_list is, for each node, the
        ! list of nodes that I have already
        ! detected interaction with
        do l= 1,jlngth    
          if (row_fill_list(jnd,l) == knd) then
            ibroke = 1
            exit
          end if
        end do
      
        ! to get here k was not in  j's list so add it     
        if (ibroke == 0) then
          jlngth = jlngth + 1 ! lengthen list
          if (jlngth > 50*nshl_SHL(i)) then
            write(*,*) 'increase overflow factor in genadj'
            stop
          end if        
          row_fill_list(jnd,jlngth)=knd ! add unique entry to list
        end if
      end do    ! finished checking all the k's for this j

      ! update the counter
      adjcnt(jnd)=jlngth  

      end do    ! done with j's
      
!*    endif
    
  end do    ! done with elements in this block



end subroutine Asadj_STRU
