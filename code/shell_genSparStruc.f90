subroutine shell_genSparStruc(FUN, NSD, SH)

  use defs_shell
  implicit none

  type(shell_bld), intent(inout) :: SH
  type(mesh),      intent(in)    :: FUN
  integer,         intent(in)    :: NSD
  integer :: ier

  allocate(SH%col(FUN%NNODE+1), SH%row(FUN%NNODE*50*FUN%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'    
  SH%col = 0; SH%row = 0
    
  call genSparStruc_shell(FUN%NEL, FUN%NNODE, FUN%maxNSHL, FUN%IEN, &
                          FUN%NSHL, SH%col, SH%row, SH%icnt)

end subroutine shell_genSparStruc




subroutine genSparStruc_shell(NEL, NNODE, mNSHL, IEN, NSHL, &
                              colm, rowp, icnt)
  implicit none
  integer, intent(in)  :: NEL, NNODE, mNSHL, IEN(NEL,mNSHL), &
                          NSHL(NEL)
  integer, intent(out) :: colm(NNODE+1), rowp(NNODE*50*mNSHL)
  integer, allocatable :: row_fill_list(:,:)
  integer :: tmpr(NNODE), adjcnt(NNODE), mloc(1)
  integer :: i, j, k, imin, icnt, ibig, ncol

  allocate(row_fill_list(NNODE,50*mNSHL))

  row_fill_list = 0 
  adjcnt = 0

  ! compute sparse matrix data structures     
  call Asadj_shell(NEL, NNODE, mNSHL, IEN, NSHL, &
                   row_fill_list, adjcnt)

  ! build the colm array     
  colm(1) = 1
     
  do i = 1, NNODE
    colm(i+1) = colm(i) + adjcnt(i)
  end do

  ! sort the rowp into increasing order     
  ibig = 10*NNODE
  ! icnt: total nonzero entries
  icnt = 0
  do i = 1, NNODE
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
end subroutine genSparStruc_shell




subroutine Asadj_shell(NEL, NNODE, mNSHL, IEN, NSHL, row_fill_list, adjcnt)
  implicit none
  integer, intent(in)  :: NEL, NNODE, mNSHL, IEN(NEL,mNSHL), NSHL(NEL)
  integer, intent(out) :: row_fill_list(NNODE,50*mNSHL), adjcnt(NNODE)
  integer :: ndlist(mNSHL)
  integer :: k, i, j, l, ibroke, knd, jnd, jlngth, ni, nj, nk

  do i = 1, NEL

    ! gen list of global "nodes" for this element
    do j = 1, nshl(i)
      ndlist(j) = IEN(i,j)  
    end do
      
    do j = 1, nshl(i)
      ! jnd is the global "node" we are working on
      jnd = ndlist(j)   
      jlngth = adjcnt(jnd)  ! current length of j's list
      do k = 1, nshl(i)
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
          if (jlngth > 50*nshl(i)) then
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
end subroutine Asadj_shell
