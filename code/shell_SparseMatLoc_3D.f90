subroutine SparseMatLoc_3D_shell(list, n, goal, locat)  
  implicit none
  integer, intent(in) :: n, list(n), goal
  integer, intent(out):: locat
  integer :: rowvl, rowvh, rowv

  ! Initialize     
  rowvl = 1
  rowvh = n + 1
     
  ! do a binary search
  do 
    if (rowvh-rowvl > 1) then
      rowv = (rowvh + rowvl) / 2
      if (list(rowv) > goal) then
        rowvh = rowv
      else
        rowvl = rowv
      end if
    else
      exit
    end if
  end do
  locat = rowvl    
end subroutine SparseMatLoc_3D_shell
