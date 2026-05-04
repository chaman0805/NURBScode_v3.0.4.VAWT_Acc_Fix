!------------------------------------------------------------------------
! This subroutine generates the IEN matrix, which relates element numbers
! and local node numbers to the appropriate global node numbers. The
! routine also generates the INN matrix, which relates global node
! number to the "NURBS coordinates" of the node.
!------------------------------------------------------------------------
subroutine genIEN_INN_shell(p, q, nshl, nnodz, nel, mcp, ncp, INN, IEN)
  implicit none
  integer, intent(in) :: p, q, nshl, nnodz, nel, mcp, ncp
  integer, intent(out):: INN(nnodz,2), IEN(nel,nshl)
  integer :: i, j, k, l, g, e, gtemp, ln

  ! Loop through control points assigning global node
  ! numbers and filling out IEN and INN as we go
  g = 0
  e = 0
  do j = 1, ncp    ! loop through control points in V direction
    do i = 1, mcp  ! loop through control points in U direction
      g = g + 1
      INN(g,1) = i
      INN(g,2) = j
      if ((i >= (p+1)).and.(j >= (q+1))) then
        e = e + 1
        do l = 0, Q
          do k = 0, p
            gtemp = g - mcp*l - k
            ln = (p+1)*l + k + 1
            IEN(e,ln) = gtemp
          end do
        end do
      end if
    end do
  end do 
end subroutine genIEN_INN_shell
