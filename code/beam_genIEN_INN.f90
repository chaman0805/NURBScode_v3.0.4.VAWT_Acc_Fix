!------------------------------------------------------------------------
! This subroutine generates the IEN matrix, which relates element numbers
! and local node numbers to the appropriate global node numbers. The
! routine also generates the INN matrix, which relates global node
! number to the "NURBS coordinates" of the node.
!------------------------------------------------------------------------
subroutine genIEN_INN_beam(p, nshl, nnodz, nel, mcp, INN, IEN)
  implicit none
  integer, intent(in) :: p, nshl, nnodz, nel, mcp
  integer, intent(out):: INN(nnodz), IEN(nel,nshl)
  integer :: i, k, g, e, gtemp, ln

  ! Loop through control points assigning global node
  ! numbers and filling out IEN and INN as we go
  g = 0
  e = 0
  do i = 1, mcp  ! loop through control points in U direction
    g = g + 1
    INN(g) = i
    if (i >= (p+1)) then
      e = e + 1
        do k = 0, p
          gtemp = g - mcp*0 - k
          ln = (p+1)*0 + k + 1
          IEN(e,ln) = gtemp
        end do
      end if
  end do 
end subroutine genIEN_INN_beam
