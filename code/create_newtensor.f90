subroutine create_KDelta(KDelta)
  implicit none
  real(8), intent(out) :: KDelta(3,3)
  integer :: i

KDelta = 0.0d0
do i = 1, 3
    KDelta(i,i) = 1.0d0
enddo
   
end subroutine create_KDelta

subroutine create_PermSymb(PermSymb)
  implicit none
  real(8), intent(out) :: PermSymb(3,3,3)
  integer :: i

PermSymb = 0.0d0
PermSymb(1,2,3) = 1.0d0
PermSymb(2,3,1) = 1.0d0
PermSymb(3,1,2) = 1.0d0
PermSymb(1,3,2) = -1.0d0
PermSymb(2,1,3) = -1.0d0
PermSymb(3,2,1) = -1.0d0
   
end subroutine create_PermSymb



