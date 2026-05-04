subroutine tensor_prod(vt1,vt2,mat)
  implicit none
  real(8), intent(in)  :: vt1(3), vt2(3) 
  real(8), intent(out) :: mat(3,3)
  integer :: i, j

do i = 1, 3
  do j = 1, 3
    mat(i,j) = vt1(i)*vt2(j)
  enddo
enddo
   
end subroutine tensor_prod

subroutine Mat_prod(vt1,mat2,vt)
  implicit none
  real(8), intent(in)  :: vt1(3), mat2(3,3) 
  real(8), intent(out) :: vt(3)
  integer :: i, j

vt=0.0d0
do i = 1, 3
  do j = 1, 3
    vt(i) = vt(i)+mat2(i,j)*vt1(j)
  enddo
enddo
   
end subroutine Mat_prod


subroutine cross_prod(vt1,vt2,vt)
  implicit none
  real(8), intent(in)  :: vt1(3), vt2(3) 
  real(8), intent(out) :: vt(3)
  
  vt(1) = vt1(2)*vt2(3)-vt1(3)*vt2(2)
  vt(2) = vt2(1)*vt1(3)-vt1(1)*vt2(3)
  vt(3) = vt1(1)*vt2(2)-vt1(2)*vt2(1)
end subroutine cross_prod

subroutine dot_prod(vt1,vt2,sl)
  implicit none
  real(8), intent(in)  :: vt1(3), vt2(3) 
  real(8), intent(out) :: sl
  
  sl = vt1(1)*vt2(1)+vt1(2)*vt2(2)+vt1(3)*vt2(3)
end subroutine dot_prod

