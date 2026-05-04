subroutine SparseProd_shell(nnode, nshl, nsd, icnt, colvec, rowvec, &
                            LHS, rhs, prod)
  implicit none
  integer, intent(in)  :: nnode, nshl, nsd, icnt, colvec(nnode+1), &
                          rowvec(nnode*50*nshl)
  real(8), intent(in)  :: LHS(nsd*nsd,icnt), rhs(nnode,nsd)
  real(8), intent(out) :: prod(nnode,nsd)
  integer :: aa, bb, cc
  real(8) :: tmp(nsd)
  
  ! clear the vector
  prod = 0.0d0
  do aa = 1, nnode    
    tmp = 0.0d0     
    do bb = colvec(aa), colvec(aa+1)-1
      cc = rowvec(bb)
      tmp(1) = tmp(1) + LHS(1,bb)*rhs(cc,1) + &
                        LHS(2,bb)*rhs(cc,2) + &
                        LHS(3,bb)*rhs(cc,3)
      tmp(2) = tmp(2) + LHS(4,bb)*rhs(cc,1) + &
                        LHS(5,bb)*rhs(cc,2) + &
                        LHS(6,bb)*rhs(cc,3)
      tmp(3) = tmp(3) + LHS(7,bb)*rhs(cc,1) + &
                        LHS(8,bb)*rhs(cc,2) +&
                        LHS(9,bb)*rhs(cc,3)     
    end do
    prod(aa,:) = prod(aa,:) + tmp(:)
  end do
end subroutine SparseProd_shell


