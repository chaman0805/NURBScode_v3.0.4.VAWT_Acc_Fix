!-------------------------------------------------------------------------
! Assemble global stiffness matrix from local stiffness matrix, taking 
! advantage of sparsity. Currently assumes 3d.
!-------------------------------------------------------------------------
subroutine FillSparseMat_3D_shell(nsd, nshl, lIEN, nnode, maxNSHL, &
                                  icnt, col, row, xKebe, LHSK_SH)
  implicit none
  integer, intent(in)    :: nsd, nshl, lIEN(nshl), nnode, maxNSHL, &
                            icnt, row(nnode*50*maxNSHL), col(nnode+1)
  real(8), intent(in)    :: xKebe(NSD*NSD,NSHL,NSHL)
  real(8), intent(inout) :: LHSK_SH(NSD*NSD,icnt)
  integer :: a, b, c, n, k, locn, i
  
  do a = 1, NSHL
    i = lIEN(a)
    c = col(i)
    n = col(i+1) - c
    do b = 1, NSHL
!     call SparseMatLoc_3D_shell(row(c), n, lIEN(b), locn)
      call SparseMatLoc_3D_shell(row(c:c+n-1), n, lIEN(b), locn)
      k = locn + c - 1
      LHSK_SH(:,k) = LHSK_SH(:,k) + xKebe(:,a,b)
    end do
  end do
end subroutine FillSparseMat_3D_shell
