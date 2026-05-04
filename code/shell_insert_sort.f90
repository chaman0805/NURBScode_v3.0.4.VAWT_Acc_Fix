!======================================================================
!    Insertion Sorting Algorithm
!----------------------------------------------------------------------
!    Sort an array and make the same interchanges in an auxiliary array. 
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      IY - array to be carried with X (all swaps of X elements are
!          matched in IY .  After the sort IY(J) contains the original
!          postition of the value X(J) in the unsorted X array.
!      N - number of values in array X to be sorted
!======================================================================
subroutine In_Sort(X, IY, N)

  implicit none

  ! Scalar Arguments ..
  integer, intent(in) :: N
  ! Array Arguments ..
  real(8), intent(inout) :: X(N)
  integer, intent(inout) :: IY(N)
  ! Local Scalars ..
  real(8) :: TEMP
  integer :: I, J, K, ITEMP

  do I = 2, N
    if ( X(I) < X(I-1) ) then
      do J = I-2, 1, -1
        if (X(I) > X(J)) exit
      end do
      ! note: if the above do loop did not "exit" by the check, 
      ! J = 0, which is exactly what we want.

      TEMP  = X(I)
      ITEMP = IY(I)
        
      do K = I, J+2, -1
        IY(K) = IY(K-1)
        X(K)  = X(K-1)
      end do
      X(J+1)  = TEMP
      IY(J+1) = ITEMP
    end if
  end do
end subroutine In_Sort
