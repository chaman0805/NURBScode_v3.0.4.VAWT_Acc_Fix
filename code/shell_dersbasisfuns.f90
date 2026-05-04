!------------------------------------------------------------------------
! This subroutine evaluates the basis functions and first derivatives 
! functions at a given parameter value u.
!
! Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: 
! Berlin 1995; pp. 72-73.
!------------------------------------------------------------------------
subroutine dersbasisfuns_shell(i, pl, nuk, u, nders, u_knotl, ders)
  
  implicit none

  !--------------VARIABLE DECLARATIONS--------------------------------
  ! knot span, degree of curve, number of control points, counters
  integer :: i, pl, nuk, j, r, k, j1, j2, s1, s2, rk, pk,nders
  ! parameter value, vector of knots, derivative matrix
  real(8) :: u,u_knotl(nuk),ders(nders+1,pl+1),ndu(pl+1,pl+1),d
  real(8) :: left(pl+1), right(pl+1), saved, temp, a(2,pl+1)

  ndu(1,1) = 1
  do j = 1, pl
    left(j+1) = u - u_knotl(i+1-j)
    right(j+1) = u_knotl(i+j) - u
    saved = 0
    do r = 0, j-1
      ndu(j+1,r+1) = right(r+2) + left(j-r+1)
      temp = ndu(r+1,j)/ndu(j+1,r+1)
      ndu(r+1,j+1) = saved + right(r+2)*temp
      saved = left(j-r+1)*temp
    end do
    ndu(j+1,j+1) = saved
  end do
  
  ! load basis functions
  do j = 0, pl
    ders(1,j+1) = ndu(j+1,pl+1)
  end do

  ! compute derivatives
  do r = 0, pl ! loop over function index
    s1 = 0
    s2 = 1         ! alternate rows in array a
    a(1,1) = 1

    ! loop to compute kth derivative
    do k = 1,nders
      d = 0d+0
      rk = r-k
      pk = pl-k
      if (r >= k) then
        a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1)
        d = a(s2+1,1)*ndu(rk+1,pk+1)
      end if
      if (rk >= -1) then
        j1 = 1
      else 
        j1 = -rk
      end if
      if ((r-1) <= pk) then
        j2 = k-1
      else 
        j2 = pl-r
      end if
      do j = j1, j2
        a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1)
        d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1)
      end do
      if (r <= pk) then
        a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1)
        d = d + a(s2+1,k+1)*ndu(r+1,pk+1)
      end if
      ders(k+1,r+1) = d
      j = s1
      s1 = s2
      s2 = j      ! switch rows
    end do
  end do
  
  ! Multiply through by the correct factors
  r = pl
  do k = 1, nders
    do j = 0, pl
      ders(k+1,j+1) = ders(k+1,j+1)*r
    end do
    r = r*(pl-k)
  end do

end subroutine dersbasisfuns_shell
