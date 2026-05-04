subroutine SparseCG_BDIAG_shell(mRDNDZ, mNSHL, NSD, icnt,  &
                                colvec, rowvec, LHSK_SH, RHSG_SH, &
                                soln, tol, O_Flag)
  use mpi
  implicit none
  integer, intent(in)  :: mRDNDZ, mNSHL, NSD, icnt, O_Flag, & 
                          colvec(mRDNDZ+1), rowvec(mRDNDZ*50*mNSHL)
  real(8), intent(in)  :: LHSK_SH(NSD*NSD,icnt), RHSG_SH(mRDNDZ,NSD), tol
  real(8), intent(out) :: soln(mRDNDZ,NSD)
  real(8), allocatable :: lhsKBdiag(:,:,:)
  integer :: n, i, j, k, iter
  real(8) :: rhstmp(mRDNDZ,NSD), prodtmp(mRDNDZ,NSD), rr(NSD), &
             rr0, rr1, pp, alpha, tmprr, beta, pv(mRDNDZ,NSD), &
             tauk, taukm1, zs(mRDNDZ,NSD), Binv(NSD,NSD), tmp

  integer, parameter :: maxiter = 3000!200000

!A  write(*,*) RHSG_SH(:,1:3)

  allocate(lhsKBdiag(3,3,mRDNDZ))
  
  ! extract block-diagonal matrix vector
  
  lhsKBdiag = 0.0d0
  do i = 1, mRDNDZ
    do j = colvec(i), colvec(i+1)-1
      n = rowvec(j)
      if (n == i) then
    
        lhsKBdiag(1,1,i) = LHSK_SH(1,j)
        lhsKBdiag(2,2,i) = LHSK_SH(5,j)
        lhsKBdiag(3,3,i) = LHSK_SH(9,j)
    
        lhsKBdiag(1,2,i) = LHSK_SH(2,j)
        lhsKBdiag(1,3,i) = LHSK_SH(3,j)
        lhsKBdiag(2,1,i) = LHSK_SH(4,j)
        lhsKBdiag(2,3,i) = LHSK_SH(6,j)
        lhsKBdiag(3,1,i) = LHSK_SH(7,j)
        lhsKBdiag(3,2,i) = LHSK_SH(8,j)
    
!*    lhsKBdiag(1,1,i) = 1d+0
!*    lhsKBdiag(2,2,i) = 1d+0
!*    lhsKBdiag(3,3,i) = 1d+0    
      end if
    end do
  end do
  
  ! invert block-diagonal
  Binv = 0.0d0     
  do i = 1, mRDNDZ    
    Binv(1,1) = lhsKBdiag(2,2,i)*lhsKBdiag(3,3,i) &
              - lhsKBdiag(3,2,i)*lhsKBdiag(2,3,i)
    Binv(1,2) = lhsKBdiag(3,2,i)*lhsKBdiag(1,3,i) &
              - lhsKBdiag(1,2,i)*lhsKBdiag(3,3,i)
    Binv(1,3) = lhsKBdiag(1,2,i)*lhsKBdiag(2,3,i) &
              - lhsKBdiag(1,3,i)*lhsKBdiag(2,2,i)
    tmp = 1.0d0/( Binv(1,1)*lhsKBdiag(1,1,i)  &
                + Binv(1,2)*lhsKBdiag(2,1,i)  &
                + Binv(1,3)*lhsKBdiag(3,1,i) )
    Binv(1,1) = Binv(1,1)*tmp
    Binv(1,2) = Binv(1,2)*tmp
    Binv(1,3) = Binv(1,3)*tmp
    Binv(2,1) = (lhsKBdiag(2,3,i)*lhsKBdiag(3,1,i) & 
              - lhsKBdiag(2,1,i)*lhsKBdiag(3,3,i))*tmp
    Binv(2,2) = (lhsKBdiag(1,1,i)*lhsKBdiag(3,3,i) &
              - lhsKBdiag(3,1,i)*lhsKBdiag(1,3,i))*tmp
    Binv(2,3) = (lhsKBdiag(2,1,i)*lhsKBdiag(1,3,i) &
              - lhsKBdiag(1,1,i)*lhsKBdiag(2,3,i))*tmp
    Binv(3,1) = (lhsKBdiag(2,1,i)*lhsKBdiag(3,2,i) &
              - lhsKBdiag(2,2,i)*lhsKBdiag(3,1,i))*tmp
    Binv(3,2) = (lhsKBdiag(3,1,i)*lhsKBdiag(1,2,i) &
              - lhsKBdiag(1,1,i)*lhsKBdiag(3,2,i))*tmp
    Binv(3,3) = (lhsKBdiag(1,1,i)*lhsKBdiag(2,2,i) &
              - lhsKBdiag(1,2,i)*lhsKBdiag(2,1,i))*tmp
    
    lhsKBdiag(:,:,i) = Binv(:,:)    
  end do

  rhstmp = 0.0d0  
  do i = 1, 3
    rhstmp(:,i) = RHSG_SH(:,i)
  end do
  
  soln = 0.0d0
  rr   = 0.0d0
  rr1  = 0.0d0  
  do n = 1, mRDNDZ
    rr(:) = rr(:) + rhstmp(n,:)**2
  end do
  
  rr1 = sum(rr)
  rr0 = rr1

  do iter = 1, maxiter
    
    if (ismaster .and. O_Flag == 1) then
      if (mod(iter,1000) == 0) then
        write(*,'(8X,I8,ES14.6)') iter, sqrt(rr1/rr0)
      end if
    end if

    ! Premultiply residual by inverse of Bdiag    
    do i = 1, mRDNDZ
      zs(i,1) = lhsKBdiag(1,1,i)*rhstmp(i,1) +  &
                lhsKBdiag(1,2,i)*rhstmp(i,2) +  &
                lhsKBdiag(1,3,i)*rhstmp(i,3)
      zs(i,2) = lhsKBdiag(2,1,i)*rhstmp(i,1) +  &
                lhsKBdiag(2,2,i)*rhstmp(i,2) +  &
                lhsKBdiag(2,3,i)*rhstmp(i,3)
      zs(i,3) = lhsKBdiag(3,1,i)*rhstmp(i,1) +  &
                lhsKBdiag(3,2,i)*rhstmp(i,2) +  &
                lhsKBdiag(3,3,i)*rhstmp(i,3)
    end do
    
    tauk = 0.0d0
    do i = 1, mRDNDZ
      do j = 1, 3
    tauk = tauk + zs(i,j)*rhstmp(i,j)
      end do
    end do
    
    if (iter == 1) then
      beta = 0.0d0
      pv   = zs
!*      do i = 1,mRDNDZ
!*    write(*,*) pv(i,:)
!*      enddo
    else
      beta = tauk/taukm1
      pv(:,:) = zs(:,:) + beta*pv(:,:)
    end if

    ! Product
    
    call SparseProd_shell(mRDNDZ, mNSHL, NSD, icnt, colvec, &
                          rowvec, LHSK_SH, pv, prodtmp)
         
    pp = 0.0d0
    rr = 0.0d0     
    do n = 1, mRDNDZ
      rr(:) = rr(:) + pv(n,:)*prodtmp(n,:)
    end do
    
    pp = sum(rr)
    
    alpha = tauk/pp
    
    ! calculate the next guess
    do i = 1, 3
      soln(:,i) = soln(:,i) + alpha*pv(:,i)
    end do
    
    ! calculate new res. and dot prod.     
    do i = 1, 3    
      rhstmp(:,i) = rhstmp(:,i) - alpha*prodtmp(:,i)
    end do
    
    tmprr = rr1
    
    rr  = 0.0d0
    rr1 = 0.0d0     
    do n = 1, mRDNDZ
      rr(:) = rr(:) + rhstmp(n,:)**2
    end do
    
    rr1 = sum(rr)	
    
    ! check for convergence 

    taukm1 = tauk

    if (ismaster .and. O_Flag == 1) then
      if (mod(iter,5000) == 0) then
        write(*,'(7X,A,ES14.6,I8)') "res. reduction|iters:", &
                                    sqrt(rr1/rr0), iter
      endif
    end if
    
    if (sqrt(rr1/rr0) < tol) exit

   
  end do
          
  deallocate(lhsKBdiag)
 
  if (ismaster .and. O_Flag == 1) then
    write(*,'(5X,A,ES14.6,I8)') 'SpaCG: Res. reduction|Iters:', &
      sqrt(rr1/rr0), iter
!    write(*,*) 'SpaCG: number of iterations:', iter
!    write(*,*) 'SpaCG:   residual reduction:', sqrt(rr1/rr0)
!    write(*,9000) iter, sqrt(rr1/rr0)
  end if
  
  ! Check if we solved the system....
  ! Product
  call SparseProd_shell(mRDNDZ, mNSHL, NSD, icnt, colvec, &
                        rowvec, LHSK_SH, soln, prodtmp)
  prodtmp = prodtmp - RHSG_SH

  if (ismaster .and. O_Flag == 1) then    
    write(*,'(5X,A,ES14.6)') "SpaCG: Linear algebra check:", &
      sqrt(sum(prodtmp(:,1)**2 + prodtmp(:,2)**2 + prodtmp(:,3)**2))
    write(*,*)
  end if
  
! 9000 format(20x,'  number of iterations:', i10,/, &
!             20x,'    residual reduction:', 2x,e10.2)
  
end subroutine SparseCG_BDIAG_shell
