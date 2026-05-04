      subroutine Waterline(istep)
      
      use aAdjKeep
      use mpi
      use commonvars

      implicit none
      
      integer :: b,ifac,i,j,wlp,wlf,istep, nc,pc
      real(8) :: xl(NSHLb,NSD),phil(NSHLb)
      real(8) :: xp(NSD),xn(NSD), php, phn, alpha
      real(8),allocatable :: xw(:,:)
      character(len=30) :: fname
      character(len=10) :: cname     

      ! loop over faces
      wlp = 0      
      do b = 1, NBOUND
        if (bound(b)%FACE_ID.eq.7) then
           wlp = wlp + bound(b)%NFACE
        endif
      enddo
      
      allocate(xw(wlp,NSD))
      
      wlp = 0          
      do b = 1, NBOUND
      if (bound(b)%FACE_ID.eq.7) then
      
        do ifac = 1,bound(b)%NFACE
    
       
          do i = 1, NSHLb
            j = bound(b)%FACE_IEN(ifac,i)
            xl(i,:)    = xg(j,:)+dg(j,:)
            phil(i)    = phig(j)
          end do

          if ((maxval(phil).ge.0d0).and.(minval(phil).le.0d0))  then   
	    wlp = wlp + 1
	    
	    php = 0d0
	    phn = 0d0
	    
	    nc = 0
	    pc = 0
	    
	    xp = 0d0
	    xn = 0d0
	    
	    do i = 1, NSHLb
	      if (phil(i).ge.0d0) then
	        pc = pc + 1
	        php = php + phil(i)
	        xp = xp + xl(i,:)
	      else
	        nc = nc + 1
	        phn = phn + phil(i)
	        xn = xn + xl(i,:)
	      endif  
	    enddo 
	    
	    phn = phn/real(nc)
	    xn  = xn /real(nc)
	    
	    php = php/real(pc)
	    xp  = xp /real(pc)

	    alpha     =  php/(php - phn)
	    
	    xw(wlp,:) =  alpha*xn + (1d0-alpha)*xp
	    
          endif
        end do     
      end if
      end do

      wlf = 123
      fname = trim('wl' // cname(istep)) // cname (myid+1)
      open(wlf, file=fname, status='replace')

      write(wlf,*) wlp      
      do i = 1, wlp
        write(wlf,*) xw(i,:)   
      enddo
      close(wlf)
      
      deallocate(xw)

      end subroutine Waterline
