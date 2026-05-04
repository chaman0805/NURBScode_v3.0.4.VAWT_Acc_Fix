
      subroutine solveLevelset(inewt,convres0,converged)
      
      use aAdjKeep
      use mpi
      use commonvars

      implicit none
	  
      real(8) convres0
      logical converged
      
      real(8) ugAlpha(NNODE,NSD),rphigAlpha(NNODE),
     &       phigAlpha(NNODE),dgAlpha(NNODE,NSD),
     &       ugmAlpha(NNODE,NSD)
     
      integer inewt
      real(8) convres, convresl                               
  
      ugAlpha   = ugold   + alfi*(ug - ugold)
      ugmAlpha  = ugmold  + alfi*(ugm - ugmold)
      dgAlpha   = dgold   + alfi*(dg - dgold)
     	 
! Get quantities at alpha levels:	           
       call setBCs()
       phigAlpha  = phigold  + alfi*(phig - phigold)
       rphigAlpha = rphigold + almi*(rphig - rphigold)
	 
! Assemble  
       LHSls  = 0d0     
       RHSGls = 0d0
       call IntElmAss_convdist(
     &        phigAlpha ,rphigAlpha, 
     &        ugAlpha, 
     &        ugmAlpha ,dgAlpha)

       call  FaceAssembly_convdist(
     &        phigAlpha ,ugAlpha ,dgAlpha,ugmAlpha)

! Compute Residual Norm
       if (numnodes.gt.1) call commu(RHSGls,1,'in ')    	                                   
       convres = sum(RHSGls*RHSGls)
         
       if (numnodes.gt.1) then	
	    convresl = convres	 
            call MPI_ALLREDUCE (convresl, convres, 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, 
     &           MPI_COMM_WORLD,mpi_err)               
       endif

       convres = sqrt(convres)
	 	 
! Print Residual Norm
       if ((inewt.eq.1).and.(convres0.lt.0d0)) convres0 = convres
       if (ismaster) then
           write(*,'(I2,x,a,x,ES12.4,x,F12.6)')
     &       inewt, ") Convection Res. Norm = ",
     &      convres,	  1d2*convres/convres0
       endif
	 
! Check Residual Norm
       if (convres.lt.LSC_NL_tol*convres0) then                     
           Converged = .true.       
       else
           Converged = .false.
       endif
       if (isnan(convres)) stop
 	 
! Solve linear system
       if (.not.Converged) then	 
         rphigAlpha   = 0d0         
       
         call SparseGMRES_ls_diag(
     &        LHSls,         
     &        LSC_GMRES_tol, col, row,
     &        RHSGls, 
     &        rphigAlpha, 
     &        LSC_GMRES_itermax, LSC_GMRES_itermin, 
     &        NNODE, maxNSHL, icnt, NSD)  

         rphig   = rphig  + rphigAlpha
         phig    = phig   + gami*Delt*rphigAlpha
       endif

       return
       end
