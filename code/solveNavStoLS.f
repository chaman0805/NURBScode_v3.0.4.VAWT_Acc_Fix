      subroutine solveNavStoLS(inewt,momres0,conres0,NSConverged,
     &	                                convres0,LSConverged)
      
      use aAdjKeep
      use mpi
      use commonvars

      implicit none 

      integer :: inewt
      real(8) :: momres0,conres0,convres0
      logical :: NSConverged, LSConverged
      
      real(8) :: dgAlpha(NNODE,NSD),
     &           ugAlpha(NNODE,NSD), ugmAlpha(NNODE,NSD),
     &           acgAlpha(NNODE,NSD), acgmAlpha(NNODE,NSD),
     &           pgAlpha(NNODE),phigAlpha(NNODE),rphigAlpha(NNODE)
      
      real(8) :: momres,conres,convres,momresl,conresl,convresl
            
 ! Get quantities at alpha levels:
      call setBCs()

      acgAlpha  = acgold  + almi*(acg - acgold)
      acgmAlpha = acgmold + almi*(acgm - acgmold)
      ugAlpha	= ugold   + alfi*(ug - ugold)
      ugmAlpha  = ugmold  + alfi*(ugm - ugmold)
      dgAlpha	= dgold   + alfi*(dg - dgold)
      pgAlpha	= pg
	
      phigAlpha = phigold + alfi*(phig - phigold)
      rphigAlpha = rphigold + almi*(rphig - rphigold)

! Assemble Navier-Stokes/Convection LHS/RHS
      RHSGu  = 0.0d0
      RHSGm  = 0.0d0 
      RHSGp  = 0.0d0 
      RHSGls = 0.0d0 
           
      LHSK11 = 0.0d0
      LHSK12 = 0.0d0
      LHSK22 = 0.0d0
      
      LHSG   = 0.0d0
      LHSD1  = 0.0d0
      LHSD2  = 0.0d0
      LHSM   = 0.0d0
        
      LHSls  = 0.0d0     
      LHSlsu = 0.0d0
      LHSPls = 0.0d0 
      LHSUls = 0.0d0 
	
      call FaceAssembly_NS_conv(dgAlpha, ugAlpha, ugmAlpha, 
     &                          acgAlpha, acgmAlpha, pgAlpha,
     &                          phigAlpha)
        

      call IntElmAss_NS_conv(dgAlpha, ugAlpha, ugmAlpha, 
     &	                     acgAlpha,acgmAlpha, pgAlpha,
     &                       phigAlpha,rphigAlpha)

! Compute residuals               
      if (numnodes.gt.1) then
        call commu(RHSGls,1,   'in ')
        call commu(RHSGp, 1  , 'in ')	  
        call commu(RHSGu, NSD, 'in ')		
      end if

      momresl = sum(RHSGu(:,1)*RHSGu(:,1)) 
     &        + sum(RHSGu(:,2)*RHSGu(:,2)) 
     &        + sum(RHSGu(:,3)*RHSGu(:,3))

      conresl  = sum(RHSGp*RHSGp)
      convresl = sum(RHSGls*RHSGls)
      
      momres  = momresl
      conres  = conresl
      convres = convresl

      if (numnodes.gt.1) then		
         call MPI_ALLREDUCE (momresl, momres, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)
         call MPI_ALLREDUCE (conresl, conres, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)		
         call MPI_ALLREDUCE (convresl, convres, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)  

      endif

      momres  = sqrt(momres)
      conres  = sqrt(conres)
      convres = sqrt(convres)

! Print residual      
      if (inewt.eq.1) then
        momres0  = momres	  
        conres0  = conres 
        if (convres0.lt.0d0) convres0 = convres
      endif
       
      if (ismaster) then
        write(*,'(I2,x,a,x,ES12.4,x,F12.6)')
     &   inewt,") Total Mom. Res. Norm = ", 
     &   momres,	1d2*momres/momres0

        write(*,'(I2,x,a,x,ES12.4,x,F12.6)')
     &   inewt,") Continuity Res. Norm = ",
     &   conres,	1d2*conres/conres0

        write(*,'(I2,x,a,x,ES12.4,x,F12.6)')
     &   inewt,") Convection Res. Norm = ",
     &     convres,1d2*convres/convres0
     
        write(*,*)
      end if

      if (isnan(momres+conres+convres)) then
        write(*,*) '!=== Norm is NaN ===!'
        stop
      end if

! Check convergence
      if  (((momres/momres0).lt.NS_NL_Utol).and.
     &     ((conres/conres0).lt.NS_NL_Ptol) )then
	 NSConverged = .true.
      else
	 NSConverged = .false.      
      endif	
      
      if ((convres/convres0).lt.LSC_NL_tol)  then 
	 LSConverged = .true.
      else
	 LSConverged = .false.      
      endif	 


! Solve appropriate linear system
      if  ((.not.NSConverged).and.(.not.LSConverged)) then
         acgAlpha   = 0d0
         pgAlpha    = 0d0
	 rphigAlpha = 0d0  
	 
     	 call SparseGMRES_NS_conv(col, row,
     &  	IBC, IPER, D_FLAG, P_FLAG, rhsGu, rhsGp, RHSGls,
     &  	acgAlpha, pgAlpha,rphigAlpha,
     &  	lhsK11, lhsG, lhsD1, lhsM,LHSls,LHSlsu,LHSPls,LHSUls,
     &  	icnt, NS_GMRES_tol,NS_GMRES_itermax,NS_GMRES_itermin,
     &  	NNODE, maxNSHL, NSD)

         acg = acg + acgAlpha
         ug  = ug  + gami*Delt*acgAlpha   
         pg  = pg  + alfi*gami*Delt*pgAlpha
	 
         rphig = rphig + rphigAlpha
         phig  = phig  + gami*Delt*rphigAlpha	
      else if (.not.NSConverged) then
         acgAlpha   = 0d0
         pgAlpha    = 0d0
	 
	 call SparseGMRES_up(col, row,
     &  	IBC, IPER, D_FLAG, P_FLAG, rhsGu, rhsGp,
     &  	acgAlpha, pgAlpha,
     &  	lhsK11, lhsG, lhsD1, lhsM,icnt, 
     &  	NS_GMRES_tol, NS_GMRES_itermax,NS_GMRES_itermin,
     &  	NNODE, maxNSHL, NSD)

         acg = acg + acgAlpha
         ug  = ug  + gami*Delt*acgAlpha   
         pg  = pg  + alfi*gami*Delt*pgAlpha
      else if (.not.LSConverged) then
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
      
      if (numnodes > 1) call MPI_BARRIER (MPI_COMM_WORLD,mpi_err)

      end subroutine solveNavStoLS
