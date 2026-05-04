!======================================================================
!
!======================================================================
subroutine redistance(istep)  
  use aAdjKeep
  use mpi
  use commonvars

  implicit none
  
  integer :: istep
  
  real(8) :: phig0(NNODE),dphig(NNODE)  
  integer :: iter
  real(8) :: rdres0,rdres,tmpl
    
  if (ismaster) then
    write(*,*)"##################################################"
    write(*,*) "Redistancing Levelset  --  SUPG"
    write(*,*)"##################################################"
  end if

  phig0 = phig  
    
  do iter = 1, LSRD_NL_itermax 
	 
    ! Assemble     
    LHSls  = 0d0     
    RHSGls = 0d0
    call IntElmAss_redist(phig0)

    ! Compute norm
    if (numnodes .gt. 1) call commu(RHSGls ,1,'in ')     	      
    rdres  = sum(RHSGls *RHSGls )
	 
    if (numnodes.gt.1) then	
      tmpl=rdres	 
      call MPI_ALLREDUCE(tmpl,rdres, 1,    &     
       			 MPI_DOUBLE_PRECISION,MPI_SUM, &
                         MPI_COMM_WORLD,mpi_err)        
    end if
    rdres  = sqrt(rdres)
	 
    if (iter.eq.1) rdres0=rdres
	 
    ! Print norm  	 
    if (ismaster) then       
      write(*,'(I2,x,a,x,ES12.4,x,F12.6)') &
     	     iter,") Redistance: R, |R/R0|(%)= ", &
      rdres,1d2*rdres/rdres0
    end if

    ! Check convergence
    if (rdres.lt.LSRD_NL_tol*rdres0) exit
	 
    ! Solve linear system	 
    dphig  = 0d0	    
   
    call SparseGMRES_ls_diag(LHSls, LSRD_GMRES_tol, col, &
                             row, RHSGls, dphig, &
			     LSRD_GMRES_itermax,LSRD_GMRES_itermin, &
			     NNODE, maxNSHL, icnt, NSD)  

    phig = phig + dphig
    rphig   = rphig  + dphig/(gami*Delt)
		  
  end do
	  
end subroutine redistance
  

!======================================================================
!
!======================================================================
subroutine IntElmAss_redist(phig0)
      
  use aAdjKeep
  use mpi
  use commonvars

  implicit none  
	  
  real(8) :: phig0(NNODE)

  ! Local variables

  integer :: nshl
	  
  integer :: iel, igauss, jgauss, kgauss, hess_flag, &
             ni, nj, nk, i, j, k,idx, aa,bb, NGAUSS
	 
  real(8), allocatable :: shlu(:), shgradgu(:,:), shconvggu(:), &
                          tmp(:),  shhessgu(:,:,:)
	 
  real(8), allocatable :: xMebe(:,:),rhs(:)  
  real(8), allocatable :: phil (:),phi0l(:),dl(:,:), xl(:,:),wl(:)   
	  
  real(8) :: gwt, da,dxidx(NSD,NSD)   
  real(8), allocatable :: gp(:,:), gw(:) 
     
  real(8) :: ai(NSD), phi, phi0, dphidxi (NSD), dphidxi0 (NSD)
  real(8) :: He,Hep,Se, tauM, kdc, res, grab, h, kdc_tmp(NSD,NSD), &
    	     hglobm1,pseudoDTGL
  real(8) :: tmp1(numnodes), tmp2(numnodes)
  real(8) :: nrm, maxnrm(2), minnrm(2),tmpl(2)

  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)
   
  maxnrm = 0d0
  minnrm = 9d9

  NGAUSS = -1
  NSHL   = -1
    
  ! loop over elements
  do iel = 1, NELEM
	    
    if (NSHL /= ELMNSHL(iel)) then     
		 
      if (NSHL >= 0) then		  
        deallocate(shlu, shgradgu, shconvggu, &
           tmp, shhessgu, &
           xMebe,rhs, &
           phil ,phi0l,dl, xl,wl,&
		     gp, gw)      	  
      end if
     
      NSHL   = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)	   
      allocate(shlu(NSHL), shgradgu(NSHL,NSD), shconvggu(NSHL), &
          tmp(NSHL), shhessgu(NSHL,NSD,NSD), &
          xMebe(NSHL,NSHL),rhs(NSHL), &
          phil (NSHL),phi0l(NSHL),dl(NSHL,NSD), xl(NSHL,NSD),wl(NSHL),&
		     gp(NGAUSS,NSD), gw(NGAUSS))   
      call genGPandGW(gp, gw, NGAUSS)
    end if
	
    xMebe = 0d0    
    rhs   = 0d0 
    do i = 1, NSHL
      idx = IEN(iel,i)
      phil(i)  = phig (idx)    	   
      phi0l(i) = phig0(idx) 
      xl(i,:)  = xg(idx,:)
      dl(i,:)  = dg(idx,:)
      wl(i)    = wg(idx)
    end do
	    
    ! Loop over integration points
    do igauss = 1, NGAUSS
   
      ! Get Element Shape functions and their gradients        
      shlu     = 0d0   
      shgradgu = 0d0
      hess_flag = 0

      call eval_shape(nshl, iel,gp(igauss,:),xl,dl,wl, shlu, &
	              shgradgu,shhessgu,dxidx,Gij,Ginv, hess_flag)
     
      ! Interpolate
      phi  = sum(phil*shlu)
      phi0 = sum(phi0l*shlu)
     
      do i = 1, NSD
        dphidxi(i)  = sum(phil*shgradgu(:,i))
        dphidxi0(i) = sum(phi0l*shgradgu(:,i))
      end do 
	   	   	           
      ! Get sign and delta functions	          
      call getElemSize(h,dxidx,dphidxi0,Ginv) ! h based on old phi
      pseudoDTGL=LSRD_pseudoDTGL_fac/h
      call getHeps2(He,Hep,phi0,h)   ! Get Heaviside     
      
      Se   = 2d0*(He-0.5d0) ! Get Sign function

      ! Convective speed
      ai = Se*dphidxi/sqrt(sum(dphidxi*dphidxi)+1d-15)   

      ! Compute tau,kdc	   
      res = Se*(sqrt(sum(dphidxi*dphidxi)) - 1d0)

      call e3STAB_TAU_nt(ai, Gij, tauM)	  	   
      call e3DC_CAUs(dphidxi,ai,Ginv,res,tauM,kdc)

      kdc = LSRD_kdc*kdc

      do aa = 1, NSHL
        shconvggu(aa) = sum(ai*shgradgu(aa,:))
        tmp(aa) = shlu(aa) + tauM*shconvggu(aa)		 	
      end do

      ! Calculate residual   
      do aa = 1, NSHL
         rhs(aa) = rhs(aa) - ( &
                 + tmp(aa)*res &
                 + kdc*sum(shgradgu(aa,:)*dphidxi(:)) &
                 + shlu(aa)*LSRD_penalty_fac*Hep*(phi-phi0) &
                 )*DetJ*gw(igauss)	  
      end do

      ! Calculate Jacobian
      do bb = 1, NSHL
	do aa = 1, NSHL
	  xMebe(aa,bb) = xMebe(aa,bb) +   (  			  &
            tmp(aa)*shconvggu(bb) &
           + kdc*sum(shgradgu(aa,:)*shgradgu(bb,:)) &
           + shlu(aa)*LSRD_penalty_fac*Hep*shlu(bb) &
           )*DetJ*gw(igauss)  
        enddo
     
        xMebe(bb,bb) = xMebe(bb,bb) + &
        shlu(bb)*pseudoDTGL*DetJ*gw(igauss)    ! Lumped Pseudo time relaxation 
      enddo 
	   
      ! Calculate norm 
      nrm = sqrt(sum(dphidxi*dphidxi))
	    
      if (phi.lt.0d0) then
        maxnrm(1) = max(nrm, maxnrm(1))
        minnrm(1) = min(nrm, minnrm(1))
      else
        maxnrm(2) = max(nrm, maxnrm(2))
        minnrm(2) = min(nrm, minnrm(2))	   
      end if
 	       
    end do    
   
    !...  assemble into the Sparse Global Stiffness Matrix and Rhs Vector  
    !    call BCLhs_ls(iel, xMebe, RHS)

    call LocaltoGlobal_ls(nshl, iel, rhs)
    call FillSparseMat_ls(nshl, iel, xMebe)     

  end do
  
  ! ...  Communicate norms
  if (numnodes.gt.1) then
    tmpl=maxnrm
    call MPI_ALLREDUCE(tmpl,maxnrm,2, &
    	       MPI_DOUBLE_PRECISION,MPI_MAX, &
          MPI_COMM_WORLD,mpi_err) 
    tmpl=minnrm
    call MPI_ALLREDUCE(tmpl,minnrm,2, &
          MPI_DOUBLE_PRECISION,MPI_MIN, &
          MPI_COMM_WORLD,mpi_err) 
  end if  
  
  ! Print 
  if (ismaster) then   
    write(*,'(F15.6,a,F9.6)') minnrm(1)," <||grad phi||<",maxnrm(1)
    write(*,'(F15.6,a,F9.6)') minnrm(2)," <||grad phi||<",maxnrm(2)
  end if   
      
  if (NSHL >= 0) then	    
    deallocate(shlu, shgradgu, shconvggu, &
    	   tmp, shhessgu, &
    	   xMebe,rhs, &
    	   phil ,phi0l,dl, xl,wl,&
	  	   gp, gw)		    
  endif	
	    	  
end subroutine IntElmAss_redist
