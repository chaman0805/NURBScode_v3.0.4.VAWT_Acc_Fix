subroutine IntElmAss_NS_conv(dgAlpha, ugAlpha, ugmAlpha, &
     	                     acgAlpha,acgmAlpha, pgAlpha, &
                             phigAlpha,rphigAlpha)
      
      use aAdjKeep
      use commonvars
      use mpi

      implicit none
	  
      integer :: nshl

      real(8) :: dgAlpha(NNODE,NSD), ugAlpha(NNODE,NSD), &
                 ugmAlpha(NNODE,NSD), acgAlpha(NNODE,NSD), & 
                 acgmAlpha(NNODE,NSD), pgAlpha(NNODE), &
                 phigAlpha(NNODE),rphigAlpha(NNODE)     
          
      ! Local variables
      integer :: iel, igauss,  i, j, k, hess_flag,aa, bb, idx, NGAUSS
            
      real(8), allocatable :: shlu(:), shgradgu(:,:), &
                              shhessgu(:,:,:),shconvggu(:)
				 
      real(8), allocatable :: fl(:,:), dl(:,:), ul(:,:), wl(:), &
                              acl(:,:), uml(:,:), acml(:,:), &
                              pl(:), dlold(:,:), xl(:,:), &
                              phil(:), rphil(:), &
                              philold(:),ulold(:,:)			 
				     
      real(8), allocatable :: xKebe11(:,:,:), xGebe(:,:,:), &
                              xDebe1(:,:,:), xDebe2(:,:,:), xMebe(:,:), &
                              xMebels(:,:), xLSUebe(:,:,:), &
                              xPLSebe(:,:), xULSebe(:,:,:)
      
      real(8), allocatable :: Rhsu(:,:), Rhsm(:,:), Rhsp(:), Rhsls(:)
       
      real(8), allocatable :: gp(:,:), gw(:)
      real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), &
                 fi(NSD), ddidxi(NSD,NSD),phi,dphidxi(NSD), &
                 duidxi(NSD,NSD), dumidxi(NSD,NSD), pri,xi(NSD), &
                 duidxixj(NSD,NSD,NSD),duiolddxi(NSD,NSD),uadvi(NSD), &
                 dxidx(NSD,NSD), dpridxi(NSD), dphiolddxi(NSD), rphi
     
      real(8) :: res,DuDt(NSD),rLi(NSD)
     
      real(8) :: tauM, tauC,tauBar,tau_ls, kdc_ns,kdc_ls(NSD,NSD)
      real(8) :: fact1,fact2,tmpl,cfl(2),gcfl(2), cfl_loc(2)
   
      real(8) :: dphidxi0(NSD),h0,He0,rho0,gwt
      real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)

      
!...  Factors
      fact1 = almi
      fact2 = alfi*gami*Delt 
      cfl = 0d0
      

      
!...  Loop over elements
      do iel = 1, NELEM
        
	 NSHL   = ELMNSHL(iel)
         NGAUSS = ELMNGAUSS(iel)
	 allocate(shlu(NSHL), shgradgu(NSHL,NSD), &
                  shhessgu(NSHL,NSD,NSD),shconvggu(NSHL),&
                  fl(NSHL,NSD), dl(NSHL,NSD), ul(NSHL,NSD),wl(NSHL), &
                  acl(NSHL,NSD), uml(NSHL,NSD), acml(NSHL,NSD), &
                  pl(NSHL), dlold(NSHL,NSD), xl(NSHL,NSD), &
                  phil(NSHL),rphil(NSHL), &
                  philold(NSHL),ulold(NSHL,NSD),&
                  xKebe11(NSD*NSD,NSHL,NSHL), &
                  xGebe(NSD,NSHL,NSHL), xDebe1(NSD,NSHL,NSHL), &
                  xDebe2(NSD,NSHL,NSHL), xMebe(NSHL,NSHL), &
                  xMebels(NSHL,NSHL),  xLSUebe(NSD,NSHL,NSHL), &
                  xPLSebe(NSHL,NSHL),  xULSebe(NSD,NSHL,NSHL),&
                  Rhsu(NSD,NSHL), Rhsm(NSD,NSHL), &
                  Rhsp(NSHL),Rhsls(NSHL),&
                  gp(NGAUSS,NSD), gw(NGAUSS))	
		
       !...  Get Gaussian points and weights             
       call genGPandGW(gp,gw,NGAUSS)

!... Initialize local matrices/vectors
         xKebe11 = 0d0         
         xGebe   = 0d0
         xDebe1  = 0d0
         xDebe2  = 0d0
         xMebe   = 0d0
		
	 xMebels = 0d0
	 xLSUebe = 0d0
	 xPLSebe = 0d0
	 xULSebe = 0d0
        
         Rhsu  = 0d0            
         Rhsm  = 0d0 
         Rhsp  = 0d0
         Rhsls = 0d0
	
!... Get local solution arrays 
         do i = 1, NSHL
	    idx = IEN(iel,i)
	    xl(i,:)    = xg        (idx,:)
	    dl(i,:)    = dgAlpha   (idx,:)
            dlold(i,:) = dgold     (idx,:)
            wl(i)      = wg        (idx)
            ul(i,:)    = ugAlpha   (idx,:)
            ulold(i,:) = ugold     (idx,:)
            acl(i,:)   = acgAlpha  (idx,:)
            uml(i,:)   = ugmAlpha  (idx,:)
            acml(i,:)  = acgmAlpha (idx,:)
            pl(i)      = pgAlpha   (idx)
	    phil(i)    = phigAlpha (idx)
            rphil(i)   = rphigAlpha(idx)
            philold(i) = phigold   (idx)	
         enddo
	
!...  Loop over integration points 	 
         do igauss = 1, NGAUSS
         
!...  Get Element Shape functions and their gradients
            shlu     = 0d+0
            shgradgu = 0d+0
            hess_flag = NS_hess_flag 
            
            call eval_shape(nshl, iel,gp(igauss,:),xl,dl,wl, shlu, &
	                    shgradgu,shhessgu,dxidx,Gij,Ginv, hess_flag)
     
!... Interpolate 
            phi  = sum(phil *shlu)
            rphi = sum(rphil*shlu)
            pri  = sum(pl   *shlu)
	  
            do i = 1, NSD 	    
	       dpridxi(i)    = sum(pl     *shgradgu(:,i))
	       dphidxi(i)    = sum(phil   *shgradgu(:,i))		
	       dphiolddxi(i) = sum(philold*shgradgu(:,i))
	       ui(i)         = sum(ul (:,i)*shlu)
	       umi(i)        = sum(uml(:,i)*shlu)
	       aci(i)        = sum(acl(:,i)*shlu)
	       di(i)         = sum(dl (:,i)*shlu)
	       xi(i)         = sum(xl (:,i)*shlu)
	       do j = 1, NSD
	          duidxi(i,j)    = sum(ul(:,i)   *shgradgu(:,j))
	          duiolddxi(i,j) = sum(ulold(:,i)*shgradgu(:,j))
	       enddo
	
            enddo

!... ALE convection velocity	  
            uadvi = ui - umi 
      	       
!... Forcing       	  
            dphidxi0(1) = 0d0
            dphidxi0(2) = 0d0
            dphidxi0(3) = 1d0 
            call getElemSize(h0,dxidx,dphidxi0,Ginv)	
            call getHeps(He0,water_level-xi(3)-di(3), h0)	  
            rho0 = (1d0-He0)*rhoa + He0*rhow

            call compRhoMu2(phi,dphidxi,dxidx,Ginv)
       
            fi = (rho-rho0)*gravvec
           
!... Residual, tau, kdc for Navier-Stokes             	  
            DuDt(:) = aci(:) + duidxi(:,1)*uadvi(1) &
                      + duidxi(:,2)*uadvi(2) &
                      + duidxi(:,3)*uadvi(3)
     
            rLi(:) = rho*DuDt + dpridxi(:) - fi(:)- &
            mu*(duidxixj(:,1,1)+duidxixj(:,2,2)+duidxixj(:,3,3))
     
            call e3STAB_3D(dxidx, Gij, Ginv, uadvi, rLi, tauM, tauC, &
                           cfl_loc)
	  
            cfl(1)  = max(cfl(1), cfl_loc(1))  
            cfl(2)  = max(cfl(2), cfl_loc(2))  
	      
            call e3DC_beta2(uadvi,duiolddxi,Gij,rLi,dxidx, kdc_ns)	  

!... Residual, tau, kdc for Convection          	  
            res = rphi + sum(uadvi*dphidxi) 

            do aa = 1, NSHL
		shconvggu(aa) = sum(shgradgu(aa,:)*uadvi)
	    enddo 

	    call e3STAB_TAU(uadvi, Gij, tau_ls)  	      
	    call e3DC_CAU(dphiolddxi,uadvi,Gij,res,tau_ls,kdc_ls)
	  
	    kdc_ls = LSC_kdc*kdc_ls  	  
	  	      
!...  Calculate Navier-Stokes residual
	    call e3Rhs_3D_fluid(nshl, ui, umi, aci, acmi, uadvi, pri, rLi, fi, &
		                duidxi, ddidxi, tauM, tauC,  kdc_ns, gw(igauss), &
                                shlu, shgradgu, shgradgu, Rhsu, Rhsp)

!...  Calculate Navier-Stokes Jacobian
            call e3LHS_3D_fluid(nshl, ui, umi, aci, pri, duidxi, dpridxi, rLi, &
                                tauM,  tauC, kdc_ns, gw(igauss), shlu, shgradgu, &
                                shgradgu, shhessgu, xKebe11, xGebe, &
                                xDebe1, xDebe2, xMebe)
     	      
!...  Calculate Convection residual     	  	    
            do aa = 1, NSHL
               rhsls(aa) = rhsls(aa) - (  &
		           (shlu(aa)+shconvggu(aa)*tau_ls)*res &
                           +shgradgu(aa,1)*sum(kdc_ls(1,:)*dphidxi) &
                           +shgradgu(aa,2)*sum(kdc_ls(2,:)*dphidxi) &
                           +shgradgu(aa,3)*sum(kdc_ls(3,:)*dphidxi) &
                           )*DetJ*gw(igauss)
            enddo

!...  Calculate Convection Jacobian             
         do bb = 1, NSHL
            do aa = 1, NSHL
                xMebels(aa,bb) = xMebels(aa,bb) +   (   &
                                (shlu(aa)+shconvggu(aa)*tau_ls)* &
                                (fact1*shlu(bb)+ fact2*shconvggu(bb)) &
                                + fact2*(shgradgu(aa,1)*sum(kdc_ls(1,:)*shgradgu(bb,:)) &
                                + shgradgu(aa,2)*sum(kdc_ls(2,:)*shgradgu(bb,:)) &
                                + shgradgu(aa,3)*sum(kdc_ls(3,:)*shgradgu(bb,:)) &
                                ))*DetJ*gw(igauss)  
            enddo
         enddo 
	  
!...  Calculate NS-Conv Jacobian (For NS Only Galerkin terms)	  
         do bb = 1, NSHL
            do aa = 1, NSHL
               xLSUebe(:,aa,bb) = xLSUebe(:,aa,bb) +   (   &
                                  (shlu(aa)+shconvggu(aa)*tau_ls)*fact2*shlu(bb)*dphidxi &
                                  +  fact2*shgradgu(aa,:)*tau_ls*res &
                                  )*DetJ*gw(igauss)  

               xULSebe(:,aa,bb) = xULSebe(:,aa,bb) +   (   &
                                  shlu(aa)*(DuDt-gravvec)*fact2*drhodphi*shlu(bb) &
                                  )*DetJ*gw(igauss)  

            enddo
         enddo 
      enddo      
	  
!...  Assemble into the Sparse Global Stiffness Matrix and Rhs Vector                   
      call BCLhs_NS_conv(nshl, iel, xKebe11, xGebe, xDebe1, &
                         xDebe2, xMebe, xMebels,xLSUebe,xPLSebe,xULSebe, &
                         Rhsu, Rhsp,RHSls)
	       
      call FillSparseMat_NS_conv(nshl, iel, xKebe11, xGebe, xDebe1, &
                                  xDebe2, xMebe, xMebels,xLSUebe, &
                                  xPLSebe,xULSebe)
	  
      call LocaltoGlobal_NS_conv(nshl, iel, Rhsu, Rhsp, rhsls)

      deallocate(shlu, shgradgu, shhessgu, shconvggu)
				 
      deallocate(fl, dl, ul, wl, acl, uml, acml, &
                 pl, dlold, xl, phil, rphil, &
                 philold, ulold)
			     
      deallocate(xKebe11, xGebe, xDebe1, xDebe2, xMebe, &
                 xMebels, xLSUebe, xPLSebe, xULSebe)
      
      deallocate(Rhsu, Rhsm, Rhsp, Rhsls, gp, gw)

   enddo
      
! ...  Find largest CFL-number and output to screen
      if (numnodes.gt.1) then
          gcfl = cfl
          call MPI_ALLReduce(gcfl, cfl,2, &
     	                     MPI_DOUBLE_PRECISION,MPI_MAX, &
                             MPI_COMM_WORLD,mpi_err) 
      endif

      if (ismaster) then 
         write(*,*)            "-  -  -  -  -  -  -  -  -  -  -  -"
         if ((minval(cfl).gt.0.1d0).and.(maxval(cfl).lt.10d0)) then
            write(*,'(a,x,2F7.4)') "           CFL = ", cfl
         else
            write(*,*)             "           CFL = ", cfl
         endif
      endif 
	  
end subroutine IntElmAss_NS_conv
