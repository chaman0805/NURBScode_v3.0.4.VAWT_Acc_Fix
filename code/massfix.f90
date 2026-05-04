      subroutine massfix(istep)
      
      use aAdjKeep
      use mpi
      use commonvars

      implicit none
      real(8) phip

      integer  inewt,istep
      real(8) ugAlpha(NNODE,NSD),dgAlpha(NNODE,NSD), &
             ugmAlpha(NNODE,NSD)
      real(8) ltmp, rhs0
      real(8) rhs,lhm,mass(2)      
              
      if (ismaster) then
        write(*,*)"##################################################"
        write(*,*) "Levelset mass fix  -- Constant update -- No h"
        write(*,*)"##################################################"
      endif
      
      ugAlpha  = 5d-1*(ug  + ugold)	      
      dgAlpha  = 5d-1*(dg  + dgold)
      ugmAlpha = 5d-1*(ugm + ugmold)
      phip = 0d0         

      do inewt = 1, Mass_NL_itermax
      
      rhs = 0d0
      lhm = 0d0
      	 
! Assemble  

        call IntElmAss_massfix(phip,rhs,lhm,mass)
    
        call FaceAssembly_massfix( &
              dgAlpha, ugAlpha,ugmAlpha, phip,rhs,lhm)    
  
!----------------------------------------------------------------------------
! Compute Momentum Residual Norm
!----------------------------------------------------------------------------
         if (numnodes.gt.1) then
	        ltmp=rhs	
            call MPI_ALLREDUCE (ltmp, rhs, 1, &
                 MPI_DOUBLE_PRECISION,MPI_SUM, &
                 MPI_COMM_WORLD,mpi_err)     
            ltmp=lhm
            call MPI_ALLREDUCE (ltmp, lhm, 1, &
                 MPI_DOUBLE_PRECISION,MPI_SUM, &
                 MPI_COMM_WORLD,mpi_err)                    
         endif

	 if (inewt == 1) rhs0 = abs(rhs)
	 
!---------------------------------------------------------------------------
         if (ismaster) then
	   write(*,'(I2,x,a,x,E12.4,x,F12.6)') &
              inewt, ") Mass. Residual Norm     = ", &
      	       rhs , 1d2*abs(rhs)/rhs0 	         
         endif
	 
!---------------------------------------------------------------------------	 
         phip = phip - rhs/lhm

         if (abs(rhs).lt. Mass_NL_tol*rhs0) exit

      enddo  ! end of newton iter
   
      phig = phig + phip
      rphig   = rphig  + phip/(gami*Delt)
      
      if (Mass_init.le.0d0) Mass_init = mass(1)      

      if (ismaster) then    
!        write(mfile,"(5ES16.8)") time,mass, 
!     &	                         1d2*(mass(2)-Mass_init)/Mass_init,
!     &	                         1d2*(mass(2)-mass(1))/mass(1)
        write(*,"(a,2ES12.4)") "Mass loss    = ", mass(2), &
                                1d2*(mass(2)-Mass_init)/Mass_init   
        write(*,"(a, ES12.4)") "Mass fix phi = ", phip
        write(*,*)"##################################################"
      endif    
        
      if (isnan(rhs)) then
        write(*,*) '!=== Norm is NaN ===!'
        stop
      end if	 
    
      return
      end
      
!##########################################################################
      subroutine IntElmAss_massfix(phip, rhs, lhm, mass)
          
      use aAdjKeep
      use mpi
      use commonvars

      implicit none  
      real(8) phip, rhs, lhm  
	  
	  integer :: nshl

!...  Local variables
      integer iel, igauss, hess_flag,i,idx, NGAUSS
	  
      real(8), allocatable :: shlu(:),shgradgu(:,:), &
                              shhessgu(:,:,:)
	  real(8), allocatable :: dl0(:,:), phil0(:), wl(:), &
                              dl1(:,:), phil1(:), xl(:,:)  	  
			  
			  
      real(8) dxidx(NSD,NSD), DetJ0, DetJ1 
      real(8), allocatable :: gp(:,:), gw(:)

      real(8) phi0,dphi0dxi(NSD),phi1,dphi1dxi(NSD)
      real(8) rhse, lhme
      real(8) eMass_init,emass1, mass(2),ltmp(2)           
      real(8) h,He0,He1,Hep0, Hep1  
      real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)


      

      mass = 0d0
!-----------------------------------------------------------------------------
!     Compute Tangent contribution
!-----------------------------------------------------------------------------
      
!...  loop over elements
      do iel = 1, NELEM
	  
        NSHL   = ELMNSHL(iel)
        NGAUSS = ELMNGAUSS(iel)
        allocate(shlu(NSHL),shgradgu(NSHL,NSD), &
                 shhessgu(NSHL,NSD,NSD),&
                 dl0(NSHL,NSD),phil0(NSHL), wl(NSHL), &
                 dl1(NSHL,NSD),phil1(NSHL), xl(NSHL,NSD),&
                 gp(NGAUSS,NSD), gw(NGAUSS))

        !...  get Gaussian points and weights 
        call genGPandGW(gp,gw,NGAUSS)
	  
	  
        rhse = 0d0
	    lhme = 0d0
	    eMass_init = 0d0
	    emass1 = 0d0          
	    do i = 1, NSHL
	       idx = IEN(iel,i)
               phil0(i) = phigold (idx)
               phil1(i) = phig    (idx) 

	       xl(i,:)  = xg   (idx,:)
               wl(i)    = wg   (idx)
               dl0(i,:) = dgold(idx,:)
	       dl1(i,:) = dg   (idx,:)	  
 	    enddo
	    	    
!...  Loop over integration points (NGAUSS in each direction)
            do igauss = 1, NGAUSS
       
!...  Get Element Shape functions and their gradients
                    
               shlu     = 0d0   ! initialize
               shgradgu = 0d0
               shhessgu = 0d0
                     
               hess_flag = 0   ! Not
         
               call eval_shape(nshl, iel,gp(igauss,:),xl,dl0,wl, shlu, &
			                   shgradgu,shhessgu,dxidx,Gij, Ginv,hess_flag)
     
               DetJ0 = DetJ

               phi0 = sum(phil0*shlu)
               do i = 1, NSD
                 dphi0dxi(i) = sum(phil0*shgradgu(:,i))
	           enddo
	       
               call getElemSize(h,dxidx,dphi0dxi,Ginv)
	          
               call getHeps2(He0,Hep0,phi0,h) 
	       
               call compRhoMu(phi0,dphi0dxi,dxidx,Ginv)	       	       
               eMass_init = eMass_init + rho*DetJ0*gw(igauss)
	       	       	                                       
!...  Get Element Shape functions and their gradients
               shlu     = 0d0   ! initialize
               shgradgu = 0d0
               shhessgu = 0d0
                     
               hess_flag = 0   ! Not
    
               call eval_shape(nshl, iel,gp(igauss,:),xl,dl1,wl, shlu, &
			                   shgradgu,shhessgu,dxidx,Gij, Ginv,hess_flag)

	       DetJ1 = DetJ

               phi1 = sum(phil1*shlu) + phip

               do i = 1, NSD
                 dphi1dxi(i) = sum(phil1*shgradgu(:,i))
	           enddo          

               call getElemSize(h,dxidx,dphi1dxi,Ginv)
	          
               call getHeps2(He1,Hep1,phi1,h) 
	       	       
!...  Calculate residual & Jacobian
	       
               rhse = rhse + (He1*DetJ1-He0*DetJ0)*Dtgl*gw(igauss) 	
               lhme = lhme +  Hep1*Dtgl*DetJ1*gw(igauss)  

               call compRhoMu(phi1,dphi1dxi,dxidx,Ginv)	       	       
               emass1 = emass1 + rho*DetJ1*gw(igauss)

            enddo                          
	  
            mass(1) = mass(1) + eMass_init               	                   
            mass(2) = mass(2) + emass1  
            rhs = rhs + rhse
            lhm = lhm + lhme
			
	    deallocate(shlu,shgradgu, shhessgu)
	    deallocate(dl0, phil0, wl, dl1, phil1, xl, gp, gw)
	  
      enddo
      
      if (numnodes.gt.1) then
            ltmp= mass
            call MPI_ALLREDUCE (ltmp,mass, 2, &
                 MPI_DOUBLE_PRECISION,MPI_SUM, &
                 MPI_COMM_WORLD,mpi_err)
      endif
               
      return
      end

!-----------------------------------------------------------------------------
!     Compute Tangent contribution
!-----------------------------------------------------------------------------
      subroutine FaceAssembly_massfix( &
              dgAlpha, ugAlpha,umgAlpha, phip, rhs, lhm)
      
      use aAdjKeep
      use mpi
      use commonvars
      implicit none
	  
	  integer :: nshl
      
      real(8) phip, rhs, lhm
      real(8) dgAlpha(NNODE,NSD), ugAlpha(NNODE,NSD), &
              umgAlpha(NNODE,NSD)
      
      integer b,ifac,idx,iel,igauss, i

      real(8) gp(NGAUSSb, 2), gw(NGAUSSb), mgp(NGAUSSb,3)
	  
	  
      real(8), allocatable :: phil0(:), phil1(:),  wl(:), &
              ul(:,:), uml(:,:), xl(:,:), dl(:,:)
     
      real(8), allocatable :: shlu(:), shgradgu(:,:)
	  
	  real(8) :: nor(NSD), dxidx(NSD,NSD)
       
      real(8) ui(NSD),umi(NSD),unor,phi,dphidxi(NSD)
           
      real(8) h,He,Hep, rhse, lhme
     
      real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)
 
!     shb will be the shape function array while shbg will hold the
!     gradients of the shape functions
 
      call genGPandGWb(gp,gw,NGAUSSb)
      
!...  Loop over Faces          

      do b = 1, NBOUND
      do ifac = 1,bound(b)%NFACE
      
         call genGPMap(NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)

         iel = bound(b)%F2E(ifac)
		 
		 NSHL = ELMNSHL(iel)
		 allocate(phil0(NSHL), phil1(NSHL),  wl(NSHL), &
              ul(NSHL,NSD),uml(NSHL,NSD), xl(NSHL,NSD), dl(NSHL,NSD))
     
         allocate(shlu(NSHL), shgradgu(NSHL,NSD))
                                                                                
	 do i = 1, NSHL
	    idx = IEN(iel,i)
            phil0(i) = phigold (idx)
            phil1(i) = phig    (idx) 
	    ul(i,:)  = ugAlpha (idx,:)
	    uml(i,:) = umgAlpha(idx,:)
            xl(i,:)  = xg      (idx,:)
            dl(i,:)  = dgAlpha (idx,:)
            wl(i)    = wg      (idx)
 	 enddo
            rhse = 0d0
	    lhme = 0d0              
!... Loop over integration points            
         do igauss = 1, NGAUSSb
                                 
            call eval_faceshape(nshl, iel,gp(igauss,:),mgp(igauss,:), &
                                bound(b)%FACE_OR(ifac), &
                                xl,dl,wl, &
                                shlu,shgradgu,dxidx,Gij,Ginv,nor)
     
!... Interpolate	     
	       
            phi  = 5d-1*(sum(phil0*shlu) + sum(phil1*shlu) + phip)

            do i = 1, NSD
              dphidxi(i) = 5d-1*( sum(phil0*shgradgu(:,i)) &
                                 +sum(phil1*shgradgu(:,i)) )
	    enddo 
	            
	    do i = 1, NSD 
              ui(i)  = sum(ul(:,i) *shlu(:))
              umi(i) = sum(uml(:,i)*shlu(:))
	    enddo
	         
            unor = sum((ui-umi)*nor)  

!... Get element size & interface thickness	
            call getElemSize(h,dxidx,dphidxi,Ginv)	          
            call getHeps2(He,Hep,phi,h) 		  
 
!... Residual & tangent	       
            rhse = rhse + unor*He *DetJb*gw(igauss) 	
            lhme = lhme + unor*Hep*DetJb*gw(igauss)   	    
         enddo
	    
	 rhs = rhs + rhse	
         lhm = lhm + lhme 
		 
      deallocate(phil0, phil1, wl, ul, uml, xl, dl)
     
      deallocate(shlu, shgradgu)
	    
      enddo      
      enddo   
        
      
      return
      end          
