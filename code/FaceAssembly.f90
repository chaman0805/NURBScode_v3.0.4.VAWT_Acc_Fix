!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                           acgmAlpha, pgAlpha, phigAlpha)
  
   use aAdjKeep
   use mpi
   use commonvars

   implicit none
  
   integer :: nshl
          
   real(8), intent(in) :: dgAlpha(NNODE,NSD),   ugAlpha(NNODE,NSD), &
                          acgAlpha(NNODE,NSD),  ugmAlpha(NNODE,NSD), &
                          acgmAlpha(NNODE,NSD), pgAlpha(NNODE), &
                          phigAlpha(NNODE)

   ! Local variables
   integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb
  
   real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)


   real(8) :: gwt, rDetJb, fact1, fact2

   real(8) :: dphidx(NSD), dxidx(NSD,NSD)
  
   real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:)
                  
   real(8), allocatable :: dl(:,:), ul(:,:), acl(:,:), uml(:,:), &
                           acml(:,:), pl(:), dlold(:,:), xl(:,:), &
                           phil(:), plold(:), wl(:)
			 
   real(8), allocatable :: xKebe11(:,:,:), xMebe(:,:), xGebe(:,:,:), &
                           xDebe1(:,:,:), xDebe2(:,:,:)
			 
   real(8), allocatable :: Rhsu(:,:), Rhsp(:)     
     
   real(8) :: ui(NSD), pri, duidxi(NSD,NSD), phi, dphidxi(NSD), &
              Gij(NSD,NSD), Ginv(NSD,NSD), ti(NSD)
   real(8) :: nor(NSD), bui(NSD), umi(NSD)
  
   real(8) :: tauB, tauNor, gi(NSD), unor, upos, uneg  
   real(8) :: ixint(NSD), bxint(NSD), atime  
   real(8) :: gneg, wave_u, wave_phi

!   real(8), allocatable :: gp(:,:), gw(:), mgp(:,:)
    
   ! "shb" will be the shape function array while "shbg" will hold the
   ! gradients of the shape functions  
   fact1 = almi
   fact2 = alfi*gami*Delt

   gi = 0.0d0
        
   atime = time -(1.0d0-alfi)*Delt

   NSHL = -1

   ! loop over faces
   do b = 1, NBOUND

      do ifac = 1, bound(b)%NFACE

         iel = bound(b)%F2E(ifac)

         if (NSHL.ne.ELMNSHL(iel)) then
            if (NSHL.ge.0) then
               deallocate(shlu, shgradgu, shhessgu)
               deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, plold, wl)
               deallocate(xKebe11, xMebe, xGebe, xDebe1, xDebe2)
               deallocate(Rhsu, Rhsp)
            endif
  
!            call genGPMap(NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)

            if (ifac==10) then
               write(*,*) bound(b)%FACE_OR(ifac)
            endif
 

            NSHL = ELMNSHL(iel)

            allocate(shlu(NSHL), shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD))
                 
            allocate(dl(NSHL,NSD), ul(NSHL,NSD), acl(NSHL,NSD), uml(NSHL,NSD), &
                     acml(NSHL,NSD), pl(NSHL), dlold(NSHL,NSD), xl(NSHL,NSD), &
                     phil(NSHL), plold(NSHL), wl(NSHL))
			 
            allocate(xKebe11(NSD*NSD,NSHL,NSHL), xMebe(NSHL,NSHL), &
                     xGebe(NSD,NSHL,NSHL), xDebe1(NSD,NSHL,NSHL), &
                     xDebe2(NSD,NSHL,NSHL))
			 
            allocate(Rhsu(NSD,NSHL), Rhsp(NSHL))

            xDebe2  = 0d0          

            call genGPandGWb(gp, gw, NGAUSSb)   	  
         end if


         ! Get local solution arrays   	       	     
         do i = 1, NSHL
            j = IEN(iel,i)
            xl   (i,:) = xg       (j,:)
            dl   (i,:) = dgAlpha  (j,:)
            dlold(i,:) = dgold    (j,:)
            ul   (i,:) = ugAlpha  (j,:)
            acl  (i,:) = acgAlpha (j,:)
            uml  (i,:) = ugmAlpha (j,:)
            acml (i,:) = acgmAlpha(j,:)
            pl   (i)   = pgAlpha  (j)
            plold(i)   = pgold    (j)
            phil (i)   = phigAlpha(j)
            wl   (i)   = wg       (j)
         end do
     
         ! initialize local resistance matrix	
         xKebe11 = 0.0d0         
         xGebe   = 0.0d0
         xDebe1  = 0.0d0
         xMebe   = 0.0d0
    
         ! initialize local load vector
         Rhsu    = 0.0d0    
         Rhsp    = 0.0d0   
  
         ! Loop over integration points (NGAUSSb in each direction) 
         do igauss = 1, NGAUSSb
 
            call eval_faceshape(nshl, iel, gp(igauss,:), mgp(igauss,:), &
                                bound(b)%FACE_OR(ifac), &
                                xl, dl, wl, shlu, shgradgu, dxidx, &
                                Gij, Ginv, nor)
          
            ! Interpolate	    
            pri = sum(pl*shlu)
            phi = sum(phil*shlu)

            do i = 1, NSD
               ui   (i) = sum(ul(:,i)*shlu)
               umi  (i) = sum(uml(:,i)*shlu)
               ixint(i) = sum((xl(:,i)+dl(:,i))*shlu)
            end do     
    
         do j = 1, NSD
            dphidxi(j) = sum(phil*shgradgu(:,j))
            do i = 1, NSD
               duidxi(i,j) = sum(ul(:,i)*shgradgu(:,j))
            end do
         end do

        ! Secondary quantities      
!!!        call compRhoMu2(phi, dphidxi, dxidx, Ginv)
         rho = rhow
         mu  = muw
    
         unor = sum((ui-umi)*nor(:))  
         upos = 0.5d0*(unor+abs(unor))
         uneg = 0.5d0*(unor-abs(unor)) 
    
         !------------------------------------------- 
         ! Prescribed boundary values
         !------------------------------------------- 
         if (bound(b)%Face_ID >= 7) then  ! Object e.g. Hull
            gneg = 0.0d0
            gi   = umi
         else             ! Domain Boundaries
            call getWave(gi, wave_phi, ixint, atime)
            gi   = gi + umi
            gneg = phi - wave_phi
         end if
      
         !------------------------------------------- 
         ! Select appropriate formulation
         !-------------------------------------------
         ! Do nothing BC - zero traction BC 	
         !-------------------------------------------
         if (BCtype(bound(b)%Face_ID) == 0) then
      
            ! Akkerman fix
            do aa = 1, NSHL     
               Rhsu(:,aa) = Rhsu(:,aa) + &
                            shlu(aa)*uneg*rho*(ui-gi)*DetJb*gw(igauss)
            end do

            do bb = 1, NSHL    
               do aa = 1, NSHL 
                  xKebe11(1,aa,bb) = xKebe11(1,aa,bb) - &
                            fact2*shlu(aa)*uneg*rho*shlu(bb)*DetJb*gw(igauss)  
                  xKebe11(5,aa,bb) = xKebe11(5,aa,bb) - &
                            fact2*shlu(aa)*uneg*rho*shlu(bb)*DetJb*gw(igauss)   
                  xKebe11(9,aa,bb) = xKebe11(9,aa,bb) - &
                            fact2*shlu(aa)*uneg*rho*shlu(bb)*DetJb*gw(igauss)
               end do
            end do

            !------------------------------------------- 
            ! Weak dirichlet BC     
            !-------------------------------------------
         else if (BCtype(bound(b)%Face_ID) == 2) then

            call e3bSTAB_weak(tauB, tauNor, ui-umi, nor, dxidx)

            call e3bRHS_weak(nshl, nor, tauB, tauNor, gw(igauss), &
                             shlu, shgradgu, ui, umi, pri, duidxi, gi,&
                             Rhsu, Rhsp)

            call e3bLHS_weak(nshl, ui, umi, duidxi, tauB, tauNor, gw(igauss), &
                             shlu, shgradgu,  &
                             xKebe11, xGebe, xDebe1, nor)
         end if    

         end do

         call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsp)
         call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe) 
         call LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp)    
         
	  
      end do
   end do
         
   if (NSHL.ge.0) then
      deallocate(shlu, shgradgu, shhessgu)
      deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, plold, wl)
      deallocate(xKebe11, xMebe, xGebe, xDebe1, xDebe2)
      deallocate(Rhsu, Rhsp)
   endif

end subroutine FaceAssembly_NS


!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS_conv(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                                acgmAlpha, pgAlpha, phigAlpha)
  
   use aAdjKeep
   use mpi
   use commonvars

   implicit none
   
   real(8), intent(in) :: dgAlpha(NNODE,NSD), ugAlpha(NNODE,NSD), &
                          acgAlpha(NNODE,NSD), ugmAlpha(NNODE,NSD), &
                          acgmAlpha(NNODE,NSD), pgAlpha(NNODE), &
                          phigAlpha(NNODE)

   ! Local variables
   integer :: b,ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl
   
   real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)
   real(8) :: gwt,rDetJb,fact1,fact2

   real(8) :: dxidx(NSD,NSD), dphidx(NSD) 
  
   real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:)
     
   real(8), allocatable :: dl (:,:), ul(:,:), acl(:,:), uml(:,:), &
                           acml(:,:), pl(:),  dlold(:,:), xl(:,:),  &
                           phil(:), plold(:), wl(:)
			 
   real(8), allocatable :: xKebe11(:,:,:), xMebe(:,:), xGebe(:,:,:), &
                           xDebe1(:,:,:), xDebe2(:,:,:), xLSebe(:,:), &
                           xLSUebe(:,:,:), xPLSebe(:,:), xULSebe(:,:,:)
   real(8), allocatable :: Rhsu(:,:), Rhsp(:), Rhsls(:)   
!   real(8), allocatable :: gp(:,:), gw(:), mgp(:,:) 
     
   real(8) :: ui(NSD),pri,duidxi(NSD,NSD), phi, dphidxi(NSD), &
              Gij(NSD,NSD), Ginv(NSD,NSD)
   real(8) :: nor(NSD), bui(NSD),umi(NSD)
   
   real(8) :: tauB, tauNor, gi(NSD), unor, upos, uneg  
   real(8) :: ixint(NSD), bxint(NSD), atime  
   real(8) :: gneg, wave_u, wave_phi
    
   ! shb will be the shape function array while shbg will hold the
   ! gradients of the shape functions  

   
   fact1 = almi
   fact2 = alfi*gami*Delt

   gi = 0.0d0

   !xDebe2  = 0.0d0

   atime = time -(1.0d0-alfi)*Delt

   call genGPandGWb(gp,gw,NGAUSSb)
   

   NSHL = -1

   ! loop over faces

   do b = 1, NBOUND
      do ifac = 1,bound(b)%NFACE

         iel = bound(b)%F2E(ifac)

         if (NSHL.ne.ELMNSHL(iel)) then
            if (NSHL.ge.0) then
               deallocate(shlu, shgradgu, shhessgu)
               deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, &
	                  plold, wl)
               deallocate(xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
	                  xLSebe, xLSUebe, xPLSebe, xULSebe)
               deallocate(Rhsu,Rhsp ,Rhsls)
            endif

            NSHL = ELMNSHL(iel)
  
            allocate(shlu(NSHL), shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD))
            allocate(dl(NSHL,NSD), ul(NSHL,NSD), acl(NSHL,NSD), &
                     uml(NSHL,NSD), acml(NSHL,NSD), pl(NSHL), dlold(NSHL,NSD), &
                     xl(NSHL,NSD), phil(NSHL), plold(NSHL), wl(NSHL))
			 
            allocate(xKebe11(NSD*NSD,NSHL,NSHL), xMebe(NSHL,NSHL), &
	             xGebe(NSD,NSHL,NSHL), xDebe1(NSD,NSHL,NSHL), &
                     xDebe2(NSD,NSHL,NSHL), xLSebe(NSHL,NSHL), &
                     xLSUebe(NSD,NSHL,NSHL), xPLSebe(NSHL,NSHL), &
                     xULSebe(NSD,NSHL,NSHL))
            allocate(Rhsu(NSD,NSHL), Rhsp(NSHL), Rhsls(NSHL))

            xDebe2 = 0.0d0
 
!            call genGPandGWb(gp,gw,NGAUSSb)
           call genGPMap(NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)
            
         endif	  

         ! Get local solution arrays   	       	     
         do i = 1, NSHL
            j = IEN(iel,i)
            xl   (i,:) = xg(j,:)
            dl   (i,:) = dgAlpha(j,:)
            dlold(i,:) = dgold(j,:)
            wl   (i)   = wg(j)
            ul   (i,:) = ugAlpha(j,:)
            acl  (i,:) = acgAlpha(j,:)
            uml  (i,:) = ugmAlpha(j,:)
            acml (i,:) = acgmAlpha(j,:)
            pl   (i)   = pgAlpha(j)
            plold(i)   = pgold(j)
            phil (i)   = phigAlpha(j)
         end do
     
         xKebe11 = 0d0    ! initialize local resistance matrix	     
         xGebe   = 0d0
         xDebe1  = 0d0
         xMebe   = 0d0
         xLSebe  = 0d0
         xLSUebe = 0d0
         xPLSebe = 0d0
         xULSebe = 0d0
    
         Rhsu    = 0d0    ! initialize local load vector
         Rhsp    = 0d0
         Rhsls   = 0d0
	  
         ! Loop over integration points (NGAUSSb in each direction) 
         do igauss = 1, NGAUSSb
           
            ! Evaluate shapes
            call eval_faceshape(nshl, iel, gp(igauss,:), mgp(igauss,:), &
                                bound(b)%FACE_OR(ifac), &
                                xl,dl,wl, &
                                shlu,shgradgu,dxidx,Gij,Ginv,nor)

          
            ! Interpolate	    
            pri = sum(pl*shlu)
            phi = sum(phil*shlu)
	     
            do i = 1, NSD
               ui(i) = sum(ul(:,i)*shlu)
	       umi(i) = sum(uml(:,i)*shlu)
               ixint(i) = sum((xl(:,i)+dl(:,i))*shlu)
            end do     
	    
            do j = 1, NSD
              dphidxi(j) = sum(phil*shgradgu(:,j))
              do i = 1, NSD
                 duidxi(i,j) = sum(ul(:,i)*shgradgu(:,j))
              end do
           end do
         
           ! Secondary quantities      
           call compRhoMu2(phi,dphidxi,dxidx,Ginv)
	    
           unor = sum((ui-umi)*nor(:))  
           upos = 5d-1*(unor+abs(unor))
           uneg = 5d-1*(unor-abs(unor)) 
       
           ! Convection
           if (bound(b)%Face_ID.ge.7) then ! Hull
              gneg = 0d0
           else
              call getWave(wave_u, wave_phi, ixint,atime)
              gneg = phi - wave_phi
           end if
	  
           do aa = 1, NSHL
	      RHSls(aa) = RHSls(aa) &
                          + uneg*gneg*shlu(aa)*DetJb*gw(igauss)	
           end do
		  
           do bb = 1, NSHL
              do aa = 1, NSHL
	         xLSebe(aa,bb) = xLSebe(aa,bb) &
       	               - fact2*uneg*shlu(aa)*shlu(bb)*DetJb*gw(igauss)	  
              end do
           end do  
    
!    do bb = 1, NSHL
!      do aa = 1, NSHL
!       xLSUebe(:,aa,bb) = xLSUebe(:,aa,bb) 
!     &  	   + fact2*max(-1d0,uneg)*nor*shlu(aa)*
!     &        gneg*shlu(aa)*DetJb*gw(igauss) 	 
!  	    	enddo
!	    enddo     
         
        ! Select appropriate formulation
         if (BCtype(bound(b)%Face_ID).eq.0) then  ! Do nothing BC - zero traction BC 	  
            do aa = 1, NSHL     
               Rhsu(:,aa) = Rhsu(:,aa) + &
               shlu(aa)*uneg*rho*ui*DetJb*gw(igauss)
            end do
	  	
            do bb = 1, NSHL    
               do aa = 1, NSHL 
	          xKebe11(1,aa,bb) = xKebe11(1,aa,bb) -    &
     	                     fact2* shlu(aa)*uneg*rho*shlu(bb) *DetJb*gw(igauss)  
	  	  xKebe11(5,aa,bb) = xKebe11(5,aa,bb) -    &
     	                     fact2* shlu(aa)*uneg*rho*shlu(bb) *DetJb*gw(igauss)   
                  xKebe11(9,aa,bb) = xKebe11(9,aa,bb) -   & 
                             fact2* shlu(aa)*uneg*rho*shlu(bb) *DetJb*gw(igauss) 	  
               end do
            end do

         else if (BCtype(bound(b)%Face_ID).eq.2) then  ! WeakBC
     
            call e3bSTAB_weak(tauB, tauNor, ui-umi, nor, dxidx)
  
            ! Integrate
            if (bound(b)%Face_ID.ge.7) then ! Hull
	       gi = umi
            else
               call getWave(gi,phi,ixint,atime)
               gi = gi + umi
            end if
	  
            call e3bRHS_weak(nshl, nor, tauB, tauNor, gw(igauss), &
		             shlu, shgradgu, ui, umi, pri, duidxi, gi,  & 
                             Rhsu, Rhsp)

            call e3bLHS_weak(nshl, ui,umi, duidxi, tauB, tauNor, gw(igauss),&
       		             shlu, shgradgu, &
                             xKebe11, xGebe, xDebe1, nor)
      
            do bb = 1, NSHL
               do aa = 1, NSHL
                  xULSebe(:,aa,bb) = xULSebe(:,aa,bb) -&
       	                    fact2*drhodphi*shlu(aa)*uneg*shlu(bb)*(ui-gi)&
                            *DetJb*gwt 	 
  	       end do
            end do     
        	    	    
            endif    

         end do
        
         ! Assemble in global vector/matrix after application of BCs        
         call BCLhs_NS_conv(nshl, iel, xKebe11, xGebe, xDebe1, &
                            xDebe2, xMebe, xLSebe,xLSUebe,xPLSebe,xULSebe,&
                            Rhsu, Rhsp,RHSls)
	   
         call FillSparseMat_NS_conv(nshl, iel, xKebe11,&
                                    xGebe, xDebe1, xDebe2, xMebe, &
                                    xLSebe,xLSUebe,xPLSebe,xULSebe)
	  
         call LocaltoGlobal_NS_conv(nshl, iel, Rhsu, Rhsp, Rhsls)    
    	
      end do
   end do
   
   if (NSHL.ge.0) then
      deallocate(shlu, shgradgu, shhessgu)
      deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, &
                 plold, wl)
      deallocate(xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                 xLSebe, xLSUebe, xPLSebe, xULSebe)
      deallocate(Rhsu,Rhsp ,Rhsls)
   endif

end subroutine FaceAssembly_NS_conv



!======================================================================
!
!======================================================================
subroutine ForceAssembly_3D(dgAlpha, ugAlpha, ugmAlpha,&
                            acgAlpha, acgmAlpha, pgAlpha, phigAlpha,&
                            gforce, gmom, gforceJ, gmomJ,alphaRB)
  
   use aAdjKeep
   use mpi
   use commonvars

   implicit none
  
   integer :: nshl
       
   ! Arguments	   
   real(8) :: dgAlpha(NNODE,NSD), ugAlpha(NNODE,NSD), &
              acgAlpha(NNODE,NSD), ugmAlpha(NNODE,NSD), &
              acgmAlpha(NNODE,NSD), pgAlpha(NNODE), &
              phigAlpha(NNODE),alphaRB

   real(8) :: gforce(NSD), gmom(NSD),gforceJ(NSD,NSD), gmomJ(NSD,NSD)

   ! Local variables
   integer :: b,ifac, iel, igauss, i, j, k, hess_flag, aa, bb 
  
   real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)
   real(8) :: gwt,rDetJb,fact1,fact2

   real(8) :: dxidx(NSD,NSD), dphidx(NSD)
  
   real(8), allocatable :: shlu(:), shgradgu(:,:),shhessgu(:,:,:)
     
   real(8), allocatable :: dl(:,:), ul(:,:), acl(:,:), uml(:,:), &
                           acml(:,:), pl(:), dlold(:,:), xl(:,:), &
                           phil(:), plold(:), wl(:)
			 
   real(8), allocatable :: xKebe11(:,:,:), xKebe22(:,:,:), xMebe(:,:), &
                           xGebe(:,:,:), xDebe1(:,:,:), xDebe2(:,:,:) 
  
   real(8), allocatable :: Rhsu(:,:), Rhsp(:)     
     
   real(8) :: ui(NSD),pri,duidxi(NSD,NSD), phi, dphidxi(NSD), &
              Gij(NSD,NSD), Ginv(NSD,NSD)
   real(8) :: nor(NSD), bui(NSD),umi(NSD)
  
   real(8) :: tauB, tauNor, gi(NSD),h,He0,rho0,dpriddF3, xm(NSD)    
   real(8) :: ixint(NSD), bxint(NSD),di(NSD),xi(NSD), rot(NSD)
   real(8) :: eforce(NSD), emom(NSD),eforceJ(NSD,NSD), emomJ(NSD,NSD)
       
   fact1 = almi
   fact2 = alfi*gami*Delt

   call genGPandGWb(gp,gw,NGAUSSb)
   
   NSHL    = -1 
   gforce  = 0d0
   gmom    = 0d0
   gforceJ = 0d0
   gmomJ   = 0d0

   ! loop over faces
   do b = 1, NBOUND
      do ifac = 1,bound(b)%NFACE

         iel = bound(b)%F2E(ifac)

         if (bound(b)%Face_ID.ge.7) then

!            call genGPMap(NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)
           
            if (NSHL .ne. ELMNSHL(iel)) then
               if (NSHL.ge.0) then
                  deallocate(shlu, shgradgu, shhessgu)
                  deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, &
                             plold, wl)
		  deallocate(xKebe11, xKebe22, xMebe, xGebe, xDebe1, xDebe2)
                  deallocate(Rhsu, Rhsp)
               end if
  
               NSHL = ELMNSHL(iel)
  		
               allocate(shlu(NSHL), shgradgu(NSHL,NSD),shhessgu(NSHL,NSD,NSD))
               allocate(dl(NSHL,NSD), ul(NSHL,NSD), acl(NSHL,NSD), uml(NSHL,NSD), &
                        acml(NSHL,NSD), pl(NSHL), dlold(NSHL,NSD), xl(NSHL,NSD), &
                        phil(NSHL),  plold(NSHL),wl(NSHL))
			 
               allocate(xKebe11(NSD*NSD,NSHL,NSHL), xKebe22(NSD*NSD,NSHL,NSHL), &
                        xMebe(NSHL,NSHL), xGebe(NSD,NSHL,NSHL), &
                        xDebe1(NSD,NSHL,NSHL),xDebe2(NSD,NSHL,NSHL))
               allocate(Rhsu(NSD,NSHL), Rhsp(NSHL))

           
            end if 
    	  
            ! Get local solution arrays   	   
            do i = 1, NSHL
               j = IEN(iel,i)
               xl(i,:)    = xg(j,:)
               dl(i,:)    = dgAlpha(j,:)
               dlold(i,:) = dgold(j,:)
               wl(i)  = wg(j)
               ul(i,:)    = ugAlpha(j,:)
               acl(i,:)   = acgAlpha(j,:)
               uml(i,:)   = ugmAlpha(j,:)
               acml(i,:)  = acgmAlpha(j,:)
               pl(i)  = pgAlpha(j)
               plold(i)   = pgold(j)
               phil(i)    = phigAlpha(j)
            end do

            eforce = 0d0
            emom   = 0d0
            eforceJ = 0d0
            emomJ   = 0d0
	   
            ! Loop over integration points
            do igauss = 1, NGAUSSb
      
               call eval_faceshape(nshl, iel,gp(igauss,:),mgp(igauss,:), &
                                   bound(b)%FACE_OR(ifac), &
                                   xl,dl,wl, &
                                   shlu,shgradgu,dxidx,Gij,Ginv,nor)
          
               ! Interpolate	    
               pri = sum(pl*shlu)   
               phi = sum(phil*shlu)
  	     
               do i = 1, NSD
                  ui(i)  = sum(ul(:,i)*shlu)
                  umi(i) = sum(uml(:,i)*shlu)  
                  xi(i)  = sum(xl(:,i)*shlu) 
                  di(i)  = sum(dl(:,i)*shlu)   
               enddo     

               do j = 1, NSD
                  dphidxi(j) = sum(phil*shgradgu(:,j))
                  do i = 1, NSD
                     duidxi(i,j) = sum(ul(:,i)*shgradgu(:,j))
                  enddo
               enddo
		  
               ! Derived stuff
               call compRhoMu(nshl, phi,dphidxi,dxidx,Ginv)
               call e3bSTAB_weak(tauB, tauNor, ui-umi, nor, dxidx)
	    	    
               ! Correct pressure for hydrostatic component
               dphidxi(1) = 0d0
               dphidxi(2) = 0d0
               dphidxi(3) = 1d0  
               call getElemSize(h,dxidx,dphidxi,Ginv)  
               call getHeps(He0,water_level-xi(3)-di(3), h)
               rho0 = (1d0-He0)*rhoa + He0*rhow
         
               pri = pri + rho0*gravity*(water_level-xi(3)-di(3))

               xm = (xi+di) - xcg - alphaRB*(dbn1+dbn0)  
     
               ! Integrate
               gi=umi

	       call e3bRHS_3D_force(nshl,nor, tauB, tauNor, gw(igauss), &
       	                            shlu, shgradgu, ui,umi, pri, duidxi, gi,&
                                    xm, eforce, emom)
            end do      

            gmom   = gmom   + emom
            gforce = gforce + eforce
		
         end if
      end do
   end do

   if (NSHL.ge.0) then
      deallocate(shlu, shgradgu, shhessgu)
      deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, &
                 plold, wl)
      deallocate(xKebe11, xKebe22, xMebe, xGebe, xDebe1, xDebe2)
      deallocate(Rhsu, Rhsp)
   end if

    
   if (numnodes.gt.1) then
      eforce   = gforce 
      call MPI_ALLREDUCE(eforce, gforce, NSD, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,  &
                         MPI_COMM_WORLD,mpi_err)
       
      emom   = gmom  
      call MPI_ALLREDUCE(emom, gmom, NSD, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,  &
                         MPI_COMM_WORLD,mpi_err)
   end if    

end subroutine ForceAssembly_3D
