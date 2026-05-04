!======================================================================
!
!======================================================================
subroutine IntElmAss_convdist(phigAlpha ,rphigAlpha, ugAlpha, &
                              ugmAlpha, dgAlpha)     
  use aAdjKeep
  use commonvars
  implicit none  

  real(8), intent(in) :: phigAlpha (NNODE),rphigAlpha (NNODE), &
                         ugAlpha(NNODE,NSD),ugmAlpha(NNODE,NSD), &
                         dgAlpha(NNODE,NSD)
    
  real(8), allocatable :: shlu(:), shgradgu(:,:), &
                          shconvggu(:), shhessgu(:,:,:)
  real(8), allocatable :: phil(:), rphil(:), ul(:,:), dl(:,:), &
                          xl(:,:), uml(:,:), wl(:), philold (:)
  real(8), allocatable :: xMebe(:,:), rhs(:)  
  real(8), allocatable :: gp(:,:), gw(:) 
       
  integer :: iel, igauss,hess_flag, i, j, k,idx, aa,bb, NGAUSS, nshl  
  real(8) :: ui(NSD),umi(NSD),di(NSD), phi, rphi, dphidxi (NSD)
  real(8) :: dxidx(NSD,NSD)   
  real(8) :: dphiolddxi(NSD),uadvi(NSD)      
  real(8) :: fact1,fact2,tau,kdc(NSD,NSD),res           
  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)


  ! Factors
  fact1 = almi
  fact2 = alfi*gami*Delt 

  ! Get Gaussian points and weights     
  NGAUSS = -1
  NSHL   = -1
    
  ! Loop over elements
  do iel = 1, NELEM
  
    if (NSHL /= ELMNSHL(iel)) then   
		 
      if (NSHL >= 0) then	
        deallocate(shlu, shgradgu, shconvggu, shhessgu, &
                   phil ,rphil, ul, dl, xl, uml, wl, &
                   philold, xMebe, rhs, gp, gw)
      end if	     

      NSHL   = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
      allocate(shlu(NSHL), shgradgu(NSHL,NSD), shconvggu(NSHL), &
               shhessgu(NSHL,NSD,NSD),&   
               phil (NSHL),rphil (NSHL),ul(NSHL,NSD),  &
               dl(NSHL,NSD), xl(NSHL,NSD),uml(NSHL,NSD),wl(NSHL),&
               philold(NSHL), xMebe(NSHL,NSHL),rhs(NSHL),&
               gp(NGAUSS,NSD), gw(NGAUSS))

      call genGPandGW(gp,gw,NGAUSS)
    end if
 
    xMebe = 0d0	     
    rhs   = 0d0
	
    ! Get local solution vectors	    
    do i = 1, NSHL
      idx = IEN(iel,i)	   
      xl(i,:)   = xg(idx,:)
      dl(i,:)   = dgAlpha(idx,:)
      uml(i,:)  = ugmAlpha(idx,:)	
      phil(i)   = phigAlpha (idx)
      rphil(i)  = rphigAlpha(idx)   
      ul(i,:)   = ugAlpha   (idx,:)  	   
      philold(i)= phigold   (idx)
      wl(i)   = wg(idx)
    end do
	
    ! Loop over integration points 
    do igauss = 1, NGAUSS
     
      ! Get Element Shape functions and their gradients
      shlu   = 0d0   
      shgradgu = 0d0
      shhessgu = 0d0
      hess_flag = 0
      call eval_shape(nshl, iel,gp(igauss,:),xl,dl,wl, &
                      shlu,shgradgu,shhessgu,dxidx,Gij,Ginv, &
                      hess_flag)
     	     
      ! Interpolate	          
      phi  = sum(phil *shlu)
      rphi = sum(rphil*shlu)
      do i = 1, NSD   
        ui(i)     = sum(ul(:,i)*shlu)
        umi(i)    = sum(uml(:,i)*shlu)
        di(i)     = sum(dl(:,i)*shlu)
        dphidxi(i)  = sum(phil   *shgradgu(:,i))	
        dphiolddxi(i) = sum(philold*shgradgu(:,i))
      end do  
	    		
      ! Compute tau and kdc
      uadvi = ui - umi	   
	    		 
      do aa = 1, NSHL
        shconvggu(aa) = sum(shgradgu(aa,:)*uadvi)
      end do  

      res = rphi + sum(uadvi*dphidxi)

      call e3STAB_TAU(uadvi, Gij, tau) 			   
      call e3DC_CAU(dphiolddxi,uadvi,Gij,res,tau,kdc)
      kdc = LSC_kdc*kdc  
	   
      ! Calculate residual
      do aa = 1, NSHL
        rhs(aa) = rhs(aa) - (  &
               (shlu(aa)+shconvggu(aa)*tau)*res &
               +shgradgu(aa,1)*kdc(1,1)*dphidxi(1) &
               +shgradgu(aa,1)*kdc(1,2)*dphidxi(2) &
               +shgradgu(aa,1)*kdc(1,3)*dphidxi(3) &
               +shgradgu(aa,2)*kdc(2,1)*dphidxi(1) &
               +shgradgu(aa,2)*kdc(2,2)*dphidxi(2) &
               +shgradgu(aa,2)*kdc(2,3)*dphidxi(3) &
               +shgradgu(aa,3)*kdc(3,1)*dphidxi(1) &
               +shgradgu(aa,3)*kdc(3,2)*dphidxi(2) &
               +shgradgu(aa,3)*kdc(3,3)*dphidxi(3) &  
               )*DetJ*gw(igauss)
      end do

      ! Calculate Jacobian
      do bb = 1, NSHL
        do aa = 1, NSHL
          xMebe(aa,bb) = xMebe(aa,bb) +	(   &
     	             (shlu(aa)+shconvggu(aa)*tau)*(fact1*shlu(bb) &
       	        + fact2*shconvggu(bb)) &
                + fact2*(   &
                +shgradgu(aa,1)*kdc(1,1)*shgradgu(bb,1) &
                +shgradgu(aa,1)*kdc(1,2)*shgradgu(bb,2) &
                +shgradgu(aa,1)*kdc(1,3)*shgradgu(bb,3) &
                +shgradgu(aa,2)*kdc(2,1)*shgradgu(bb,1) &
                +shgradgu(aa,2)*kdc(2,2)*shgradgu(bb,2) &
                +shgradgu(aa,2)*kdc(2,3)*shgradgu(bb,3) &
                +shgradgu(aa,3)*kdc(3,1)*shgradgu(bb,1) &
                +shgradgu(aa,3)*kdc(3,2)*shgradgu(bb,2) &
                +shgradgu(aa,3)*kdc(3,3)*shgradgu(bb,3) & 
      	        ))*DetJ*gw(igauss)  
        end do
      end do 
        	    
    end do
    
    ! Assemble into the Sparse Global Stiffness Matrix and Rhs Vector
    call BCLhs_ls(nshl, iel, xMebe, RHS)   
    call LocaltoGlobal_ls(nshl, iel, rhs)
    call FillSparseMat_ls(nshl, iel, xMebe)  

  end do
  
  deallocate(shlu, shgradgu, shconvggu, shhessgu)
  deallocate(phil, rphil, ul, dl, xl,uml, wl)
  deallocate(philold, xMebe, rhs, gp, gw)		  

end subroutine IntElmAss_convdist



!======================================================================
!
!======================================================================
subroutine FaceAssembly_convdist(phigAlpha, ugAlpha,dgAlpha,ugmAlpha)   
  use aAdjKeep
  use mpi
  use commonvars
  implicit none
          
  real(8), intent(in) :: dgAlpha(NNODE,NSD), ugAlpha(NNODE,NSD), &
                         phigAlpha(NNODE), ugmAlpha(NNODE,NSD)

  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl
  real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:)
   
  real(8), allocatable :: dl(:,:), ul(:,:), acl(:,:), uml(:,:), &
               acml(:,:), pl(:), dlold(:,:), xl(:,:),  &
               phil(:), plold(:), wl(:)  
  real(8), allocatable :: xKebe11(:,:,:), xKebe22(:,:,:), &
               xMebe(:,:), xGebe(:,:,:), xDebe1(:,:,:), &
               xDebe2(:,:,:), Rhsls(:)   
	  
  real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)
  real(8) :: rDetJb,fact1,fact2

  real(8) :: dxidx(NSD,NSD),  dphidx(NSD)
  real(8) :: ui(NSD),pri,duidxi(NSD,NSD), phi, dphidxi(NSD)
  real(8) :: nor(NSD), bui(NSD),umi(NSD)
  real(8) :: tauB, gi(NSD), unor,upos,uneg    
  real(8) :: ixint(NSD), bxint(NSD), gpos,gneg,atime
    
  real(8) :: wave_u(NSD), wave_phi
  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD) 
    
  ! shb will be the shape function array while shbg will hold the
  ! gradients of the shape functions
  fact1 = almi
  fact2 = alfi*gami*Delt 

  atime = time -(1d0-alfi)*Delt

  ! Loop over Faces

  call genGPandGWb(gp,gw,NGAUSSb)

  do b = 1, NBOUND
    do ifac = 1, bound(b)%NFACE
 
      call genGPMap(NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)

      ! get Gauss Point/Weight Arrays
      iel = bound(b)%F2E(ifac)
		 
      NSHL = ELMNSHL(iel)
      allocate(shlu(NSHL), shgradgu(NSHL,NSD),shhessgu(NSHL,NSD,NSD))
      allocate(dl(NSHL,NSD), ul(NSHL,NSD), acl(NSHL,NSD), &
           uml(NSHL,NSD), acml(NSHL,NSD), pl(NSHL), &
          dlold(NSHL,NSD), xl(NSHL,NSD), phil(NSHL), &
          plold(NSHL),wl(NSHL))  
      allocate(xKebe11(NSD*NSD,NSHL,NSHL), xKebe22(NSD*NSD,NSHL,NSHL), &
          xMebe(NSHL,NSHL),xGebe(NSD,NSHL,NSHL), &
          xDebe1(NSD,NSHL,NSHL),xDebe2(NSD,NSHL,NSHL))
      allocate(Rhsls(NSHL))
 	 
      ! Get local solution arrays		 
      do i = 1, NSHL
        j = IEN(iel,i)
        xl(i,:)  = xg(j,:)
        dl(i,:)  = dgAlpha(j,:)
        dlold(i,:) = dgold(j,:)
        wl(i)    = wg(j)
        ul(i,:)  = ugAlpha(j,:)
        uml(i,:)   = ugmAlpha(j,:)
        phil(i)  = phigAlpha(j)	 
      end do
     
      ! initialize local stiffness matrix  	     
      xMebe = 0d+0		 
      Rhsls = 0d+0	    ! initialize local load vector
           
      ! Loop over integration points
      do igauss = 1, NGAUSSb
     
        ! Evaluate shapes   
        call eval_faceshape(nshl, iel,gp(igauss,:),mgp(igauss,:), &
                bound(b)%FACE_OR(ifac), &
                xl,dl,wl, &
                shlu,shgradgu,dxidx,Gij,Ginv,nor)
      
        ! Interpolate	  
        phi = sum(phil*shlu)
   
        do i = 1, NSD
          ui(i) = sum(ul(:,i)*shlu)
	  umi(i) = sum(uml(:,i)*shlu)	    
          ixint(i) = sum((xl(:,i)+dl(:,i))*shlu)
        end do
	  	  
        ! Get normal velocity                  
        unor = sum((ui-umi)*nor)
        upos = 5d-1*(unor + abs(unor))
        uneg = 5d-1*(unor - abs(unor))
		 		  
        call getWave(wave_u, wave_phi, ixint,atime)

        ! BC value
        if (bound(b)%Face_ID.ge.7) then ! Hull
	  gneg = 0d0
        else
          call getWave(wave_u, wave_phi, ixint,atime)
          gneg = phi - wave_phi
        end if
	  
        ! Residual vector      
	do aa = 1, NSHL
	  RHSls(aa) = RHSls(aa)  &
    	       + uneg*gneg*shlu(aa)*DetJb*gw(igauss)	
        end do

        ! Tangent matrix		  
        do bb = 1, NSHL
          do aa = 1, NSHL
	    xMebe(aa,bb) = xMebe(aa,bb) &
    	          - fact2*uneg*shlu(aa)*shlu(bb)*DetJb*gw(igauss)	    
          end do
        end do
    
      end do
      
      ! assemble into the Sparse Global Stiffness Matrix and Rhs Vector        
      call BCLhs_ls(nshl, iel, xMebe, Rhsls)
      call FillSparseMat_ls(nshl, iel, xMebe) 
      call LocaltoGlobal_ls(nshl, iel, Rhsls) 
    
      deallocate(shlu, shgradgu,shhessgu)
      deallocate(dl, ul, acl, uml, acml, pl, dlold, xl, phil, plold, wl)  
      deallocate(xKebe11, xKebe22, xMebe, xGebe, xDebe1,xDebe2)
      deallocate(Rhsls)
		
    end do
  end do
      
end subroutine FaceAssembly_convdist 
