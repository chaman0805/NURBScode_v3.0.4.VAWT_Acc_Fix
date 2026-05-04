!======================================================================
!
!======================================================================
subroutine compRhoMu(phi,dphidxi,dxidx,Ginv)
     
  use commonvars
  implicit none

  real(8), intent(in) :: Ginv(NSD,NSD)
  real(8) phi,dphidxi(NSD), dxidx(NSD,NSD)
  real(8) h, He  
  
  call getElemSize(h,dxidx,dphidxi,Ginv)
  call getHeps(He,phi,h)
  
  rho = (1.0d0-He)*rhoa + He*rhow
  mu  = (1.0d0-He)*mua  + He*muw     

end subroutine compRhoMu
  


!======================================================================
!
!======================================================================    
subroutine compRhoMu2(phi,dphidxi,dxidx,Ginv)
     
  use commonvars
  implicit none

  real(8), intent(in) :: Ginv(NSD,NSD)
  real(8) phi,dphidxi(NSD), dxidx(NSD,NSD)
  real(8) h, He, Hep  

  call getElemSize(h,dxidx,dphidxi,Ginv)
  call getHeps2(He, Hep, phi, h) 
  
  rho = (1d0-He)*rhoa + He*rhow
  mu  = (1d0-He)*mua  + He*muw  
  
  drhodphi = (rhow-rhoa)*Hep 
  dmudphi  = (muw -mua )*Hep 	   
	    
end subroutine compRhoMu2
    


!======================================================================
!
!======================================================================
subroutine getElemSize(h,dxidx,dphidxi,Ginv)
     
  use commonvars
  implicit none
  
  real(8), intent(in) :: Ginv(NSD,NSD)
  real(8) h, dxidx(NSD,NSD), dphidxi(NSD)
  real(8) temp1

  integer i,j
	   	    
  temp1 = 0d+0 ! h in the direction of the gradient

  do i = 1, NSD
     do j = 1, NSD
    temp1 = temp1 + dphidxi(i)*Ginv(i,j)*dphidxi(j)
     enddo
  enddo

  h = sqrt( temp1 /(sum(dphidxi*dphidxi)  +1d-10))
  
end subroutine getElemSize
  


!======================================================================
!
!======================================================================
subroutine getHeps(He, phi, h) 
  use commonpars
  use commonvars
  implicit none
  
  real(8) He,phi,h, eps
    	  
  eps = mp_eps*h

  if (phi.lt.-eps) then
    He  = 0d0
  else if (phi.gt.eps) then
    He  = 1d0
  else
    He  = 5d-1*(1d0 + (phi/eps) + (sin(pi*phi/eps)/pi))
  endif   
  
end subroutine getHeps



!======================================================================
!
!======================================================================
subroutine getHeps2(He, Hep, phi, h) 
  use commonpars
  use commonvars
  implicit none
  
  real(8) He,Hep,phi,h, eps
  eps = mp_eps*h

  if (phi.lt.-eps) then
    He  = 0d0
    Hep = 0d0
  else if (phi.gt.eps) then
    He  = 1d0  
    Hep = 0d0   
  else
    He  = 5d-1*(1d0 + phi/eps + sin(pi*phi/eps)/pi)
    Hep = 5d-1*(	  1d0/eps + cos(pi*phi/eps)/eps)
  endif
  
end subroutine getHeps2



!======================================================================
!
!======================================================================
subroutine getHpH0(HpH0,phi,h) 
  use commonpars
  use commonvars
  implicit none
  
  real(8) HpH0,phi,h, eps
  eps = mp_eps*h

  if (phi.lt.-eps) then
    HpH0 = 0d0
  else if (phi.gt.eps) then  
    HpH0 = 0d0   
  else
    HpH0 = 5d-1*(1d0 + cos(pi*phi/eps))
  endif
    
end subroutine getHpH0
