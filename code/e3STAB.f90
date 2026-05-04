!======================================================================
! Stabilization parameters (tauM, tauC) for VMS
!====================================================================== 
subroutine e3STAB_3D(Gij, Ginv, AD_VEL_L, rLi, tauM, tauP, &
                     tauC, tauBar, uprime, cfl)   
  use aAdjKeep
  use commonvars
  implicit none
      
  real(8), intent(in)  :: Gij(NSD,NSD), Ginv(NSD,NSD), &
                          AD_VEL_L(NSD), rLi(NSD)
  real(8), intent(out) :: tauM, tauC, tauP, tauBar, &
                          cfl(2), uprime(NSD)

  real(8) :: taua, taut, taud, m_k, dtfact, nu, taus, gij2
  integer :: i, j

  dtfact = 4.0d0            
!!!  dtfact = 16.0d0
!!!  dtfact = 64.0d0
!!!  dtfact = 256.0d0

  m_k = 3.0d0
!!!  m_k    = 36.0d0!*real(max(Pu**2,Qu**2,Ru**2))

  cfl = 0.0d0         
!  ! get cfl(2)
!  taua  = 0.0d0
!  do j = 1, NSD
!    do i = 1, NSD
!      taua = taua + AD_VEL_L(i)*Ginv(i,j)*AD_VEL_L(j)
!    end do
!  end do   
!  cfl(2)  = Delt*sum(AD_VEL_L*AD_VEL_L)/sqrt(taua+1d-8)    

  ! get tauM, tauC, cfl(1)
  taua = 0.0d0
  taut = 0.0d0
  taud = 0.0d0
  gij2 = 0.0d0
  do j = 1, NSD
    do i = 1, NSD
      taua = taua + AD_VEL_L(i)*Gij(i,j)*AD_VEL_L(j)
      gij2 = gij2 + Gij(i,j)*Gij(i,j)
    end do
  end do

  cfl(1) = Delt*sqrt(taua)     

  nu   = mu/rho
  taud = nu*nu*m_k*gij2
  taut = dtfact*Dtgl*Dtgl

  taus = taua + taud
      
  ! Get tauM - so far a scalar
  tauM = 1.0d0/(rho*sqrt(taus+taut))
  tauP = tauM

  ! Get tauC
  tauC = rho*sqrt(taus)/(Gij(1,1)+Gij(2,2)+Gij(3,3))
!!!  tauC = rho*sqrt(taus)/sqrt(gij2)
!!!  tauC = 1.0d0/(tauM*(Gij(1,1)+Gij(2,2)+Gij(3,3)))

  ! Tayfun's tauC
!!!  tauC = rho**2*tauM*sum(AD_VEL_L*AD_VEL_L)

  uprime(:) = -tauM*rLi(:) 

  ! Get TauBar
  tauBar = 0.0d0
  do i = 1, NSD
    do j = 1, NSD
      tauBar = tauBar + uprime(i)*Gij(i,j)*uprime(j)
    end do
  end do
  tauBar = 1.0d0/sqrt(tauBar+1.0d-15)

end subroutine e3STAB_3D

  

!======================================================================
!
!======================================================================
subroutine e3DC_beta2(ui, duidxi, Gij, res, dxidx, kdc)

  use aAdjKeep
  use commonvars
  
  implicit none
  
  real(8), intent(in) :: Gij(NSD,NSD)
  real(8) :: ui(NSD),duidxi(NSD,NSD),res(NSD),He,dxidx(NSD,NSD) 
  integer :: i,j
  real(8) :: h4,resnorm,kdc

  He  = 1d0!!!(rho - rhoa)/(rhow - rhoa)

  ! Convection based meshsize   
  h4 = 0d0
  do j = 1, NSD
    do i = 1, NSD
      h4  = h4 + Gij(i,j)*Gij(i,j)
    end do
  end do
  
  ! Residual norm
  resnorm = sqrt(sum(res(:)*res(:)))
       
  ! KDC parameter  
  kdc = resnorm*(He*NS_kdc_w  + (1d0-He)*NS_kdc_a)/sqrt(h4)

end subroutine e3DC_beta2 
  


!======================================================================
!
!======================================================================
subroutine e3DC_shakibNS(ui,duidxi,Ginv,res,dxidx, KAPPA_DC)

  use aAdjKeep
  use commonvars  
  implicit none
  
  real(8), intent(in) :: Ginv(NSD,NSD)
  real(8) :: ui(NSD), duidxi(NSD,NSD), res(NSD), He, dxidx(NSD,NSD) 
  real(8) :: KAPPA_DC(NSD,NSD)
  integer :: i,j
  real(8) :: resnorm, temp(NSD)

  He  = (rho - rhoa)/(rhow - rhoa)

  ! Grad norm in reference coordinates    
  temp = 0.0d0
  do j = 1, NSD
    do i = 1, NSD
      temp  = temp + duidxi(:,i)*Ginv(i,j)*duidxi(:,j)
    end do
  end do
  
  ! Residual norm
  resnorm = sqrt(sum(res(:)*res(:)))    

  ! Kdc  
  KAPPA_DC = Ginv*(He*NS_kdc_w  + (1d0-He)*NS_kdc_a)*&
             resnorm/(sqrt(sum(temp)) + 1d-3)

end subroutine e3DC_shakibNS


!======================================================================
! Stabilization parameter (tauB) for weak BC
!======================================================================
subroutine e3bSTAB_weak(tauB, tauNor, ui, nor, dxidx)

  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  real(8), intent(in)  :: ui(NSD), nor(NSD), dxidx(NSD,NSD)
  real(8), intent(out) :: tauB, tauNor
  
  real(8) :: temp(NSD), hn, unor, upos, uneg  
  integer :: i

  do i = 1, NSD
    temp(i) = sum(dxidx(i,:)*nor)
  end do

!!!  hn = 1.0d0/(sqrt(sum(temp*temp)))  ! hn for tets
  hn = 2.0d0/(sqrt(sum(temp*temp)))  ! hn for other elements

  tauB = 4.0d0*mu/hn
  
  unor = sum(ui*nor(:))  
  upos = 0.5d0*(unor+abs(unor))
  uneg = 0.5d0*(unor-abs(unor))
  
  tauNor = tauB !+ rho*upos

end subroutine e3bSTAB_weak


!======================================================================
! Stabilization parameter (tauB) for DG
!======================================================================
subroutine e3bSTAB_DG(tauB, ui, nor, dxidx)
  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  real(8), intent(in)  :: ui(NSD), nor(NSD), dxidx(NSD,NSD)
  real(8), intent(out) :: tauB
  
  real(8) :: temp(NSD), hn, unor, uneg
  integer :: i

  do i = 1, NSD
    temp(i) = sum(dxidx(i,:)*nor)
  end do

  hn   = 1.0d0/(sqrt(sum(temp*temp)))  
  tauB = 16.0d0*mu/hn

  unor = sum(ui*nor(:))

  uneg = 0.5d0*(unor-abs(unor))
  tauB = tauB - rho*uneg

!!!  tauB = tauB + rho*abs(unor)

end subroutine e3bSTAB_DG
