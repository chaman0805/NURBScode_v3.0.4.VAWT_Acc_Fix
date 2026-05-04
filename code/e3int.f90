!======================================================================
!
!======================================================================
subroutine e3int_fluid(nshl, xl, dl, ul, acl, uml, acml, pl, fl, phil, &
                       shlu, shgradlu, shgradgu, shhessgu, dxidx,      &
                       Ginv, di, ui, aci, umi, acmi, pri, fi, ddidxi,  &
                       duidxi, duidxixj, dpridxi, phi, dphidxi, xi, rLi)
  use aAdjKeep
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl  
  real(8), intent(in) :: Ginv(NSD,NSD)

  integer :: i,j,k
  real(8) :: shlu(NSHL), shgradlu(NSHL,NSD), &
             shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD),dxidx(NSD,NSD)
  real(8) :: dl(NSHL,NSD), ul(NSHL,NSD),xl(NSHL,NSD), &
             acl(NSHL,NSD), uml(NSHL,NSD), acml(NSHL,NSD), pl(NSHL), &
             fl(NSHL,NSD),phil(NSHL)
  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), pri, &
             fi(NSD), ddidxi(NSD,NSD), &
             duidxi(NSD,NSD), duidxixj(NSD,NSD,NSD), &
             dpridxi(NSD),phi,dphidxi(NSD),  rLi(NSD), xi(NSD)
  real(8) :: uadv(NSD),He0,rho0, h,dphidxi0(NSD)
 
  xi   = 0.0d0
  di   = 0.0d0
  ui   = 0.0d0
  aci  = 0.0d0
  umi  = 0.0d0
  acmi = 0.0d0
  pri  = 0.0d0
  fi   = 0.0d0
  phi  = 0.0d0
  
  ddidxi   = 0.0d0
  duidxi   = 0.0d0
  duidxixj = 0.0d0
  dpridxi  = 0.0d0
  dphidxi  = 0.0d0
  
  rLi = 0.0d0 
       
  ! Get quantities and their grads and hessians  
  do i = 1, NSHL
    xi(:) = xi(:) + xl(i,:)*shlu(i)
    di(:) = di(:) + dl(i,:)*shlu(i)
    ui(:) = ui(:) + ul(i,:)*shlu(i)
    aci(:) = aci(:) + acl(i,:)*shlu(i)
    umi(:) = umi(:) + uml(i,:)*shlu(i)
    acmi(:) = acmi(:) + acml(i,:)*shlu(i)
    pri = pri + pl(i)*shlu(i)
    fi(:) = fi(:) + fl(i,:)*shlu(i)
    phi = phi + phil(i)*shlu(i)
    dpridxi(:) = dpridxi(:) + pl(i)*shgradgu(i,:)
    dphidxi(:) = dphidxi(:) + phil(i)*shgradgu(i,:)
    do j = 1, NSD
      ddidxi(:,j) = ddidxi(:,j)+dl(i,:)*shgradgu(i,j)
      duidxi(:,j) = duidxi(:,j)+ul(i,:)*shgradgu(i,j)
      do k = 1, NSD
        duidxixj(:,j,k) = duidxixj(:,j,k)+ul(i,:)*shhessgu(i,j,k)
      end do
    end do
  end do
    
!  call compRhoMu(phi, dphidxi, dxidx, Ginv)
  rho = rhow
  mu  = muw

!  ! NS-ALE PDE Residual
!  dphidxi0(1) = 0.0d0
!  dphidxi0(2) = 0.0d0
!  dphidxi0(3) = 1.0d0  
!  call getElemSize(h,dxidx,dphidxi0,Ginv)  
!  call getHeps(He0,water_level-xi(3)-di(3), h)
!  rho0 = (1.0d0-He0)*rhoa + He0*rhow   
  
!  fi = fi + (rho-rho0)*gravvec 
  
  uadv(:) = ui(:) - umi(:)
       
  rLi(:) = rho*aci(:) + &
           rho*(duidxi(:,1)*uadv(1) + &
                duidxi(:,2)*uadv(2) + &
                duidxi(:,3)*uadv(3)) + &
           dpridxi(:) - fi(:) - &
           mu*(duidxixj(:,1,1)+duidxixj(:,2,2)+duidxixj(:,3,3))

end subroutine e3int_fluid
  

!======================================================================
!
!======================================================================
subroutine e3int_struct(nshl, dl, ul, acl, fl, shlu, shgradgu, &
                        di, ui, aci, fi, ddidxi)
  use aAdjKeep
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl
  
  integer :: i,j,k
  real(8) :: shlu(NSHL), shgradgu(NSHL,NSD)
  real(8) :: dl(NSHL,NSD), ul(NSHL,NSD),  &
             acl(NSHL,NSD), fl(NSHL,NSD)
  real(8) :: di(NSD), ui(NSD), aci(NSD), fi(NSD), ddidxi(NSD,NSD)
  
  
  di  = 0.0d0
  ui  = 0.0d0
  aci = 0.0d0
  fi  = 0.0d0
  ddidxi = 0.0d0
  
  ! Get quantities and their grads and hessians
  do i = 1, NSHL
    di(:) = di(:) + dl(i,:)*shlu(i)
    ui(:) = ui(:) + ul(i,:)*shlu(i)
    aci(:) = aci(:) + acl(i,:)*shlu(i)
    fi(:) = fi(:) + fl(i,:)*shlu(i)
    do j = 1, NSD
      ddidxi(:,j) = ddidxi(:,j)+dl(i,:)*shgradgu(i,j)
    enddo
  enddo  

end subroutine e3int_struct

  
  

!======================================================================
!
!======================================================================
subroutine e3tensors_nonlin(NSD,ddidxi,mu,lambda,Ftens,Stens,Ctens)

  implicit none
  
  integer :: i,j,k,l,NSD
  real(8) :: ddidxi(NSD,NSD)
  real(8) :: Ftens(NSD,NSD), Stens(NSD,NSD), Ctens(NSD,NSD,NSD,NSD)
  real(8) :: Itens(NSD,NSD), Jdet, Finv(NSD,NSD), tmp, &
             Ct(NSD,NSD), CtInv(NSD,NSD), Jmod, trC, pt33, pt5, &
             fact1, fact2, fact3, mu, lambda, Ctemp(NSD,NSD,NSD,NSD)
    
  pt33 = 1.0d0/3.0d0
  pt5  = 0.5d0
  
  Itens = 0.0d0
  Itens(1,1) = 1.0d0
  Itens(2,2) = 1.0d0
  Itens(3,3) = 1.0d0
  
  Ftens = Itens + ddidxi    ! F = I + du/dX
  
  Jdet = 0.0d0
  Finv = 0.0d0
  tmp  = 0.0d0
         
  ! compute the inverse of deformation gradient and Jacobian  
  Finv(1,1) = Ftens(2,2) * Ftens(3,3) &
            - Ftens(3,2) * Ftens(2,3)
  Finv(1,2) = Ftens(3,2) * Ftens(1,3)  &
            - Ftens(1,2) * Ftens(3,3)
  Finv(1,3) = Ftens(1,2) * Ftens(2,3)  &
            - Ftens(1,3) * Ftens(2,2)
  tmp       = 1.0d0 / ( Finv(1,1) * Ftens(1,1) & 
                      + Finv(1,2) * Ftens(2,1)   &
                      + Finv(1,3) * Ftens(3,1) )
  Finv(1,1) = Finv(1,1) * tmp
  Finv(1,2) = Finv(1,2) * tmp
  Finv(1,3) = Finv(1,3) * tmp
  Finv(2,1) = (Ftens(2,3) * Ftens(3,1)  &
             - Ftens(2,1) * Ftens(3,3)) * tmp
  Finv(2,2) = (Ftens(1,1) * Ftens(3,3)  &
             - Ftens(3,1) * Ftens(1,3)) * tmp
  Finv(2,3) = (Ftens(2,1) * Ftens(1,3)  &
             - Ftens(1,1) * Ftens(2,3)) * tmp
  Finv(3,1) = (Ftens(2,1) * Ftens(3,2)  &
             - Ftens(2,2) * Ftens(3,1)) * tmp
  Finv(3,2) = (Ftens(3,1) * Ftens(1,2)  &
             - Ftens(1,1) * Ftens(3,2)) * tmp
  Finv(3,3) = (Ftens(1,1) * Ftens(2,2)  &
             - Ftens(1,2) * Ftens(2,1)) * tmp
  
  Jdet = 1.0d0/tmp       ! J
  Jmod = Jdet**(-2.0d0/3.0d0) ! J^-2/3
!!!  Jmod = Jdet
  
  Ct = 0.0d0
    
  do i = 1, NSD
    do j = 1, NSD
      Ct(i,j) = sum(Ftens(:,i)*Ftens(:,j))
     end do
  end do
  
     
  ! compute the inverse of the F^T F tensor     
  CtInv = 0.0d0
  tmp = 0.0d0
  trC = 0.0d0
  
  CtInv(1,1) = Ct(2,2) * Ct(3,3)  &
             - Ct(3,2) * Ct(2,3)
  CtInv(1,2) = Ct(3,2) * Ct(1,3)  &
             - Ct(1,2) * Ct(3,3)
  CtInv(1,3) = Ct(1,2) * Ct(2,3)  &
             - Ct(1,3) * Ct(2,2)
  tmp      = 1.0d0 / ( CtInv(1,1) * Ct(1,1)  &
                     + CtInv(1,2) * Ct(2,1)  &
                     + CtInv(1,3) * Ct(3,1) )
  CtInv(1,1) = CtInv(1,1) * tmp
  CtInv(1,2) = CtInv(1,2) * tmp
  CtInv(1,3) = CtInv(1,3) * tmp
  CtInv(2,1) = (Ct(2,3) * Ct(3,1) &
              - Ct(2,1) * Ct(3,3)) * tmp
  CtInv(2,2) = (Ct(1,1) * Ct(3,3) &
              - Ct(3,1) * Ct(1,3)) * tmp
  CtInv(2,3) = (Ct(2,1) * Ct(1,3) &
              - Ct(1,1) * Ct(2,3)) * tmp
  CtInv(3,1) = (Ct(2,1) * Ct(3,2) &
              - Ct(2,2) * Ct(3,1)) * tmp
  CtInv(3,2) = (Ct(3,1) * Ct(1,2) &
              - Ct(1,1) * Ct(3,2)) * tmp
  CtInv(3,3) = (Ct(1,1) * Ct(2,2) &
              - Ct(1,2) * Ct(2,1)) * tmp
    
  trC = Ct(1,1)+Ct(2,2)+Ct(3,3) ! Trace of the F^T Ftensor
  
  ! Compute Second P-K tensor
  ! for Neo-Hookean Mat w/ slight compressibility
  
  Stens = 0.0d0
  
  Stens(:,:) = 0.5d0*lambda*(Jdet*Jdet-1.0d0)*CtInv(:,:) +&
                             mu*Jmod*(Itens(:,:) - pt33*trC*CtInv(:,:)) 

  ! Compute Nonlinear Tangent
  ! Stiffness Tensor C_{IJKL}
  
  Ctens = 0.0d0
  
  fact1 = lambda*Jdet*Jdet+2d+0*pt33*pt33*mu*Jmod*trC
  fact2 = 2d+0*pt33*mu*Jmod
  fact3 = 2d+0*pt33*mu*Jmod*trC-lambda*(Jdet*Jdet-1.0d0)
  fact3 = fact3*0.5d0
  
  do i = 1, NSD
    do j = 1, NSD
      do k = 1, NSD
        do l = 1, NSD
          Ctemp(i,j,k,l) = fact1*CtInv(i,j)*CtInv(k,l)-&
                           fact2*(Itens(i,j)*CtInv(k,l)+&
                                  CtInv(i,j)*Itens(k,l))+&
                           fact3*(CtInv(i,k)*CtInv(j,l)+&
                                  CtInv(i,l)*CtInv(j,k))
        end do
      end do
    end do
  end do

  ! Build C_{ijJI} = F_{iK}C_{KIJL}F_{jL}
  do i = 1, NSD
    do j = 1, NSD 
      do k = 1, NSD
        do l = 1, NSD
          Ctens(i,j,l,k) = Ftens(i,1)*Ftens(j,1)*Ctemp(1,k,l,1) + &
                           Ftens(i,1)*Ftens(j,2)*Ctemp(1,k,l,2) + &
                           Ftens(i,1)*Ftens(j,3)*Ctemp(1,k,l,3) + &
                           Ftens(i,2)*Ftens(j,1)*Ctemp(2,k,l,1) + &
                           Ftens(i,2)*Ftens(j,2)*Ctemp(2,k,l,2) + &
                           Ftens(i,2)*Ftens(j,3)*Ctemp(2,k,l,3) + &
                           Ftens(i,3)*Ftens(j,1)*Ctemp(3,k,l,1) + &
                           Ftens(i,3)*Ftens(j,2)*Ctemp(3,k,l,2) + &
                           Ftens(i,3)*Ftens(j,3)*Ctemp(3,k,l,3)
        end do
      end do
    end do
  end do
  
end subroutine e3tensors_nonlin

  
