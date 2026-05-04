!======================================================================
!
!======================================================================
subroutine e3Rhs_3D_fluid(nshl, ui, aci, umi, acmi, uadvi, pri, &
                          rLi, fi, duidxi, ddidxi, tauM, tauP,  &
                          tauC, tauBar, kap_dc, gwt, shgu,      &
                          shgradgu, uprime, Rhsu, Rhsp)
  use aAdjKeep
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl      
  real(8), intent(in) :: ui(NSD), aci(NSD), umi(NSD), acmi(NSD),  &
                         uadvi(NSD), pri, rLi(NSD), fi(NSD),      &
                         duidxi(NSD,NSD), ddidxi(NSD,NSD), gwt,   &
                         shgu(NSHL), shgradgu(NSHL,NSD), kap_dc,  &
                         tauM, tauP, tauC, tauBar, uprime(NSD)

  real(8), intent(inout) :: Rhsu(NSD,NSHL), Rhsp(NSHL)

  integer :: aa, bb, i, j, k
  real(8) :: fact1, fact2, mupkdc, divu, ptot, advu(NSD), &
             tmp1(NSD), tmp2(NSD,NSD), tmp4(NSD)

  tmp1 = 0.0d0; tmp2 = 0.0d0; tmp4 = 0.0d0; divu = 0.0d0

!!!  uprime(:) = -tauM*rLi(:)
      
  mupkdc = mu + kap_dc

  divu = duidxi(1,1) + duidxi(2,2) + duidxi(3,3)

  advu(:) = uadvi(:) + uprime(:)
          
  ptot = pri - tauC*divu

  tmp1(:) = rho*(aci(:) + advu(1)*duidxi(:,1) &
                        + advu(2)*duidxi(:,2) &
                        + advu(3)*duidxi(:,3)) - fi(:)
            
  tmp4(:) = uprime(1)*duidxi(:,1) + &
            uprime(2)*duidxi(:,2) + &
            uprime(3)*duidxi(:,3)
          
  tmp2(1,1) = -ptot + 2.0d0*mupkdc*duidxi(1,1)   &
            - rho*(uadvi(1)+uprime(1))*uprime(1) &
            + rho*uprime(1)*tauBar*tmp4(1)
      
  tmp2(1,2) = mupkdc*(duidxi(1,2) + duidxi(2,1)) &
            - rho*(uadvi(2)+uprime(2))*uprime(1) &
            + rho*uprime(2)*tauBar*tmp4(1)
      
  tmp2(1,3) = mupkdc*(duidxi(1,3) + duidxi(3,1)) &
            - rho*(uadvi(3)+uprime(3))*uprime(1) &
            + rho*uprime(3)*tauBar*tmp4(1)
      
  tmp2(2,1) = mupkdc*(duidxi(2,1) + duidxi(1,2)) &
            - rho*(uadvi(1)+uprime(1))*uprime(2) &
            + rho*uprime(1)*tauBar*tmp4(2) 
      
  tmp2(2,2) = -ptot + 2.0d0*mupkdc*duidxi(2,2)   &
            - rho*(uadvi(2)+uprime(2))*uprime(2) &
            + rho*uprime(2)*tauBar*tmp4(2)
      
  tmp2(2,3) = mupkdc*(duidxi(2,3) + duidxi(3,2)) &
            - rho*(uadvi(3)+uprime(3))*uprime(2) &
            + rho*uprime(3)*tauBar*tmp4(2)
      
  tmp2(3,1) = mupkdc*(duidxi(3,1) + duidxi(1,3)) &
            - rho*(uadvi(1)+uprime(1))*uprime(3) &
            + rho*uprime(1)*tauBar*tmp4(3)
     
  tmp2(3,2) = mupkdc*(duidxi(3,2) + duidxi(2,3)) &
            - rho*(uadvi(2)+uprime(2))*uprime(3) &
            + rho*uprime(2)*tauBar*tmp4(3)
      
  tmp2(3,3) = -ptot + 2.0d0*mupkdc*duidxi(3,3)   &
            - rho*(uadvi(3)+uprime(3))*uprime(3) &
            + rho*uprime(3)*tauBar*tmp4(3)
      
  ! Physics Residual   
  do aa = 1, NSHL     
    do i = 1, NSD
      Rhsu(i,aa) = Rhsu(i,aa) - &
                   ( shgu(aa)*tmp1(i) + &
                     sum(shgradgu(aa,:)*tmp2(i,:)) )*DetJ*gwt
    end do
  end do
      
  do aa = 1, NSHL
    Rhsp(aa) = Rhsp(aa) - &
               ( shgu(aa)*divu - &
                 sum(shgradgu(aa,:)*uprime(:)) )*DetJ*gwt
  end do

!!!  do aa = 1, NSHL
!!!    Rhsp(aa) = Rhsp(aa) - &
!!!               ( shgu(aa)*divu + &
!!!                 sum(shgradgu(aa,:)*tauP*rLi(:)) )*DetJ*gwt
!!!  end do

end subroutine e3Rhs_3D_fluid


!======================================================================
!
!======================================================================
subroutine e3Rhs_3D_struct(nshl, aci, ui, fi, Ftens, Stens, gwt, &
                           shlu, shgradgu, Rhsu, dc)
  
  use aAdjKeep
  use commonvars

  implicit none
  
  integer, intent(in) :: nshl
  
  integer :: aa, i, j
  
  real(8) :: gwt, shlu(NSHL), shgradgu(NSHL,NSD), &
             fact1, aci(NSD), ui(NSD), fi(NSD), &
             Rhsu(NSD,NSHL), Ftens(NSD,NSD), Stens(NSD,NSD), dc
  
  real(8) :: tmp1(NSD), tmp2(NSD,NSD)
  
  tmp1 = 0d+0
  tmp2 = 0d+0
  
!  dC = 5.0d+4
!  dc = 6.0d+5
!  dC = 0.0d0

  tmp1(:) = rho*aci(:) + dC*ui(:) - fi(:)
  
  ! First P-K Stress  
  do j = 1, NSD
     do i = 1, NSD
       tmp2(i,j) = Ftens(i,1)*Stens(1,j) &
                 + Ftens(i,2)*Stens(2,j) &
                 + Ftens(i,3)*Stens(3,j) 
!                + FPKS0(i,j)
     end do
  end do
  
  ! Discrete Residual  
  do aa = 1, NSHL
     
     Rhsu(1,aa) = Rhsu(1,aa) - &
         (shlu(aa)*tmp1(1) + &
         shgradgu(aa,1)*tmp2(1,1) + &
         shgradgu(aa,2)*tmp2(1,2) + &
         shgradgu(aa,3)*tmp2(1,3))*DetJ*gwt
     
     Rhsu(2,aa) = Rhsu(2,aa) - &
         (shlu(aa)*tmp1(2) + &
         shgradgu(aa,1)*tmp2(2,1) + &
         shgradgu(aa,2)*tmp2(2,2) + &
         shgradgu(aa,3)*tmp2(2,3))*DetJ*gwt
     
     Rhsu(3,aa) = Rhsu(3,aa) - &
         (shlu(aa)*tmp1(3) + &
         shgradgu(aa,1)*tmp2(3,1) + &
         shgradgu(aa,2)*tmp2(3,2) + &
         shgradgu(aa,3)*tmp2(3,3))*DetJ*gwt
     
  end do

end subroutine e3Rhs_3D_struct



!======================================================================
!
!======================================================================
subroutine e3Rhs_3D_mesh(nshl, ddidxi, gwt, shlu, shgradgu, Rhsm)
  
  use aAdjKeep
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl
  
  integer :: aa, bb, i, j, k

  real(8) :: shlu(NSHL), shgradgu(NSHL,NSD)
  real(8) :: ui(NSD), aci(NSD), umi(NSD), acmi(NSD), uadvi(NSD),&
             pri, rLi(NSD), fi(NSD), duidxi(NSD,NSD), ddidxi(NSD,NSD),&
             uprime(NSD), advu(NSD), strn((NSD+1)*NSD/2)
  real(8) :: tauM, tauP, tauC, tauBar, gwt, gwt0
  real(8) :: Rhsu(NSD,NSHL), Rhsm(NSD,NSHL), Rhsp(NSHL) 
  real(8) :: fact1, fact2, &
             tmp1(NSD), tmp2(NSD,NSD), tmp3(NSD,NSD), &
             divu, tmp4(NSD), mbulk, mshear, vval, numod, Emod
    
  ! Poisson Ratio and Young's modulus
  ! for mesh motion
  numod = 3d-1
  Emod  = 1d0 

  ! Mesh "elastic" parameters  
  mbulk = numod*Emod/((1d+0+numod)*(1d+0-2d+0*numod))
  mshear = Emod/(2d+0*(1d+0+numod)) ! PDE
  
  tmp1 = 0d+0
  tmp2 = 0d+0
  
  strn(1) = ddidxi(1,1)
  strn(2) = ddidxi(2,2)
  strn(3) = ddidxi(3,3)
  
  strn(4) = ddidxi(1,2) + ddidxi(2,1)
  strn(5) = ddidxi(2,3) + ddidxi(3,2)
  strn(6) = ddidxi(1,3) + ddidxi(3,1)
    
  tmp1 = 0d+0
  tmp2 = 0d+0
  
  tmp1(:) = acmi(:)
  
  tmp2(1,1) = (mbulk+2d+0*mshear)*strn(1) + &
               mbulk*strn(2) + mbulk*strn(3)
  
  tmp2(1,2) = mshear*strn(4)
  
  tmp2(1,3) = mshear*strn(6)
  
  tmp2(2,1) = tmp2(1,2)
  
  tmp2(2,2) = (mbulk+2d+0*mshear)*strn(2) + &
               mbulk*strn(1) + mbulk*strn(3)
  
  tmp2(2,3) = mshear*strn(5)
  
  tmp2(3,1) = tmp2(1,3)
  
  tmp2(3,2) = tmp2(2,3)
  
  tmp2(3,3) = (mbulk+2d+0*mshear)*strn(3) + &
               mbulk*strn(1) + mbulk*strn(2)  
    
  do aa = 1, NSHL
     
    Rhsm(1,aa) = Rhsm(1,aa) - &
         (shgradgu(aa,1)*tmp2(1,1) + &
          shgradgu(aa,2)*tmp2(1,2) + &
          shgradgu(aa,3)*tmp2(1,3))*gwt
    
    Rhsm(2,aa) = Rhsm(2,aa) - &
         (shgradgu(aa,1)*tmp2(2,1) + &
          shgradgu(aa,2)*tmp2(2,2) + &
          shgradgu(aa,3)*tmp2(2,3))*gwt
     
    Rhsm(3,aa) = Rhsm(3,aa) - &
         (shgradgu(aa,1)*tmp2(3,1) + &
          shgradgu(aa,2)*tmp2(3,2) + &
          shgradgu(aa,3)*tmp2(3,3))*gwt     
  end do
    
end subroutine e3Rhs_3D_mesh



!======================================================================
! RHS for weak BC
!======================================================================
subroutine e3bRHS_weak(nshl, nor, tauB, tauNor, gwt, &
                       shlu, shgradgu, ui, umi, pri, duidxi, gi, &
                       Rhsu, Rhsp, ti, tmp1, tmp2)  
  use aAdjKeep  
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl

  real(8), intent(inout) :: Rhsu(NSD,NSHL), Rhsp(NSHL)
  real(8), intent(out)   :: ti(NSD), tmp1(NSD), tmp2(NSD,NSD)
  real(8), intent(in)    :: shlu(NSHL), shgradgu(NSHL,NSD), &
                            nor(NSD), tauB, tauNor, gwt, &
                            ui(NSD), umi(NSD), pri, duidxi(NSD,NSD), &
                            gi(NSD)

  integer :: aa, bb, i  
  real(8) :: fact1, fact2, upos, uneg, unor, tr, pt33, gmul, gnor
  
  tmp1 = 0.0d0
  tmp2 = 0.0d0
  
  ti(:) = -pri*nor(:) + mu*(duidxi(:,1)*nor(1) + &
                            duidxi(:,2)*nor(2) + &
                            duidxi(:,3)*nor(3) + &
                            duidxi(1,:)*nor(1) + &
                            duidxi(2,:)*nor(2) + &
                            duidxi(3,:)*nor(3))
  
  ! Relative normal velocity for convective term
  unor = sum((ui-umi)*nor)  ! u \cdot n
  upos = 0.5d0*(unor+abs(unor))
  uneg = 0.5d0*(unor-abs(unor))  

  ! Absolute normal for the rest
  gnor = sum(gi*nor)
  unor = sum(ui*nor)  ! u \cdot n
 
  tmp1(:) = - ti(:) &
            + tauB*(ui(:)-gi(:)) &           
            - uneg*rho*(ui(:)-gi(:))      &
            + (tauNor-tauB)*unor*nor(:)

  ! gmul =  1.0d0 => sym
  ! gmul = -1.0d0 => skew       
  gmul = 1.0d0 
  do aa = 1, NSD
    do bb = 1, NSD
      tmp2(aa,bb) = -gmul*mu*((ui(aa)-gi(aa))*nor(bb) &
                             +(ui(bb)-gi(bb))*nor(aa))
    end do
  end do

  ! gnor = 0.0d0
  do aa = 1, NSHL
    do i = 1, NSD
      Rhsu(i,aa) = Rhsu(i,aa) - &
                 ( shlu(aa)*tmp1(i) + &
                   sum(shgradgu(aa,:)*tmp2(i,:)) )*DetJb*gwt
    end do

    Rhsp(aa) = Rhsp(aa) + shlu(aa)*(unor-gnor)*DetJb*gwt     
  end do

end subroutine e3bRHS_weak



!======================================================================
!
!======================================================================
subroutine e3bRHS_3D_force(nshl, nor, tauB, tauNor, gwt, shlu, shgradgu, &
                           ui, umi, pri, duidxi, gi, xi, efor, emom)
  
  use aAdjKeep  
  use commonvars  
  
  integer, intent(in) :: nshl

  real(8), intent(inout) :: efor(NSD), emom(NSD)

  real(8), intent(in) :: shlu(NSHL), shgradgu(NSHL,NSD), &
                         nor(NSD), tauB, tauNor, gwt, &
                         ui(NSD), umi(NSD), pri, duidxi(NSD,NSD), &
                         gi(NSD), xi(NSD)

  integer :: aa, bb  
  real(8) :: fact1, fact2, tmp1(NSD), tmp2(NSD,NSD), &
             upos, uneg, unor, tr, pt33, gmul, ti(NSD)
  
  tmp1 = 0d0
  tmp2 = 0d0  
     
  ti(:) = -pri*nor(:) + mu*(duidxi(:,1)*nor(1) + &
                            duidxi(:,2)*nor(2) + &
                            duidxi(:,3)*nor(3) + &
                            duidxi(1,:)*nor(1) + &
                            duidxi(2,:)*nor(2) + &
                            duidxi(3,:)*nor(3))
  
  unor = sum((ui-umi)*nor(:))  ! u \cdot n
  upos = 0.5d0*(unor+abs(unor))
  uneg = 0.5d0*(unor-abs(unor))  
 
  tmp1(:) = - ti(:) &
            + tauB*(ui(:)-gi(:)) &
            - uneg*rho*(ui(:)-gi(:))    &  
            + (tauNor-tauB)*unor*nor(:)

  ! gmul =  1.0d0 => sym
  ! gmul = -1.0d0 => skew       
  gmul = 1.0d0
  do aa = 1, NSD
    do bb = 1, NSD
      tmp2(aa,bb) = -gmul*mu*((ui(aa)-gi(aa))*nor(bb) &
                             +(ui(bb)-gi(bb))*nor(aa))
    end do
  end do
     
  efor(1) = efor(1) + tmp1(1)*DetJb*gwt     
  efor(2) = efor(2) + tmp1(2)*DetJb*gwt     
  efor(3) = efor(3) + tmp1(3)*DetJb*gwt

  emom(1) = emom(1) + (xi(2)*tmp1(3) - xi(3)*tmp1(2) + &
                       tmp2(3,2) - tmp2(2,3))*DetJb*gwt
 
  emom(2) = emom(2) + (xi(3)*tmp1(1) - xi(1)*tmp1(3) + &
                       tmp2(1,3) - tmp2(3,1))*DetJb*gwt

  emom(3) = emom(3) + (xi(1)*tmp1(2) - xi(2)*tmp1(1) + &
                       tmp2(2,1) - tmp2(1,2))*DetJb*gwt

end subroutine e3bRHS_3D_force



!======================================================================
! RHS for DG/Non-matching
! for w1t1, fact=1.0; for w2t2, fact=-1.0
!======================================================================
subroutine e3bRHS_DG(nshl,shgu, shgradgu, gwt, ui, ti, nor, tauB, &
                     Rhsu, Rhsp, sgn)
  use aAdjKeep  
  use commonvars
  implicit none
  
  integer, intent(in)    :: nshl
  real(8), intent(inout) :: Rhsu(NSD,NSHL), Rhsp(NSHL)
  real(8), intent(in)    :: shgu(NSHL), shgradgu(NSHL,NSD), gwt, &
                            ui(NSD), ti(NSD),&
                            nor(NSD), tauB, sgn
  integer :: aa, bb, i
  real(8) :: tmp1(NSD), tmp2(NSD,NSD), upos, uneg, unor, gmul
  
  tmp1 = 0.0d0
  tmp2 = 0.0d0

  unor = sum(ui*nor)
 
  tmp1(:) = -0.5d0*ti(:) + tauB*ui(:)

  ! viscous contribution (in weighting functions)
  do aa = 1, NSD
    do bb = 1, NSD
      tmp2(aa,bb) = -0.5d0*mu*(ui(aa)*nor(bb)+ui(bb)*nor(aa))
    end do
  end do

  do aa = 1, NSHL
    do i = 1, NSD
      Rhsu(i,aa) = Rhsu(i,aa) - sgn*&
                 ( shgu(aa)*tmp1(i) + &
                   sum(shgradgu(aa,:)*tmp2(i,:)) )*DetJb*gwt
    end do
    Rhsp(aa) = Rhsp(aa) + sgn*0.5d0*shgu(aa)*unor*DetJb*gwt
  end do
end subroutine e3bRHS_DG
