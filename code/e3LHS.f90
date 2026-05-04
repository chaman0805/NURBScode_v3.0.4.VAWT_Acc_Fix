!======================================================================
!
!======================================================================
subroutine e3LHS_3D_fluid_Old(nshl, ui, umi, aci, pri, duidxi,       &
                              dpridxi, rLi, tauM, tauP, tauC,        &
                              tauBar, kap_dc, gwt, shgu, shgradgu,   &
                              shhessgu, xKebe11,                     &
                              xGebe, xDebe1, xMebe)    
  use aAdjKeep
  use commonvars
  implicit none  
  
  integer, intent(in) :: nshl
  
  real(8), intent(in) :: ui(NSD), umi(NSD), aci(NSD), pri, kap_dc, &
                         duidxi(NSD,NSD), dpridxi(NSD), rLi(NSD), &
                         shgu(NSHL), shgradgu(NSHL,NSD), &
                         shhessgu(NSHL,NSD,NSD), &
                         gwt, tauM, tauP, tauC, tauBar

  real(8), intent(inout) :: xKebe11(NSD*NSD,NSHL,NSHL), &
                            xGebe(NSD,NSHL,NSHL),&
                            xDebe1(NSD,NSHL,NSHL),& 
                            xMebe(NSHL,NSHL)

  integer :: aa, bb, i, j
 
  real(8) :: fact1, fact2, fact3,&
             tmp1(NSHL), tmp2(NSHL), tmp4(NSHL,NSHL),&
             advu1(NSD), advu2(NSD), mupkdc
  
  ! loop over local shape functions in each direction
  fact1 = almi
  fact2 = alfi*gami*Delt
  fact3 = alfi*beti*Delt*Delt
    
  mupkdc = mu + kap_dc

  tmp1 = 0.0d0
  tmp2 = 0.0d0
  tmp4 = 0.0d0

  advu1(:) = ui(:)-umi(:)
  advu2(:) = -tauM*rLi(:)
    
  tmp1(:) = rho*(advu1(1)*shgradgu(:,1) + &! Na,_j (u_j-um_j)
                 advu1(2)*shgradgu(:,2) + &
                 advu1(3)*shgradgu(:,3))
  
  tmp2(:) = advu2(1)*shgradgu(:,1) + &! Na,_i (-tauM*Li)
            advu2(2)*shgradgu(:,2) + &
            advu2(3)*shgradgu(:,3)
  
  
  do bb = 1, NSHL      ! Diagonal blocks of K11
    do aa = 1, NSHL
      
      tmp4(aa,bb) = fact1*(shgu(aa)*rho*shgu(bb)+&
                           tmp1(aa)*tauM*rho*shgu(bb)) +&
                           fact2*(shgu(aa)*tmp1(bb) +&
                           tmp1(aa)*tauM*tmp1(bb) +&
                           mupkdc*(shgradgu(aa,1)*shgradgu(bb,1) +&
                                   shgradgu(aa,2)*shgradgu(bb,2) +&
                                   shgradgu(aa,3)*shgradgu(bb,3)) +&
                           tmp2(aa)*rho*tauBar*tmp2(bb))
    end do
  end do
 
  ! Physics-Physics Interaction
  do bb = 1, NSHL
    do aa = 1, NSHL      
      xKebe11(1,aa,bb) = xKebe11(1,aa,bb) +&
        (tmp4(aa,bb) +&
         fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,1) +&
                shgradgu(aa,1)*tauC*shgradgu(bb,1)))*DetJ*gwt
      
      xKebe11(5,aa,bb) = xKebe11(5,aa,bb) +&
        (tmp4(aa,bb) +&
         fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,2) +&
                shgradgu(aa,2)*tauC*shgradgu(bb,2)))*DetJ*gwt
      
      xKebe11(9,aa,bb) = xKebe11(9,aa,bb) +&
        (tmp4(aa,bb) +&
         fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,3) +&
                shgradgu(aa,3)*tauC*shgradgu(bb,3)))*DetJ*gwt
      
      xKebe11(2,aa,bb) = xKebe11(2,aa,bb) + &
        fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,1) +&
               shgradgu(aa,1)*tauC*shgradgu(bb,2))*DetJ*gwt
      
      xKebe11(4,aa,bb) = xKebe11(4,aa,bb) + &
        fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,2) +&
               shgradgu(aa,2)*tauC*shgradgu(bb,1))*DetJ*gwt
      
      xKebe11(3,aa,bb) = xKebe11(3,aa,bb) + &
        fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,1) +&
               shgradgu(aa,1)*tauC*shgradgu(bb,3))*DetJ*gwt
      
      xKebe11(7,aa,bb) = xKebe11(7,aa,bb) + &
        fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,3) +&
               shgradgu(aa,3)*tauC*shgradgu(bb,1))*DetJ*gwt
      
      xKebe11(6,aa,bb) = xKebe11(6,aa,bb) + &
        fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,2) +&
               shgradgu(aa,2)*tauC*shgradgu(bb,3))*DetJ*gwt
      
      xKebe11(8,aa,bb) = xKebe11(8,aa,bb) + &
        fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,3) +&
               shgradgu(aa,3)*tauC*shgradgu(bb,2))*DetJ*gwt
    end do    
  end do

  ! Physics-Mesh
  ! Divergence Matrix  
  do bb = 1, NSHL   
    do aa = 1, NSHL      
      xGebe(1,aa,bb) = xGebe(1,aa,bb) + &
        fact2*(-shgradgu(aa,1)*shgu(bb) +&
               tmp1(aa)*tauM*shgradgu(bb,1))*DetJ*gwt
      xGebe(2,aa,bb) = xGebe(2,aa,bb) + &
        fact2*(-shgradgu(aa,2)*shgu(bb) +&
               tmp1(aa)*tauM*shgradgu(bb,2))*DetJ*gwt
      xGebe(3,aa,bb) = xGebe(3,aa,bb) + &
        fact2*(-shgradgu(aa,3)*shgu(bb) +&
               tmp1(aa)*tauM*shgradgu(bb,3))*DetJ*gwt
    end do    
  end do

  ! Physics-Physics
  ! Divergence Matrix  
  do bb = 1, NSHL
    do aa = 1, NSHL      
      xDebe1(1,aa,bb) = xDebe1(1,aa,bb) +&
        (fact2*(shgu(aa)*shgradgu(bb,1)+&
                shgradgu(aa,1)*tauP*tmp1(bb)) +&
         fact1*(shgradgu(aa,1)*tauP*rho*shgu(bb)))*DetJ*gwt

      xDebe1(2,aa,bb) = xDebe1(2,aa,bb) +&
        (fact2*(shgu(aa)*shgradgu(bb,2)+&
                shgradgu(aa,2)*tauP*tmp1(bb)) +&
         fact1*(shgradgu(aa,2)*tauP*rho*shgu(bb)))*DetJ*gwt
      
      xDebe1(3,aa,bb) = xDebe1(3,aa,bb) +&
        (fact2*(shgu(aa)*shgradgu(bb,3)+&
                shgradgu(aa,3)*tauP*tmp1(bb)) +&
         fact1*(shgradgu(aa,3)*tauP*rho*shgu(bb)))*DetJ*gwt   
    end do    
  end do   

  ! Mass Matrix
  do bb = 1, NSHL
    do aa = 1, NSHL
      xMebe(aa,bb) = xMebe(aa,bb) +&
        fact2*tauP*(shgradgu(aa,1)*shgradgu(bb,1) +&
                    shgradgu(aa,2)*shgradgu(bb,2) +&
                    shgradgu(aa,3)*shgradgu(bb,3))*DetJ*gwt
    end do
  end do
  
end subroutine e3LHS_3D_fluid_Old


!======================================================================
!
!======================================================================
subroutine e3LHS_3D_fluid(nshl, ui, umi, aci, pri, duidxi, &
                          dpridxi, rLi, tauM,  tauC, kdc, &
                          gwt, shlu, shgradlu, shgradgu, shhessgu, &
                          xKebe11, xGebe, xDebe1, xMebe)
  
  use aAdjKeep
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl
    
  integer aa, bb, i, j, k

  real(8)  gwt, gwt0, tauM, tauP, tauC, tauBar, vval
  
  real(8)  shlu(NSHL), shgradlu(NSHL,NSD), &
   shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD)

  
  real(8) xKebe11(NSD*NSD,NSHL,NSHL), &
   xKebe21(NSD*NSD,NSHL,NSHL), xKebe22(NSD*NSD,NSHL,NSHL), &
   xGebe(NSD,NSHL,NSHL), &
   xDebe1(NSD,NSHL,NSHL),  &
   xMebe(NSHL,NSHL)

  real(8) ui(NSD), umi(NSD), aci(NSD), pri, duidxi(NSD,NSD), &
   dpridxi(NSD), rLi(NSD)
  
  real(8)  fact1, fact2, fact3,kdc,mupkdc, &
       tmp1(NSHL), tmp2(NSHL),tmp3(NSHL), tmp4(NSHL,NSHL), &
       tmp5(NSHL,NSD), &
       advu1(NSD),  advu2(NSD), &
       mbulk, mshear, divu, numod, Emod,DetJgwt
  
  DetJgwt = 1d0*DetJ*gwt !! q-weight
!...  loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt

  mupkdc = mu+kdc
  
  advu1(:) = ui(:)-umi(:)
  advu2(:) = -tauM*rLi(:)

  tmp1(:) = rho*(advu1(1)*shgradgu(:,1) + & ! Na,_j (u_j-um_j)
                 advu1(2)*shgradgu(:,2) + &
                 advu1(3)*shgradgu(:,3))
  
  tmp2(:) = rho*(advu2(1)*shgradgu(:,1) + & ! Na,_i (-tauM*Li)
                 advu2(2)*shgradgu(:,2) + &
                 advu2(3)*shgradgu(:,3))
  
  tmp3(:) = -mu*(shhessgu(:,1,1) &
               + shhessgu(:,2,2) + shhessgu(:,3,3) )
    
  
  do bb = 1, NSHL      ! Diagonal blocks of K11
    do aa = 1, NSHL
      
       tmp4(aa,bb) = fact1*(shlu(aa)*     rho*shlu(bb)+ &
               tmp1(aa)*tauM*rho*shlu(bb)) + &
           fact2*(shlu(aa)*(tmp1(bb) + tmp2(bb)) +     &
          (tmp2(aa)+tmp1(aa))*tauM*(tmp1(bb) + tmp3(bb)) + &
           mupkdc*(shgradgu(aa,1)*shgradgu(bb,1) + &
           shgradgu(aa,2)*shgradgu(bb,2) + &
           shgradgu(aa,3)*shgradgu(bb,3)))
       
    end do
  end do


  do bb = 1, NSHL
    do aa = 1, NSHL
      
       xKebe11(1,aa,bb) = xKebe11(1,aa,bb) + &
       (tmp4(aa,bb) + &
       fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,1) + &
       shgradgu(aa,1)*tauC*shgradgu(bb,1)))*DetJ*gwt
      
      xKebe11(5,aa,bb) = xKebe11(5,aa,bb) + &
       (tmp4(aa,bb) + &
       fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,2) + &
       shgradgu(aa,2)*tauC*shgradgu(bb,2)))*DetJ*gwt
      
      xKebe11(9,aa,bb) = xKebe11(9,aa,bb) + &
       (tmp4(aa,bb) + &
       fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,3) + &
       shgradgu(aa,3)*tauC*shgradgu(bb,3)))*DetJ*gwt
      
      xKebe11(2,aa,bb) = xKebe11(2,aa,bb) +  &
       fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,1) + &
       shgradgu(aa,1)*tauC*shgradgu(bb,2))*DetJ*gwt
      
      xKebe11(4,aa,bb) = xKebe11(4,aa,bb) +  &
       fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,2) + &
       shgradgu(aa,2)*tauC*shgradgu(bb,1))*DetJ*gwt

      
      xKebe11(3,aa,bb) = xKebe11(3,aa,bb) +  &
       fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,1) + &
       shgradgu(aa,1)*tauC*shgradgu(bb,3))*DetJ*gwt
      
      xKebe11(7,aa,bb) = xKebe11(7,aa,bb) +  &
       fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,3) + &
       shgradgu(aa,3)*tauC*shgradgu(bb,1))*DetJ*gwt
      
      xKebe11(6,aa,bb) = xKebe11(6,aa,bb) +  &
       fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,2) + &
       shgradgu(aa,2)*tauC*shgradgu(bb,3))*DetJ*gwt
      
      xKebe11(8,aa,bb) = xKebe11(8,aa,bb) +  &
       fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,3) + &
       shgradgu(aa,3)*tauC*shgradgu(bb,2))*DetJ*gwt
      
    enddo    
  enddo


  do bb = 1, NSHL     
    do aa = 1, NSHL
    
    xGebe(1,aa,bb) = xGebe(1,aa,bb) + &
            fact2*(-shgradgu(aa,1)*shlu(bb) + &
            tmp1(aa)*tauM*shgradgu(bb,1))*DetJ*gwt
    xGebe(2,aa,bb) = xGebe(2,aa,bb) + &
            fact2*(-shgradgu(aa,2)*shlu(bb) +&
            tmp1(aa)*tauM*shgradgu(bb,2))*DetJ*gwt
    xGebe(3,aa,bb) = xGebe(3,aa,bb) + &
            fact2*(-shgradgu(aa,3)*shlu(bb) +&
            tmp1(aa)*tauM*shgradgu(bb,3))*DetJ*gwt
          
    end do     
  end do
  

  do bb = 1, NSHL     
    do aa = 1, NSHL
      
    xDebe1(1,aa,bb) = xDebe1(1,aa,bb) + &
            (fact2*(shlu(aa)*shgradgu(bb,1) + &
            shgradgu(aa,1)*tauM*(tmp1(bb) + tmp3(bb) )) + &
            fact1*(shgradgu(aa,1)*tauM*rho*shlu(bb)))*DetJgwt
    
    xDebe1(2,aa,bb) = xDebe1(2,aa,bb) + &
            (fact2*(shlu(aa)*shgradgu(bb,2) + &
            shgradgu(aa,2)*tauM*(tmp1(bb) + tmp3(bb) )) + &
            fact1*(shgradgu(aa,2)*tauM*rho*shlu(bb)))*DetJgwt
    
    xDebe1(3,aa,bb) = xDebe1(3,aa,bb) + &
            (fact2*(shlu(aa)*shgradgu(bb,3) + &
            shgradgu(aa,3)*tauM*(tmp1(bb) + tmp3(bb) )) + &
            fact1*(shgradgu(aa,3)*tauM*rho*shlu(bb)))*DetJgwt
        
    end do     
  end do


  do bb = 1, NSHL
    do aa = 1, NSHL
      xMebe(aa,bb) = xMebe(aa,bb) + &
            fact2*tauM*(shgradgu(aa,1)*shgradgu(bb,1) + &
                        shgradgu(aa,2)*shgradgu(bb,2) + &
                        shgradgu(aa,3)*shgradgu(bb,3))*DetJgwt
    end do
  end do

end subroutine e3LHS_3D_fluid



!======================================================================
!
!======================================================================
subroutine e3LHS_3D_struct(nshl, gwt, shg, shgradg, Ftens, &
                           Stens, Ctens, xKebe, dc)
  
  use aAdjKeep
  use commonvars

  implicit none
  
  integer, intent(in) :: nshl
    
  integer aa, bb, i,j,k,l,m,o
  
  real(8)  gwt, shg(NSHL), shgradg(NSHL,NSD), &
          xKebe(NSD*NSD,NSHL,NSHL), fact1, fact2, pt5
  
  real(8) Ftens(NSD,NSD), Stens(NSD,NSD), Ctens(NSD,NSD,NSD,NSD)
  
  real(8) temp1(NSHL,NSHL), temp2(NSD,NSD,NSD,NSD), &
          temp3(NSD,NSD,NSD,NSD), temp4(NSHL,NSHL,NSD,NSD), &
          temp5(NSHL,NSHL,NSD,NSD), temp6(NSHL,NSHL)

  real(8) dC
    
  pt5 = 0.5d+0

!  dC = 6.0d+5
!  dC = 5d+4
!  dC = 0.0d0

  fact1 = alfi*beti*Delt*Delt
  fact2 = alfi*gami*Delt

  temp1 = 0d+0
  
  do bb = 1, NSHL
     do aa = 1, NSHL   ! Mass Contribution
    temp1(aa,bb) = almi*shg(aa)*rho*shg(bb) &
                 + fact2*shg(aa)*dC*shg(bb)
     enddo
  enddo
  
! ijlo component
! Build material contribution
  
  temp4 = 0d+0
  
  do j = 1, NSD
    do i = 1, NSD
      do bb = 1, NSHL
        do aa = 1, NSHL
          temp4(aa,bb,i,j) &
              = shgradg(aa,1)*Ctens(i,j,1,1)*shgradg(bb,1) &
              + shgradg(aa,2)*Ctens(i,j,1,2)*shgradg(bb,1) &   
              + shgradg(aa,3)*Ctens(i,j,1,3)*shgradg(bb,1) &
              + shgradg(aa,1)*Ctens(i,j,2,1)*shgradg(bb,2) &
              + shgradg(aa,2)*Ctens(i,j,2,2)*shgradg(bb,2) &   
              + shgradg(aa,3)*Ctens(i,j,2,3)*shgradg(bb,2) &
              + shgradg(aa,1)*Ctens(i,j,3,1)*shgradg(bb,3) &
              + shgradg(aa,2)*Ctens(i,j,3,2)*shgradg(bb,3) &   
              + shgradg(aa,3)*Ctens(i,j,3,3)*shgradg(bb,3)
        end do
      end do
    enddo
  end do
    
  temp4 = temp4*fact1

            ! Build geometric nonlinearity

  temp6 = 0d+0
  
  do bb = 1, NSHL
    do aa = 1, NSHL
      temp6(aa,bb) = shgradg(aa,1)*Stens(1,1)*shgradg(bb,1) &
                   + shgradg(aa,2)*Stens(1,2)*shgradg(bb,1) &
                   + shgradg(aa,3)*Stens(1,3)*shgradg(bb,1) &
                   + shgradg(aa,1)*Stens(2,1)*shgradg(bb,2) &
                   + shgradg(aa,2)*Stens(2,2)*shgradg(bb,2) &
                   + shgradg(aa,3)*Stens(2,3)*shgradg(bb,2) &
                   + shgradg(aa,1)*Stens(3,1)*shgradg(bb,3) &
                   + shgradg(aa,2)*Stens(3,2)*shgradg(bb,3) &
                   + shgradg(aa,3)*Stens(3,3)*shgradg(bb,3)
    end do
  end do
  
  temp6 = temp6*fact1
  
  ! loop over elements in each direction
 
  do bb = 1, NSHL
     do aa = 1, NSHL
    
       xKebe(1,aa,bb) = xKebe(1,aa,bb) + &
            (temp1(aa,bb) +  & ! Mass
             temp4(aa,bb,1,1) +  & ! Mat Stiff
             temp6(aa,bb)  )*DetJ*gwt ! Geom Stiff
    
       xKebe(5,aa,bb) = xKebe(5,aa,bb) + &
             (temp1(aa,bb) + &! Mass
              temp4(aa,bb,2,2) + &! Mat Stiff
              temp6(aa,bb)  )*DetJ*gwt ! Geom Stiff 
    
       xKebe(9,aa,bb) = xKebe(9,aa,bb) + &
            (temp1(aa,bb) + &! Mass
            temp4(aa,bb,3,3) + &! Mat Stiff
            temp6(aa,bb)  )*DetJ*gwt ! Geom Stiff
    
       xKebe(2,aa,bb) = xKebe(2,aa,bb) + temp4(aa,bb,1,2)*DetJ*gwt ! Mat Stiff    
       xKebe(4,aa,bb) = xKebe(4,aa,bb) + temp4(aa,bb,2,1)*DetJ*gwt ! Mat Stiff    
       xKebe(3,aa,bb) = xKebe(3,aa,bb) + temp4(aa,bb,1,3)*DetJ*gwt ! Mat Stiff    
       xKebe(7,aa,bb) = xKebe(7,aa,bb) + temp4(aa,bb,3,1)*DetJ*gwt ! Mat Stiff    
       xKebe(6,aa,bb) = xKebe(6,aa,bb) + temp4(aa,bb,2,3)*DetJ*gwt ! Mat Stiff    
       xKebe(8,aa,bb) = xKebe(8,aa,bb) + temp4(aa,bb,3,2)*DetJ*gwt ! Mat Stiff
    
    end do
  end do
  
  return
  end


!======================================================================
!
!======================================================================
subroutine e3LHS_3D_mesh(nshl, gwt, shgradgu, xKebe22)
  

  use aAdjKeep
  use commonvars

  implicit none
  
  integer, intent(in) :: nshl
    
  integer aa, bb, i, j

  real(8)  gwt
  
  real(8) shgradgu(NSHL,NSD)  
  real(8) xKebe22(NSD*NSD,NSHL,NSHL)
  
  real(8)  fact1, fact2, fact3, &
       tmp1(NSHL), tmp2(NSD), &
       tmp3(NSHL), tmp4(NSHL,NSHL), &
       tmp5(NSHL,NSD), &
       advu1(NSD),  advu2(NSD), &
       mbulk, mshear, divu, numod, Emod
  
  ! loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt
  fact3 = alfi*beti*Delt*Delt

  ! Mesh "elastic" parameters
  numod = 3d-1
  Emod  = 1d0
  mbulk = numod*Emod/((1d+0+numod)*(1d+0-2d+0*numod))
  mshear = Emod/(2d+0*(1d+0+numod)) ! PDE
        
  ! K22 - Mesh-Mesh Interaction - Use Linear Elasticity  
  do bb = 1, NSHL     
     do aa = 1, NSHL
      
      xKebe22(1,aa,bb) = xKebe22(1,aa,bb) +  &
          fact3*((mbulk + 2d+0*mshear)*shgradgu(aa,1)*shgradgu(bb,1) + &  !St
          mshear*shgradgu(aa,2)*shgradgu(bb,2) + &
          mshear*shgradgu(aa,3)*shgradgu(bb,3))*gwt

      
      xKebe22(5,aa,bb) = xKebe22(5,aa,bb) + &
          fact3*(mshear*shgradgu(aa,1)*shgradgu(bb,1) + &! Stiff 
       (mbulk + 2d+0*mshear)*shgradgu(aa,2)*shgradgu(bb,2) +&
        mshear*shgradgu(aa,3)*shgradgu(bb,3))*gwt
      
      xKebe22(9,aa,bb) = xKebe22(9,aa,bb) + &
          fact3*(mshear*shgradgu(aa,1)*shgradgu(bb,1) + &! Stiff
        mshear*shgradgu(aa,2)*shgradgu(bb,2) + &
       (mbulk+2d+0*mshear)*shgradgu(aa,3)*shgradgu(bb,3))*gwt
      
      
      xKebe22(2,aa,bb) = xKebe22(2,aa,bb) +  &
       fact3*(mbulk*shgradgu(aa,1)*shgradgu(bb,2) + &
       mshear*shgradgu(aa,2)*shgradgu(bb,1))*gwt
      
      xKebe22(4,aa,bb) = xKebe22(4,aa,bb) +  &
       fact3*(mshear*shgradgu(aa,1)*shgradgu(bb,2) + &
       mbulk*shgradgu(aa,2)*shgradgu(bb,1))*gwt
      
      xKebe22(3,aa,bb) = xKebe22(3,aa,bb) +  &
       fact3*(mbulk*shgradgu(aa,1)*shgradgu(bb,3) + &
       mshear*shgradgu(aa,3)*shgradgu(bb,1))*gwt
      
      xKebe22(7,aa,bb) = xKebe22(7,aa,bb) +  &
       fact3*(mshear*shgradgu(aa,1)*shgradgu(bb,3) + &
       mbulk*shgradgu(aa,3)*shgradgu(bb,1))*gwt
      
      xKebe22(6,aa,bb) = xKebe22(6,aa,bb) +  &
       fact3*(mbulk*shgradgu(aa,2)*shgradgu(bb,3) + &
       mshear*shgradgu(aa,3)*shgradgu(bb,2))*gwt
      
      xKebe22(8,aa,bb) = xKebe22(8,aa,bb) + &
       fact3*(mshear*shgradgu(aa,2)*shgradgu(bb,3) + &
       mbulk*shgradgu(aa,3)*shgradgu(bb,2))*gwt
      
   enddo

  enddo


end subroutine e3LHS_3D_mesh



!======================================================================
! LHS for weak BC
!======================================================================
subroutine e3bLHS_weak(nshl, ui, umi, duidxi, tauB, tauNor, gwt, &
                       shlu, shgradgu, xKebe, xGebe, xDebe, nor)
  
  use aAdjKeep  
  use commonvars
  implicit none
  
  integer, intent(inout) :: nshl
  
  real(8), intent(in) :: ui(NSD), umi(NSD), duidxi(NSD,NSD), &
                         tauB, tauNor, gwt, nor(NSD), &
                         shlu(NSHL), shgradgu(NSHL,NSD)
  real(8), intent(inout) :: xKebe(NSD*NSD,NSHL,NSHL), &
                            xGebe(NSD,NSHL,NSHL), &
                            xDebe(NSD,NSHL,NSHL)

  integer :: aa, bb
  real(8) :: fact1, fact2, tmp1(NSHL), tmp2(NSHL,NSHL), &
             unor, umul, munor, gmul, uneg 
  
  ! loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt

  tmp1 = 0.0d0
  tmp2 = 0.0d0

  tmp1(:) = shgradgu(:,1)*nor(1) + shgradgu(:,2)*nor(2) &
          + shgradgu(:,3)*nor(3) 

  unor = sum((ui-umi)*nor(:))  ! u \cdot n
  uneg = 0.5d0*(unor-abs(unor))
  munor = tauNor-tauB
  
  ! gmul =  1d0 => sym
  ! gmul = -1d0 => skew  
  gmul = 1.0d0
  do bb = 1, NSHL      ! Diagonal blocks of K
    do aa = 1, NSHL
    
      tmp2(aa,bb) = -shlu(aa)*mu*tmp1(bb) &
                  - gmul*tmp1(aa)*mu*shlu(bb) &
                  + shlu(aa)*tauB*shlu(bb) &
                  - shlu(aa)*uneg*rho*shlu(bb)
    end do
  end do


  do bb = 1, NSHL    
    do aa = 1, NSHL
      
      xKebe(1,aa,bb) = xKebe(1,aa,bb) + &
        fact2*( tmp2(aa,bb)                            &
              - shlu(aa)*mu*shgradgu(bb,1)*nor(1)      &
              - gmul*shgradgu(aa,1)*nor(1)*mu*shlu(bb) &
              + shlu(aa)*nor(1)*munor*nor(1)*shlu(bb) )*DetJb*gwt
       
      xKebe(5,aa,bb) = xKebe(5,aa,bb) + &
        fact2*(  tmp2(aa,bb)                           &
              - shlu(aa)*mu*shgradgu(bb,2)*nor(2)      &
              - gmul*shgradgu(aa,2)*nor(2)*mu*shlu(bb) &
              + shlu(aa)*nor(2)*munor*nor(2)*shlu(bb) )*DetJb*gwt
       
      xKebe(9,aa,bb) = xKebe(9,aa,bb) + &
        fact2*( tmp2(aa,bb)                            &
              - shlu(aa)*mu*shgradgu(bb,3)*nor(3)      &
              - gmul*shgradgu(aa,3)*nor(3)*mu*shlu(bb) &
              + shlu(aa)*nor(3)*munor*nor(3)*shlu(bb) )*DetJb*gwt
       
      xKebe(2,aa,bb) = xKebe(2,aa,bb) + &
        fact2*(-shlu(aa)*mu*shgradgu(bb,1)*nor(2)      &
              - gmul*shgradgu(aa,2)*nor(1)*mu*shlu(bb) &
              + shlu(aa)*nor(1)*munor*nor(2)*shlu(bb))*DetJb*gwt

      xKebe(4,aa,bb) = xKebe(4,aa,bb) +  &
        fact2*(-shlu(aa)*mu*shgradgu(bb,2)*nor(1)      &
              - gmul*shgradgu(aa,1)*nor(2)*mu*shlu(bb) &
              + shlu(aa)*nor(2)*munor*nor(1)*shlu(bb))*DetJb*gwt
       
      xKebe(3,aa,bb) = xKebe(3,aa,bb) +  &
        fact2*(-shlu(aa)*mu*shgradgu(bb,1)*nor(3)      &
              - gmul*shgradgu(aa,3)*nor(1)*mu*shlu(bb) &
              + shlu(aa)*nor(1)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(7,aa,bb) = xKebe(7,aa,bb) +  &
        fact2*(-shlu(aa)*mu*shgradgu(bb,3)*nor(1)      &
              - gmul*shgradgu(aa,1)*nor(3)*mu*shlu(bb) &
              + shlu(aa)*nor(3)*munor*nor(1)*shlu(bb))*DetJb*gwt

      xKebe(6,aa,bb) = xKebe(6,aa,bb) +  &
        fact2*(-shlu(aa)*mu*shgradgu(bb,2)*nor(3)      &
              - gmul*shgradgu(aa,3)*nor(2)*mu*shlu(bb) &
              + shlu(aa)*nor(2)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(8,aa,bb) = xKebe(8,aa,bb) + &
        fact2*(-shlu(aa)*mu*shgradgu(bb,3)*nor(2)      &
              - gmul*shgradgu(aa,2)*nor(3)*mu*shlu(bb) &
              + shlu(aa)*nor(3)*munor*nor(2)*shlu(bb))*DetJb*gwt

    end do    
  end do
  
  do bb = 1, NSHL   
    do aa = 1, NSHL
      xDebe(1,aa,bb) = xDebe(1,aa,bb) - &
                       fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xDebe(2,aa,bb) = xDebe(2,aa,bb) - &
                       fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xDebe(3,aa,bb) = xDebe(3,aa,bb) - &
                       fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
    end do
  end do 


  do bb = 1, NSHL
    do aa = 1, NSHL    
      xGebe(1,aa,bb) = xGebe(1,aa,bb) +  & 
                       fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xGebe(2,aa,bb) = xGebe(2,aa,bb) +  &
                       fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xGebe(3,aa,bb) = xGebe(3,aa,bb) +  &
                       fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt    
    end do
  end do

end subroutine e3bLHS_weak



!======================================================================
! LHS for DG
!======================================================================
subroutine e3bLHS_DG(nshl, shlu, shgradgu, ui, umi, duidxi, tauB, gwt, &
                     nor, xKebe, xGebe, xDebe)  
  use aAdjKeep  
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl  
  real(8), intent(in) :: shlu(NSHL), shgradgu(NSHL,NSD), &
                         ui(NSD), umi(NSD), duidxi(NSD,NSD), &
                         tauB, gwt, nor(NSD)
                         
  real(8), intent(inout) :: xKebe(NSD*NSD,NSHL,NSHL), &
                            xGebe(NSD,NSHL,NSHL), &
                            xDebe(NSD,NSHL,NSHL)

  integer :: aa, bb
  real(8) :: fact1, fact2, tmp1(NSHL), tmp2(NSHL,NSHL), &
             unor, uneg 

  fact1 = almi
  fact2 = alfi*gami*Delt

  tmp1 = 0.0d0
  tmp2 = 0.0d0

  tmp1(:) = shgradgu(:,1)*nor(1) + shgradgu(:,2)*nor(2) &
          + shgradgu(:,3)*nor(3)
  
  ! Diagonal blocks of K

!!$  ! No viscous contributions...
!!$  do bb = 1, NSHL
!!$    do aa = 1, NSHL
!!$      tmp2(aa,bb) = shlu(aa)*tauB*shlu(bb)
!!$    end do
!!$  end do
!!$
!!$  do bb = 1, NSHL    
!!$    do aa = 1, NSHL      
!!$      xKebe(1,aa,bb) = xKebe(1,aa,bb) + fact2*(tmp2(aa,bb))*DetJb*gwt
!!$      xKebe(5,aa,bb) = xKebe(5,aa,bb) + fact2*(tmp2(aa,bb))*DetJb*gwt
!!$      xKebe(9,aa,bb) = xKebe(9,aa,bb) + fact2*(tmp2(aa,bb))*DetJb*gwt
!!$    end do    
!!$  end do

  ! with viscous contributions
  do bb = 1, NSHL
    do aa = 1, NSHL
      tmp2(aa,bb) = -0.5d0*(shlu(aa)*mu*tmp1(bb)    &
                           +tmp1(aa)*mu*shlu(bb))   &
                  + shlu(aa)*tauB*shlu(bb)
    end do
  end do

  do bb = 1, NSHL    
    do aa = 1, NSHL      
      xKebe(1,aa,bb) = xKebe(1,aa,bb) + &
        fact2*(tmp2(aa,bb) &
              -0.5d0*(mu*shgradgu(bb,1)*nor(1)*shlu(aa) &
                     +mu*shgradgu(aa,1)*nor(1)*shlu(bb)))*DetJb*gwt
       
      xKebe(5,aa,bb) = xKebe(5,aa,bb) + &
        fact2*(tmp2(aa,bb) &
              -0.5d0*(mu*shgradgu(bb,2)*nor(2)*shlu(aa) &
                     +mu*shgradgu(aa,2)*nor(2)*shlu(bb)))*DetJb*gwt
       
      xKebe(9,aa,bb) = xKebe(9,aa,bb) + &
        fact2*(tmp2(aa,bb) &
              -0.5d0*(mu*shgradgu(bb,3)*nor(3)*shlu(aa) &
                     +mu*shgradgu(aa,3)*nor(3)*shlu(bb)))*DetJb*gwt
       
      xKebe(2,aa,bb) = xKebe(2,aa,bb) + &
        fact2*(-0.5d0*(mu*shgradgu(bb,1)*nor(2)*shlu(aa) &
                      +mu*shgradgu(aa,2)*nor(1)*shlu(bb)))*DetJb*gwt

      xKebe(4,aa,bb) = xKebe(4,aa,bb) + &
        fact2*(-0.5d0*(mu*shgradgu(bb,2)*nor(1)*shlu(aa) &
                      +mu*shgradgu(aa,1)*nor(2)*shlu(bb)))*DetJb*gwt
       
      xKebe(3,aa,bb) = xKebe(3,aa,bb) + &
        fact2*(-0.5d0*(mu*shgradgu(bb,1)*nor(3)*shlu(aa) &
                      +mu*shgradgu(aa,3)*nor(1)*shlu(bb)))*DetJb*gwt

      xKebe(7,aa,bb) = xKebe(7,aa,bb) + &
        fact2*(-0.5d0*(mu*shgradgu(bb,3)*nor(1)*shlu(aa) &
                      +mu*shgradgu(aa,1)*nor(3)*shlu(bb)))*DetJb*gwt

      xKebe(6,aa,bb) = xKebe(6,aa,bb) + &
        fact2*(-0.5d0*(mu*shgradgu(bb,2)*nor(3)*shlu(aa) &
                      +mu*shgradgu(aa,3)*nor(2)*shlu(bb)))*DetJb*gwt

      xKebe(8,aa,bb) = xKebe(8,aa,bb) + &
        fact2*(-0.5d0*(mu*shgradgu(bb,3)*nor(2)*shlu(aa) &
                      +mu*shgradgu(aa,2)*nor(3)*shlu(bb)))*DetJb*gwt
    end do    
  end do

  ! weighting function and pressure
  do bb = 1, NSHL
    do aa = 1, NSHL
      xGebe(1,aa,bb) = xGebe(1,aa,bb) + & 
                       fact2*0.5d0*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xGebe(2,aa,bb) = xGebe(2,aa,bb) + &
                       fact2*0.5d0*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xGebe(3,aa,bb) = xGebe(3,aa,bb) + &
                       fact2*0.5d0*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt    
    end do
  end do

  ! continuity weighting (q) and velocity
  do bb = 1, NSHL   
    do aa = 1, NSHL
      xDebe(1,aa,bb) = xDebe(1,aa,bb) - &
                       fact2*0.5d0*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xDebe(2,aa,bb) = xDebe(2,aa,bb) - &
                       fact2*0.5d0*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xDebe(3,aa,bb) = xDebe(3,aa,bb) - &
                       fact2*0.5d0*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
    end do
  end do 

end subroutine e3bLHS_DG
