subroutine e3LRhs_KLBeam(shl, shgradl, shhessl, ptype, Ec, Rad, Dens, C_dp, &
                          xKebe, Rhs, Rhs_ext, nshl, NSD, &
                          tgt, nor, binor, &
                          tgt_u, nor_u, binor_u, &
                          xu, xd, dxdxi, ddxddxi, &
                          dxdxi_u, ddxddxi_u, &
                          ddu_alpham, du_alphaf, &
                          ddxddxi_psd, ddxddxi_u_psd, psd_switch, psd_u_switch, &
                          B_NET_SH, B_NET_SH_D, &
                          lIEN, nnode, Delta_t, &
                          alpha_m, alpha_f, gamma_t, beta_t)

  implicit none  

  integer, intent(in) :: ptype, nshl, NSD, lIEN(nshl), nnode, psd_switch, psd_u_switch  
  real(8), intent(in) :: shl(NSHL), shgradl(NSHL), shhessl(NSHL), Ec, Rad, Dens, C_dp, &
                         tgt(NSD), nor(NSD), binor(NSD), &
                         tgt_u(NSD), nor_u(NSD), binor_u(NSD), &
                         xu(NSD), xd(NSD), dxdxi(NSD), ddxddxi(nsd), &
                         dxdxi_u(NSD), ddxddxi_u(NSD), &
                         ddu_alpham(NSD), du_alphaf(NSD), &
                         ddxddxi_psd(nsd), ddxddxi_u_psd(NSD), &
                         B_NET_SH(nnode,nsd+1), B_NET_SH_D(nnode,nsd+1), &
                         Delta_t, alpha_m, alpha_f, gamma_t, beta_t
                         
  real(8), intent(out):: xKebe(NSD*NSD,NSHL,NSHL), Rhs(NSD,NSHL), Rhs_ext(NSD,NSHL)

  integer :: i, j, k, l, m, p, q, s, A, B
  real(8) :: pi, temp1, temp2, temp3, temp, temp0, temp_vec(NSD), temp_vec1(NSD)
  real(8) :: norm_G1, xKebe_mem(NSD*NSD,NSHL,NSHL), xKebe_bd(NSD*NSD,NSHL,NSHL), &
             xKebe_ddu(NSD*NSD,NSHL,NSHL), xKebe_du(NSD*NSD,NSHL,NSHL), &
             Rhs_mem(NSD,NSHL), Rhs_bd(NSD,NSHL), &
             Identity_mat(NSD,NSD), PermSymb(NSD,NSD,NSD), &
             A_mat(NSD,NSD), B_mat(NSD,NSD), &
             nk_uAi(NSD,NSD,NSHL), bl_uBj(NSD,NSD,NSHL), &
             Aks_uBj(NSD,NSD,NSD,NSHL), Blm_uBj(NSD,NSD,NSD,NSHL), &
             nk_uAiBj(NSD,NSD,NSD,NSHL,NSHL)
  real(8) :: Ec1

  pi = acos(-1.0d0)  

  Rhs=0.0d0; Rhs_ext=0.0d0; xKebe=0.0d0; xKebe_ddu=0.0d0; xKebe_du=0.0d0
  Rhs_mem=0.0d0; Rhs_bd=0.0d0; xKebe_mem=0.0d0; xKebe_bd=0.0d0  

  call create_KDelta(Identity_mat)
  norm_G1=1.0d0/sqrt(sum(dxdxi_u(:)**2.0d0))
  temp1=sum(dxdxi**2.0d0)-sum(dxdxi_u**2.0d0)
  Ec1=Ec
  !if ((temp1-0.0d0)<-1.0d-15 .and. ptype .eq. 1) Ec1=0.0d0
  if ((ptype .eq. 2) .or. (ptype .eq. 4) ) Ec1=Ec*1.0d3

  do A = 1, nshl
    do i = 1, nsd
      Rhs_ext(i,A)=-Dens*pi*(Rad**2.0d0)*shl(A)*ddu_alpham(i)- Dens*pi*(Rad**2.0d0)*shl(A)*C_dp*du_alphaf(i)!
      
    enddo
  enddo

  do A = 1, nshl
    do B = 1, nshl
      do i = 1, nsd
        do j = 1, nsd
          xKebe_ddu((i-1)*nsd+j,A,B)= Dens*pi*(Rad**2.0d0)*shl(A)*shl(B)*Identity_mat(i,j)             
        enddo
      enddo   
    enddo
  enddo

  do A = 1, nshl
    do B = 1, nshl
      do i = 1, nsd
        do j = 1, nsd
          xKebe_du((i-1)*nsd+j,A,B)=Dens*pi*(Rad**2.0d0)*C_dp*shl(A)*shl(B)*Identity_mat(i,j) 
        enddo
      enddo   
    enddo
  enddo

 
  norm_G1=1.0d0/sqrt(sum(dxdxi_u(:)**2.0d0))
  temp1=sum(dxdxi**2.0d0)-sum(dxdxi_u**2.0d0)
  !write(*,*) "norm_G1=", norm_G1
  !write(*,*) "tgt=", tgt
  !write(*,*) "nor=", nor
  !write(*,*) "dxdxi=", dxdxi
  !write(*,*) "ddxddxi=", ddxddxi
  !write(*,*) "binor=", binor
  !write(*,*) "psd_switch", psd_switch
  !write(*,*) "dxdxi_u=", dxdxi_u
  !write(*,*) "Ec=", Ec
  !write(*,*) "Rad=", Rad


  do A = 1, nshl
    do i = 1, nsd
      Rhs_mem(i,A)=Rhs_mem(i,A)-(norm_G1**4.0d0)*Ec1*pi/2.0d0*(Rad**2.0d0)*temp1*dxdxi(i)*shgradl(A) !
    enddo
  enddo

  !do i = 1, nshl
  !   write(*,*) "Rhs=", Rhs(:,i)
  !   write(*,*) "Rhs_mem=",(norm_G1**4.0d0)*Ec1*pi/2.0d0*(Rad**2.0d0)*temp1, Rhs_mem(:,i)
  !   write(*,*) "Rhs_bd=", Rhs_bd(:,i)
  !enddo

  do A = 1, nshl
    do B = 1, nshl
      do i = 1, nsd
        do j = 1, nsd
          !change sign here
          xKebe_mem((i-1)*nsd+j,A,B)= xKebe_mem((i-1)*nsd+j,A,B)+(norm_G1**4.0d0)*Ec1*pi/2.0d0*(Rad**2.0d0) &
                                      *(2.0d0*dxdxi(i)*dxdxi(j)*shgradl(A)*shgradl(B) &
                                      + temp1*shgradl(A)*shgradl(B)*Identity_mat(i,j))
          !write(*,*) A,B,(i-1)*nsd+j,'xKebe_mem=', xKebe_mem((i-1)*nsd+j,A,B)
          !write(*,*) 'shgradl=', shgradl(A), shgradl(B)
          !write(*,'(/)')
        enddo
      enddo   
    enddo
  enddo

  A_mat=0.0d0; B_mat=0.0d0
  nk_uAi=0.0d0; bl_uBj=0.0d0
  Aks_uBj=0.0d0; Blm_uBj=0.0d0
  nk_uAiBj=0.0d0

  if (psd_switch .eq. 0)  then

  call tensor_prod(nor,nor,A_mat)
  A_mat=Identity_mat-A_mat
  call cross_prod(binor,dxdxi,temp_vec1)
  A_mat=A_mat/sqrt(sum(temp_vec1(:)**2.0d0))
  !write(*,*) A_mat

  call tensor_prod(binor,binor,B_mat)
  B_mat=Identity_mat-B_mat
  call cross_prod(dxdxi,ddxddxi,temp_vec1)
  B_mat=B_mat/sqrt(sum(temp_vec1(:)**2.0d0))
  !write(*,*) B_mat


  call create_PermSymb(PermSymb)

  do A = 1, nshl
    do k = 1, nsd
      do i = 1, nsd
        do s = 1, nsd
          do l = 1, nsd
            do p = 1, nsd
              do q = 1, nsd
                do m = 1, nsd
                  nk_uAi(k,i,A)=nk_uAi(k,i,A)+A_mat(k,s)*PermSymb(s,l,q)*B_mat(l,m)*(PermSymb(m,i,p)*shgradl(A)*ddxddxi(p) &
                                +PermSymb(m,p,i)*shhessl(A)*dxdxi(p))*dxdxi(q)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do A = 1, nshl
    do k = 1, nsd
      do i = 1, nsd
        do s = 1, nsd
          do l = 1, nsd
            nk_uAi(k,i,A)=nk_uAi(k,i,A)+A_mat(k,s)*PermSymb(s,l,i)*binor(l)*shgradl(A)    
          enddo
        enddo
      enddo
    enddo
  enddo


  do B = 1, nshl
    do l = 1, nsd
      do j = 1, nsd
        do m = 1, nsd
          do p = 1, nsd
            bl_uBj(l,j,B)=bl_uBj(l,j,B)+B_mat(l,m)*(PermSymb(m,j,p)*ddxddxi(p)*shgradl(B) &
                          +PermSymb(m,p,j)*dxdxi(p)*shhessl(B))    
          enddo
        enddo
      enddo
    enddo
  enddo


  call dot_prod(binor,binor,temp1)
  call dot_prod(dxdxi,dxdxi,temp2)
  call dot_prod(binor,dxdxi,temp3)
  temp0=(temp1*temp2-temp3**2.0d0)**(0.5d0)

  do B = 1, nshl
    do k = 1, nsd
      do s = 1, nsd
        do j = 1, nsd
          temp=0.0d0
          do p = 1, nsd
            temp= temp-(bl_uBj(p,j,B)*binor(p)*temp2+ &
                  binor(p)*binor(p)*shgradl(B)*dxdxi(j)- bl_uBj(p,j,B)*dxdxi(p)*temp3- &
                  binor(j)*shgradl(B)*binor(p)*dxdxi(p))                         
          enddo
          Aks_uBj(k,s,j,B)=Aks_uBj(k,s,j,B)+ temp*(Identity_mat(k,s)-nor(k)*nor(s))/ (temp0**3.0d0)- &
                            (nk_uAi(k,j,B)*nor(s)+nk_uAi(s,j,B)*nor(k))/temp0
        enddo
      enddo
    enddo
  enddo

  call dot_prod(dxdxi,dxdxi,temp1)
  call dot_prod(ddxddxi,ddxddxi,temp2)
  call dot_prod(dxdxi,ddxddxi,temp3)
  temp0=(temp1*temp2-temp3**2.0d0)**(0.5d0)
  !write (*,*) 'temp0=',temp0

  do B = 1, nshl
    do l = 1, nsd
      do m = 1, nsd
        do j = 1, nsd
          temp=0.0d0
          do k = 1, nsd
            temp= temp-(shgradl(B)*dxdxi(j)*ddxddxi(k)*ddxddxi(k)+ &
                 shhessl(B)*ddxddxi(j)*dxdxi(k)*dxdxi(k)- shgradl(B)*ddxddxi(j)*dxdxi(k)*ddxddxi(k)- &
                             shhessl(B)*dxdxi(j)*dxdxi(k)*ddxddxi(k))                          
          enddo
          Blm_uBj(l,m,j,B)= Blm_uBj(l,m,j,B)+temp*(Identity_mat(l,m)-binor(l)*binor(m))/ (temp0**3.0d0)- &
                            (bl_uBj(l,j,B)*binor(m)+bl_uBj(m,j,B)*binor(l))/temp0
        enddo
      enddo
    enddo
  enddo

  do A = 1, nshl
    do B = 1, nshl
      do k = 1, nsd
        do i = 1, nsd
          do j = 1, nsd
            do s = 1, nsd
              do l = 1, nsd
                do p = 1, nsd
                  do q = 1, nsd
                    do m = 1, nsd
                      nk_uAiBj(k,i,j,A,B)=nk_uAiBj(k,i,j,A,B)+Aks_uBj(k,s,j,B)*PermSymb(s,l,q)*B_mat(l,m)*(PermSymb(m,i,p) &
                                          *shgradl(A)*ddxddxi(p) +PermSymb(m,p,i)*shhessl(A)*dxdxi(p))*dxdxi(q)+ &
                                          A_mat(k,s)*PermSymb(s,l,q)*Blm_uBj(l,m,j,B)*(PermSymb(m,i,p) &
                                          *shgradl(A)*ddxddxi(p) +PermSymb(m,p,i)*shhessl(A)*dxdxi(p))*dxdxi(q)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do A = 1, nshl
    do B = 1, nshl
      do k = 1, nsd
        do i = 1, nsd
          do j = 1, nsd
            do s = 1, nsd
              do l = 1, nsd
                nk_uAiBj(k,i,j,A,B)=nk_uAiBj(k,i,j,A,B)+Aks_uBj(k,s,j,B)*PermSymb(s,l,i)*binor(l)*shgradl(A)    
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo


  do A = 1, nshl
    do B = 1, nshl
      do k = 1, nsd
        do i = 1, nsd
          do j = 1, nsd
            do s = 1, nsd
              do l = 1, nsd
                do q = 1, nsd
                  do m = 1, nsd
                    nk_uAiBj(k,i,j,A,B)=nk_uAiBj(k,i,j,A,B)+A_mat(k,s)*PermSymb(s,l,q)*B_mat(l,m)*(PermSymb(m,i,j) &
                                        *shgradl(A)*shhessl(B) +PermSymb(m,j,i)*shhessl(A)*shgradl(B))*dxdxi(q)+ &
                                        A_mat(k,s)*PermSymb(s,l,j)*B_mat(l,m)*(PermSymb(m,i,q) &
                                        *shgradl(A)*ddxddxi(q) +PermSymb(m,q,i)*shhessl(A)*dxdxi(q))*shgradl(B)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo


  !write(*,*) "nk_uAiBj(k,i,j,A,B)=", nk_uAiBj(2,3,3,3,1), nk_uAiBj(2,3,3,1,3)
  
  do A = 1, nshl
    do B = 1, nshl
      do k = 1, nsd
        do i = 1, nsd
          do j = 1, nsd
            do s = 1, nsd
              do l = 1, nsd
                nk_uAiBj(k,i,j,A,B)=nk_uAiBj(k,i,j,A,B)+A_mat(k,s)*PermSymb(s,l,i)*bl_uBj(l,j,B)*shgradl(A)    
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  !write(*,*) "nk_uAiBj(k,i,j,A,B)=", nk_uAiBj(2,3,3,3,1), nk_uAiBj(2,3,3,1,3)

  else

  call tensor_prod(nor,nor,A_mat)
  A_mat=Identity_mat-A_mat
  call cross_prod(binor,dxdxi,temp_vec1)
  A_mat=A_mat/sqrt(sum(temp_vec1(:)**2.0d0))
  !write(*,*) A_mat

  call tensor_prod(binor,binor,B_mat)
  B_mat=Identity_mat-B_mat
  call cross_prod(dxdxi,ddxddxi_psd,temp_vec1)
  B_mat=B_mat/sqrt(sum(temp_vec1(:)**2.0d0))
  !write(*,*) B_mat

  call create_PermSymb(PermSymb)

  do A = 1, nshl
    do k = 1, nsd
      do i = 1, nsd
        do s = 1, nsd
          do l = 1, nsd
            do p = 1, nsd
              do q = 1, nsd
                do m = 1, nsd
                  nk_uAi(k,i,A)=nk_uAi(k,i,A)+A_mat(k,s)*PermSymb(s,l,q)*B_mat(l,m)*(PermSymb(m,i,p)*shgradl(A)*ddxddxi_psd(p) &
                                )*dxdxi(q)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do A = 1, nshl
    do k = 1, nsd
      do i = 1, nsd
        do s = 1, nsd
          do l = 1, nsd
            nk_uAi(k,i,A)=nk_uAi(k,i,A)+A_mat(k,s)*PermSymb(s,l,i)*binor(l)*shgradl(A)    
          enddo
        enddo
      enddo
    enddo
  enddo

  endif
  

  call dot_prod(ddxddxi_u,nor_u,temp1)
  !write(*,*) "ddxddxi_u=", ddxddxi_u, "nor_u", nor_u
  !write(*,*) "whole=", temp1
  call dot_prod(ddxddxi,nor,temp2)
  !write(*,*) "whole=", temp2
  temp1=temp1-temp2

  do A = 1, nshl
    do i = 1, nsd
      Rhs_bd(i,A)=Rhs_bd(i,A)+(norm_G1**4.0d0)*Ec1*pi/4.0d0*(Rad**4.0d0)*temp1*nor(i)*shhessl(A) !
      !write(*,*) "whole=", temp1
      !write(*,*) "Rhs_bd1=", Rhs_bd(i,A)
    enddo
  enddo

  do A = 1, nshl
    do i = 1, nsd
      do k= 1, nsd
        Rhs_bd(i,A)=Rhs_bd(i,A)+(norm_G1**4.0d0)*Ec1*pi/4.0d0*(Rad**4.0d0)*temp1*ddxddxi(k)*nk_uAi(k,i,A) !
        !write(*,*) "Rhs_bd2=", Rhs_bd(i,A)
        !write(*,*) "whole2=", (norm_G1**4.0d0)*Ec1*pi/4.0d0*(Rad**4.0d0)*temp1*ddxddxi(k)
        !write(*,*) "nk_uAi=", nk_uAi(k,i,A)
      enddo
    enddo
  enddo

  do A = 1, nshl
    do B = 1, nshl
      do i = 1, nsd
        do j = 1, nsd
          temp=0.0d0
          temp0=0.0d0
          do k = 1, nsd
            temp=temp+ddxddxi(k)*nk_uAi(k,j,B)
            temp0=temp0+ddxddxi(k)*nk_uAi(k,i,A)
          enddo
          !change sign here
          xKebe_bd((i-1)*nsd+j,A,B)=xKebe_bd((i-1)*nsd+j,A,B)+(norm_G1**4.0d0)*Ec1*pi/4.0d0*(Rad**4.0d0)* &
                            (shhessl(B)*nor(j)+temp)*(shhessl(A)*nor(i)+temp0)
        enddo
      enddo    
    enddo
  enddo

  !temp1=1.0
  do A = 1, nshl
    do B = 1, nshl
      do i = 1, nsd
        do j = 1, nsd
          !change sign here
          xKebe_bd((i-1)*nsd+j,A,B)=xKebe_bd((i-1)*nsd+j,A,B)-(norm_G1**4.0d0)*Ec1*pi/4.0d0*(Rad**4.0d0)*temp1* &
                              (shhessl(A)*nk_uAi(i,j,B)+shhessl(B)*nk_uAi(j,i,A))
          !write(*,*) xKebe_bd((i-1)*nsd+j,A,B)
        enddo
      enddo    
    enddo
  enddo

  do A = 1, nshl
    do B = 1, nshl
      do i = 1, nsd
        do j = 1, nsd
          do k= 1, nsd
            !change sign here
            xKebe_bd((i-1)*nsd+j,A,B)=xKebe_bd((i-1)*nsd+j,A,B)-(norm_G1**4.0d0)*Ec1*pi/4.0d0*(Rad**4.0d0)*temp1* &
                              ddxddxi(k)*nk_uAiBj(k,i,j,A,B)
            !write(*,*) xKebe_bd((i-1)*nsd+j,A,B)
          enddo
        enddo
      enddo    
    enddo
  enddo

  !write(*,*) "Xkebe_bd=", xKebe_bd(1*2,2,1),xKebe_bd(1*2,1,2)
  !Rhs_bd=0.0d0;  xKebe_bd=0.0d0  
  xKebe=(xKebe_mem+xKebe_bd)*alpha_f*beta_t*(Delta_t**2.0d0)+xKebe_du*alpha_f*gamma_t*Delta_t+xKebe_ddu*alpha_m
  Rhs=Rhs_mem+Rhs_bd+Rhs_ext
  if ((ptype .eq. 2) .or. (ptype .eq. 4) ) then
    xKebe=(xKebe_bd)*alpha_f*beta_t*(Delta_t**2.0d0)
    Rhs=Rhs_bd
  end if

  !write(*,*) "Xkebe=", xKebe(1*2,2,1),xKebe(1*2,1,2)
  !do i = 1, nshl
  !   write(*,*) "Rhs=", Rhs(:,i)
  !   write(*,*) "Rhs_mem=", Rhs_mem(:,i)
  !   write(*,*) "Rhs_bd=", Rhs_bd(:,i)
  !enddo

  !do A = 1, nshl
  !  do B = 1, nshl
  !    do i = 1, nsd*nsd
         !write(*,*) A,B,i,'xKebe=', xKebe(i,A,B)
         !write(*,'(/)')
  !    enddo   
  !  enddo
  !enddo

end subroutine e3LRhs_KLBeam
