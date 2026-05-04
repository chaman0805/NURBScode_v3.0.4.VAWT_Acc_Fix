!================================================================
! subroutine to compute basis functions for NURBS
!================================================================
subroutine eval_SHAPE_beam(u_hat, shl, shgradl, shhessl, shthdgradl, &
                            tgt, nor, binor, gradnor, gradbinor, &
                            tgt_u, nor_u, binor_u, &
                            xu, xd, dxdxi, ddxddxi, dddxdddxi, &
                            dxdxi_u, ddxddxi_u, &
                            ddu_alpham, du_alphaf, &
                            ddxddxi_psd, ddxddxi_u_psd, psd_switch, psd_u_switch, &
                            p, nsd, nshl, lIEN, nnode,  &
                            B_NET_U, B_NET_D, B_NET_Dt_alphaf, B_NET_DDt_alpham, &
                            DetJb, ni, nuk, U_KNOT)

  implicit none

  integer, intent(in) :: p, nshl, nuk, lIEN(nshl), nnode, nsd, ni
  integer, intent(out) ::psd_switch, psd_u_switch

  ! u and v coordinates of integration point in parent element
  real(8), intent(in) :: u_hat
  real(8), intent(in) :: U_KNOT(nuk)
  real(8), intent(in) :: B_NET_D(nnode,nsd+1), B_NET_U(nnode,nsd+1)
  real(8), intent(in) :: B_NET_Dt_alphaf(nnode,nsd), B_NET_DDt_alpham(nnode,nsd)
                         
  ! Vector of Local basis function values at (u_hat, v_hat), local and 
  ! global gradients
  real(8), intent(out) :: shl(nshl), shgradl(nshl), shhessl(nshl), shthdgradl(nshl), &
                          xu(NSD), xd(NSD), dxdxi(NSD), ddxddxi(NSD), dddxdddxi(NSD), &
                          dxdxi_u(NSD), ddxddxi_u(NSD), &
                          ddxddxi_psd(nsd), ddxddxi_u_psd(NSD), &
                          ddu_alpham(NSD), du_alphaf(NSD), &
                          tgt(NSD), nor(NSD), binor(NSD), gradnor(NSD), gradbinor(NSD), &
                          tgt_u(NSD), nor_u(NSD), binor_u(NSD), DetJb

  ! Local Variables
  real(8) ::temp_vec(NSD), temp_vec1(NSD), Identity_mat(NSD,NSD), temp_mat(NSD,NSD), temp
  ! 1D nonrational basis functions and derivs in u
  real(8) :: N(4,p+1)

  ! u coordinates of integration point, denominator and derivative sums
  real(8) :: u, ee, du, dv, tol1, tol2
  real(8) :: tmpshl(nshl), tmpshgradl(nshl), tmpshhessl(nshl), tmpshthdgradl(nshl)
  real(8) :: shl_sum, shgradl_sum, shhessl_sum, shthdgradl_sum

  ! NURBS coordinates, counters for loops
  integer :: i, j, ct, ii

  ! initialization
  N = 0.0d0; tol1 = 1.0d-10; tol2 = 1.0d-6  
  shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0     
  tmpshl = 0.0d0; tmpshgradl = 0.0d0; tmpshhessl = 0.0d0

  ! Get u coordinates of integration point
  u = ((U_KNOT(ni+1)-U_KNOT(ni))*u_hat + U_KNOT(ni+1)+U_KNOT(ni))/2.0d0        
        
  ! Evaluate 1D shape functions and derivatives each direction
  ! calculate in u direction
  call dersbasisfuns_shell(ni, p, nuk, u, 3, U_KNOT, N)    
  !write(*,*)ni,N(1,:)
  !------------------------------------------------------------
  ! Form basis functions and derivatives "dR/du" and "dR/dv"
  !------------------------------------------------------------
  ct = 0
    do i = 0, p
      ct = ct + 1
   !   write(*,*) B_NET_U(lIEN(ct),NSD+1),N(1,p+1-i),N(2,p+1-i),N(3,p+1-i)
      ! basis functions
      tmpshl(ct) = N(1,p+1-i)*1.0*B_NET_U(lIEN(ct),NSD+1)

      ! first derivatives
      tmpshgradl(ct) = N(2,p+1-i)*1.0*B_NET_U(lIEN(ct),NSD+1) ! ,u
  
      ! second derivatives
      tmpshhessl(ct) = N(3,p+1-i)*1.0*B_NET_U(lIEN(ct),NSD+1) ! ,uu

      ! third derivatives
      tmpshthdgradl(ct) = N(4,p+1-i)*1.0*B_NET_U(lIEN(ct),NSD+1) ! ,uuu

    end do

  shl_sum = sum(tmpshl)

  shgradl_sum = sum(tmpshgradl(:))

  shhessl_sum = sum(tmpshhessl(:))

  shthdgradl_sum = sum(tmpshthdgradl(:))

  ! Divide through by denominator
  shl = tmpshl/shl_sum

  shgradl(:) = tmpshgradl(:)/shl_sum -&
                 (tmpshl(:)*shgradl_sum)/(shl_sum**2.0d0)

  shhessl(:) = tmpshhessl(:)/shl_sum - &
                 tmpshgradl(:)*shgradl_sum/(shl_sum**2.0d0) -&
                 ((tmpshgradl(:)*shgradl_sum+tmpshl(:)*shhessl_sum)/&
                 (shl_sum**2.0d0)&
                  -2.0d0*tmpshl(:)*shgradl_sum*shgradl_sum/&
                 (shl_sum**3.0d0))

  shthdgradl(:) = tmpshthdgradl(:)/shl_sum - &
                  3.0d0*tmpshhessl(:)/(shl_sum**2.0d0)*shgradl_sum + &
                  6.0d0*tmpshgradl(:)/(shl_sum**3.0d0)*(shgradl_sum**2.0d0) - &
                  3.0d0*tmpshgradl(:)/(shl_sum**2.0d0)*shhessl_sum - &
                  6.0d0*tmpshl(:)/(shl_sum**4.0d0)*(shgradl_sum**3.0d0) - &
                  1.0d0*tmpshl(:)/(shl_sum**2.0d0)*(shthdgradl_sum) + &
                  6.0d0*tmpshl(:)/(shl_sum**3.0d0)*shgradl_sum*shhessl_sum
  !write(*,*) shthdgradl(:)


  ! Now we calculate the face Jacobian
  xu = 0.0d0
  xd = 0.0d0
  dxdxi   = 0.0d0
  ddxddxi = 0.0d0
  dddxdddxi = 0.0d0
  dxdxi_u   = 0.0d0
  ddxddxi_u = 0.0d0
  ddu_alpham = 0.0d0
  du_alphaf = 0.0d0
  ddxddxi_psd=0.0d0
  ddxddxi_u_psd= 0.0d0
  ! pseudo vector mode is by default open
  psd_switch = 1
  psd_u_switch = 1

  do i = 1, nshl
  ! write(*,*)B_NET_U(lIEN(i),:),shl(i),shgradl(i)
    do ii = 1, 3
      xu(ii) = xu(ii) + B_NET_U(lIEN(i),ii)*shl(i)
      xd(ii) = xd(ii) + B_NET_D(lIEN(i),ii)*shl(i)
      dxdxi(ii)   = dxdxi(ii)   + B_NET_D(lIEN(i),ii)*shgradl(i)
      ! write(*,*) dxdxi(:)
      ddxddxi(ii) = ddxddxi(ii) + B_NET_D(lIEN(i),ii)*shhessl(i)
      dddxdddxi(ii) = dddxdddxi(ii) + B_NET_D(lIEN(i),ii)*shthdgradl(i)

      dxdxi_u(ii)   = dxdxi_u(ii)   + B_NET_U(lIEN(i),ii)*shgradl(i)
      ddxddxi_u(ii) = ddxddxi_u(ii) + B_NET_U(lIEN(i),ii)*shhessl(i)

      ddu_alpham(ii) = ddu_alpham(ii) + B_NET_DDt_alpham(lIEN(i),ii)*shl(i)
      du_alphaf(ii) = du_alphaf(ii) + B_NET_Dt_alphaf(lIEN(i),ii)*shl(i)
    end do
  end do


  ee = dxdxi(1)**2+dxdxi(2)**2+dxdxi(3)**2
  
  tgt = 0.0d0
  nor = 0.0d0
  binor= 0.0d0
  gradnor = 0.0d0
  gradbinor = 0.0d0
  tgt_u = 0.0d0
  nor_u = 0.0d0
  binor_u = 0.0d0
  DetJb = sqrt(ee) ! Jacobian of face mapping
  !write(*,*) dxdxi, DetJb
  !tangent, normal and binormal vectors tgt(NSD), nor(NSD), binor(NSD)
  tgt(:)=dxdxi(:)/sqrt(sum(dxdxi(:)**2.0d0))
  call cross_prod(dxdxi,ddxddxi,temp_vec)
  !write(*,*) "temp_vec=", temp_vec

  if (sqrt(sum(ddxddxi(:)**2))> tol1) then
    temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi(:)**2.0d0))*sqrt(sum(ddxddxi(:)**2.0d0)))
    if (abs(temp-0.0d0)> tol2)  then
      psd_switch = 0
    endif
  endif

  if (psd_switch .eq. 1) then
      ddxddxi_psd(1)=0.0d0; ddxddxi_psd(2)=0.0d0; ddxddxi_psd(3)=1.0d0
      call cross_prod(dxdxi,ddxddxi_psd,temp_vec)
      temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi(:)**2.0d0))*1.0d0)
      if (abs(temp-0.0d0)<=tol2) then
        ddxddxi_psd(1)=0.0d0; ddxddxi_psd(2)=1.0d0; ddxddxi_psd(3)=0.0d0
        call cross_prod(dxdxi,ddxddxi_psd,temp_vec)
        temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi(:)**2.0d0))*1.0d0)
        if (abs(temp-0.0d0)<=tol2) then
          ddxddxi_psd(1)=1.0d0; ddxddxi_psd(2)=0.0d0; ddxddxi_psd(3)=0.0d0
          call cross_prod(dxdxi,ddxddxi_psd,temp_vec)
          temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi(:)**2.0d0))*1.0d0)
        endif
      endif 
  endif

  binor(:)=temp_vec(:)/sqrt(sum(temp_vec(:)**2.0d0))
  call cross_prod(binor,dxdxi,temp_vec)
  nor(:)=temp_vec(:)/sqrt(sum(temp_vec(:)**2.0d0))
  !write(*,*) "xd=", xd
  !write(*,*) "dxdxi=", dxdxi
  !write(*,*) "ddxdxdi=", ddxddxi
  !write(*,*) "tgt=", tgt
  !write(*,*) "nor1=", nor
  !write(*,*) "binor1=", binor

  !gradnor(NSD), gradbinor(NSD)
  call cross_prod(dxdxi_u,ddxddxi_u,temp_vec)
  !write(*,*) "temp_vec=", temp_vec
  temp=sqrt(sum(temp_vec(:)**2.0d0))

  if (psd_switch .eq. 0)  then
    call tensor_prod(nor,nor,temp_mat)
    call create_KDelta(Identity_mat)
    temp_mat=Identity_mat-temp_mat
    call cross_prod(binor,ddxddxi,temp_vec)
    call cross_prod(binor,dxdxi,temp_vec1)
    temp_vec(:)=temp_vec(:)/sqrt(sum(temp_vec1(:)**2.0d0))
    call Mat_prod(temp_vec,temp_mat,gradnor)

    call tensor_prod(binor,binor,temp_mat)
    temp_mat=Identity_mat-temp_mat
    call cross_prod(dxdxi,dddxdddxi,temp_vec)
    call cross_prod(dxdxi,ddxddxi,temp_vec1)
    temp_vec(:)=temp_vec(:)/sqrt(sum(temp_vec1(:)**2.0d0))
    !write(*,*) "temp_vec=", temp_vec
    call Mat_prod(temp_vec,temp_mat,gradbinor)
  else
    call tensor_prod(nor,nor,temp_mat)
    call create_KDelta(Identity_mat)
    temp_mat=Identity_mat-temp_mat
    call cross_prod(binor,ddxddxi,temp_vec)
    call cross_prod(binor,dxdxi,temp_vec1)
    temp_vec(:)=temp_vec(:)/sqrt(sum(temp_vec1(:)**2.0d0))
    call Mat_prod(temp_vec,temp_mat,gradnor)

    call tensor_prod(binor,binor,temp_mat)
    temp_mat=Identity_mat-temp_mat
    call cross_prod(ddxddxi,ddxddxi_psd,temp_vec)
    call cross_prod(dxdxi,ddxddxi_psd,temp_vec1)
    temp_vec(:)=temp_vec(:)/sqrt(sum(temp_vec1(:)**2.0d0))
    !write(*,*) "temp_vec=", temp_vec
    call Mat_prod(temp_vec,temp_mat,gradbinor)
  endif

  !verify
  call dot_prod(gradnor,dxdxi,temp)
  !write(*,*) "gn*gx=", temp
  call dot_prod(nor,ddxddxi,temp)
  !write(*,*) "n*ggx=", temp

  !tangent, normal and binormal vectors in ref config. tgt_u(NSD), nor_u(NSD), binor_u(NSD)
  tgt_u(:)=dxdxi_u(:)/sqrt(sum(dxdxi_u(:)**2.0d0))


  call cross_prod(dxdxi_u,ddxddxi_u,temp_vec)
  !write(*,*) "temp_vec=", temp_vec
  if (sqrt(sum(ddxddxi_u(:)**2))> tol1) then
    temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi_u(:)**2.0d0))*sqrt(sum(ddxddxi_u(:)**2.0d0)))
    if (abs(temp-0.0d0)> tol2)  then
      psd_u_switch = 0
    endif
  endif

  if (psd_u_switch .eq. 1) then
      ddxddxi_u_psd(1)=0.0d0; ddxddxi_u_psd(2)=0.0d0; ddxddxi_u_psd(3)=1.0d0
      call cross_prod(dxdxi_u,ddxddxi_u_psd,temp_vec)
      temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi_u(:)**2.0d0))*1.0d0)
      if (abs(temp-0.0d0)<=tol2) then
        ddxddxi_u_psd(1)=0.0d0; ddxddxi_u_psd(2)=1.0d0; ddxddxi_u_psd(3)=0.0d0
        call cross_prod(dxdxi_u,ddxddxi_u_psd,temp_vec)
        temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi_u(:)**2.0d0))*1.0d0)
        if (abs(temp-0.0d0)<=tol2) then
          ddxddxi_u_psd(1)=1.0d0; ddxddxi_u_psd(2)=0.0d0; ddxddxi_u_psd(3)=0.0d0
          call cross_prod(dxdxi_u,ddxddxi_u_psd,temp_vec)
          temp=sqrt(sum(temp_vec(:)**2.0d0))/(sqrt(sum(dxdxi_u(:)**2.0d0))*1.0d0)
        endif
      endif 
  endif

  binor_u(:)=temp_vec(:)/sqrt(sum(temp_vec(:)**2.0d0))
  call cross_prod(binor_u,dxdxi_u,temp_vec)
  nor_u(:)=temp_vec(:)/sqrt(sum(temp_vec(:)**2.0d0))
  !write(*,*) "nor_u=", nor_u

  !------------------------------------------
  ! change from dx/du to dxdxi
  !------------------------------------------
  ! Get knot span sizes
  du = U_KNOT(ni+1)-U_KNOT(ni)
  !write(*,*) dxdxi(:)
  !dxdxi(:) = dxdxi(:)*du/2.0d0

  !ddxddxi(:) = ddxddxi(:)*du*du/4.0d0

  !dddxdddxi(:) = dddxdddxi(:)*du*du*du/8.0d0

  !dxdxi_u(:) = dxdxi_u(:)*du/2.0d0

  !ddxddxi_u(:) = ddxddxi_u(:)*du*du/4.0d0



end subroutine eval_SHAPE_beam

