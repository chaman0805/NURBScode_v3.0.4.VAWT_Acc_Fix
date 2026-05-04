!================================================================
! subroutine to compute basis functions for NURBS
!================================================================
subroutine eval_SHAPE_shell(u_hat, v_hat, shl, shgradl, shhessl, &
                            nor, xu, xd, dxdxi, ddxddxi,   &
                            p, q, nsd, nshl, lIEN, nnode,  &
                            B_NET_U, B_NET_D, DetJb, ni, nj, &
                            nuk, nvk, U_KNOT, V_KNOT)

  implicit none

  integer, intent(in) :: p, q, nshl, nuk, nvk, lIEN(nshl), nnode, nsd, &
                         ni, nj

  ! u and v coordinates of integration point in parent element
  real(8), intent(in) :: u_hat, v_hat
  real(8), intent(in) :: U_KNOT(nuk), V_KNOT(nvk)
  real(8), intent(in) :: B_NET_D(nnode,nsd+1), B_NET_U(nnode,nsd+1)
                         
  ! Vector of Local basis function values at (u_hat, v_hat), local and 
  ! global gradients
  real(8), intent(out) :: shl(nshl), shgradl(nshl,2), shhessl(nshl,3),    &
                          xu(NSD), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3), &
                          nor(NSD), DetJb

  ! Local Variables
  ! 1D nonrational basis functions and derivs in u and v
  real(8) :: N(3,p+1), M(3,q+1)

  ! u and v coordinates of integration point, denominator and derivative sums
  real(8) :: u, v, ee, ff, gg, du, dv
  real(8) :: tmpshl(nshl), tmpshgradl(nshl,2), tmpshhessl(nshl,3)
  real(8) :: shl_sum, shgradl_sum(2), shhessl_sum(3)

  ! NURBS coordinates, counters for loops
  integer :: i, j, ct, ii

  ! initialization
  M = 0.0d0; N = 0.0d0     
  shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0     
  tmpshl = 0.0d0; tmpshgradl = 0.0d0; tmpshhessl = 0.0d0

  ! Get u and v coordinates of integration point
  u = ((U_KNOT(ni+1)-U_KNOT(ni))*u_hat + U_KNOT(ni+1)+U_KNOT(ni))/2.0d0        
  v = ((V_KNOT(nj+1)-V_KNOT(nj))*v_hat + V_KNOT(nj+1)+V_KNOT(nj))/2.0d0
        
  ! Evaluate 1D shape functions and derivatives each direction
  ! calculate in u and v direction
  call dersbasisfuns_shell(ni, P, nuk, u, 2, U_KNOT, N)
  call dersbasisfuns_shell(nj, Q, nvk, v, 2, V_KNOT, M)    
  
  !------------------------------------------------------------
  ! Form basis functions and derivatives "dR/du" and "dR/dv"
  !------------------------------------------------------------
  ct = 0
  do j = 0, Q
    do i = 0, P
      ct = ct + 1

      ! basis functions
      tmpshl(ct) = N(1,P+1-i)*M(1,Q+1-j)*B_NET_U(lIEN(ct),NSD+1)

      ! first derivatives
      tmpshgradl(ct,1) = N(2,p+1-i)*M(1,Q+1-j)*B_NET_U(lIEN(ct),NSD+1) ! ,u
      tmpshgradl(ct,2) = N(1,p+1-i)*M(2,Q+1-j)*B_NET_U(lIEN(ct),NSD+1) ! ,v
    end do
  end do

  ! Now, second derivatives
  ct = 0
  do j = 0, Q
    do i = 0, P
      ct = ct + 1
      tmpshhessl(ct,1) = N(3,P+1-i)*M(1,Q+1-j)*B_NET_U(lIEN(ct),NSD+1) ! ,uu
      tmpshhessl(ct,2) = N(2,P+1-i)*M(2,Q+1-j)*B_NET_U(lIEN(ct),NSD+1) ! ,uv
      tmpshhessl(ct,3) = N(1,P+1-i)*M(3,Q+1-j)*B_NET_U(lIEN(ct),NSD+1) ! ,vv
    end do
  end do

  shl_sum = sum(tmpshl)

  shgradl_sum(1) = sum(tmpshgradl(:,1))
  shgradl_sum(2) = sum(tmpshgradl(:,2))

  shhessl_sum(1) = sum(tmpshhessl(:,1))
  shhessl_sum(2) = sum(tmpshhessl(:,2))
  shhessl_sum(3) = sum(tmpshhessl(:,3))

  ! Divide through by denominator
  shl = tmpshl/shl_sum

  shgradl(:,1) = tmpshgradl(:,1)/shl_sum -&
                 (tmpshl(:)*shgradl_sum(1))/(shl_sum**2)
  shgradl(:,2) = tmpshgradl(:,2)/shl_sum -&
                 (tmpshl(:)*shgradl_sum(2))/(shl_sum**2)

  shhessl(:,1) = tmpshhessl(:,1)/shl_sum - &
                 tmpshgradl(:,1)*shgradl_sum(1)/(shl_sum**2) -&
                 ((tmpshgradl(:,1)*shgradl_sum(1)+tmpshl(:)*shhessl_sum(1))/&
                 (shl_sum**2)&
                  -2.0d0*tmpshl(:)*shgradl_sum(1)*shgradl_sum(1)/&
                 (shl_sum**3))

  shhessl(:,2) = tmpshhessl(:,2)/shl_sum - &
                 tmpshgradl(:,1)*shgradl_sum(2)/(shl_sum**2) - &
                 ((tmpshgradl(:,2)*shgradl_sum(1)+tmpshl(:)*shhessl_sum(2))/ &
                 (shl_sum**2) &
                  -2.0d0*tmpshl(:)*shgradl_sum(1)*shgradl_sum(2)/ &
                 (shl_sum**3))

  shhessl(:,3) = tmpshhessl(:,3)/shl_sum - &
                 tmpshgradl(:,2)*shgradl_sum(2)/(shl_sum**2) - &
                 ((tmpshgradl(:,2)*shgradl_sum(2)+tmpshl(:)*shhessl_sum(3))/ &
                 (shl_sum**2) &
                 - 2.0d0*tmpshl(:)*shgradl_sum(2)*shgradl_sum(2)/ &
                 (shl_sum**3))

  ! Now we calculate the face Jacobian
  xu = 0.0d0
  xd = 0.0d0
  dxdxi   = 0.0d0
  ddxddxi = 0.0d0
  do i = 1, nshl
    do ii = 1, 3
      xu(ii) = xu(ii) + B_NET_U(lIEN(i),ii)*shl(i)
      xd(ii) = xd(ii) + B_NET_D(lIEN(i),ii)*shl(i)
      dxdxi(ii,:)   = dxdxi(ii,:)   + B_NET_D(lIEN(i),ii)*shgradl(i,:)
      ddxddxi(ii,:) = ddxddxi(ii,:) + B_NET_D(lIEN(i),ii)*shhessl(i,:)
    end do
  end do

  ee = dxdxi(1,1)**2+dxdxi(2,1)**2+dxdxi(3,1)**2
  ff = dxdxi(1,1)*dxdxi(1,2)+dxdxi(2,1)*dxdxi(2,2)+dxdxi(3,1)*dxdxi(3,2)
  gg = dxdxi(1,2)**2+dxdxi(2,2)**2+dxdxi(3,2)**2

  DetJb = sqrt(ee*gg-ff**2) ! Jacobian of face mapping

  ! outward normal vector
  nor(1) = dxdxi(2,1)*dxdxi(3,2) - dxdxi(3,1)*dxdxi(2,2)
  nor(2) = dxdxi(3,1)*dxdxi(1,2) - dxdxi(1,1)*dxdxi(3,2)
  nor(3) = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)

  !------------------------------------------
  ! change from dx/du to dxdxi
  !------------------------------------------
  ! Get knot span sizes
  du = U_KNOT(ni+1)-U_KNOT(ni)
  dv = V_KNOT(nj+1)-V_KNOT(nj)

  dxdxi(:,1) = dxdxi(:,1)*du/2.0d0
  dxdxi(:,2) = dxdxi(:,2)*dv/2.0d0

  ddxddxi(:,1) = ddxddxi(:,1)*du*du/4.0d0
  ddxddxi(:,2) = ddxddxi(:,2)*du*dv/4.0d0
  ddxddxi(:,3) = ddxddxi(:,3)*dv*dv/4.0d0

end subroutine eval_SHAPE_shell





!================================================================
! subroutine to reconstruct the shape functions using 
! bezier extraction operator
!================================================================
subroutine eval_SHAPE_bez_sh(gpu, gpv, shl, shgradl, shhessl, &
                             nor, xu, xd, dxdxi, ddxddxi,  &
                             p, q, nsd, nshl, nshb, &
                             lIEN, nnode,  &
                             B_NET_U, B_NET_D, &
                             DetJb, BezExt)
  implicit none

  integer, intent(in)  :: nshl, nshb, nnode, lIEN(nshl), p, q, nsd
  real(8), intent(in)  :: gpu, gpv, BezExt(nshl,nshb), &
                          B_NET_D(nnode,nsd+1),  &
                          B_NET_U(nnode,nsd+1)
  real(8), intent(out) :: shl(nshl), shgradl(nshl,2), shhessl(nshl,3),  &                          
                          xu(NSD), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3), &
                          nor(NSD), DetJb

  real(8) :: M(3,p+1), N(3,q+1), Bern(nshb), dBern(nshb,2), &
             ddBern(nshb,3), ee, ff, gg, &
             tmpshl(nshl), tmpshgradl(nshl,2), tmpshhessl(nshl,3), &
             shl_sum, shgradl_sum(2), shhessl_sum(3), we
  integer :: i, j, k, ii, ct, kk

  M = 0.0d0
  N = 0.0d0

  shl     = 0.0d0
  shgradl = 0.0d0
  shhessl = 0.0d0

  tmpshl     = 0.0d0
  tmpshgradl = 0.0d0
  tmpshhessl = 0.0d0

  ! M(1,4), M(2,4) and M(3,4) are the Bernstein basis and
  ! its first and second derivatives, respectively, in 1D
  call Bernstein_p3(gpu, M(1,:), M(2,:), M(3,:))
  call Bernstein_p3(gpv, N(1,:), N(2,:), N(3,:))

  ! Compute the tensor product of the 2D Bernstein functions
  ct = 0
  do j = 1, q+1
    do i = 1, p+1
      ct = ct + 1

      Bern(ct)     = M(1,i)*N(1,j)

      dBern(ct,1)  = M(2,i)*N(1,j)
      dBern(ct,2)  = M(1,i)*N(2,j)

      ddBern(ct,1) = M(3,i)*N(1,j)
      ddBern(ct,2) = M(2,i)*N(2,j)
      ddBern(ct,3) = M(1,i)*N(3,j)
    end do
  end do

  ! reconstruct the T-splne functions
  tmpshl = MATMUL(BezExt,Bern)

  tmpshgradl(:,1) = MATMUL(BezExt, dBern(:,1))
  tmpshgradl(:,2) = MATMUL(BezExt, dBern(:,2))

  tmpshhessl(:,1) = MATMUL(BezExt, ddBern(:,1))
  tmpshhessl(:,2) = MATMUL(BezExt, ddBern(:,2))
  tmpshhessl(:,3) = MATMUL(BezExt, ddBern(:,3))

  ! multiply the functions with the weight
  do i = 1, nshl
    we = B_NET_U(lIEN(i),NSD+1)

    tmpshl(i)   = tmpshl(i)*we
    tmpshgradl(i,:) = tmpshgradl(i,:)*we
    tmpshhessl(i,:) = tmpshhessl(i,:)*we
  end do

  shl_sum = sum(tmpshl)

  shgradl_sum(1) = sum(tmpshgradl(:,1))
  shgradl_sum(2) = sum(tmpshgradl(:,2))

  shhessl_sum(1) = sum(tmpshhessl(:,1))
  shhessl_sum(2) = sum(tmpshhessl(:,2))
  shhessl_sum(3) = sum(tmpshhessl(:,3))


  shl = tmpshl/shl_sum

  shgradl(:,1) = tmpshgradl(:,1)/shl_sum - &
                 (tmpshl(:)*shgradl_sum(1))/(shl_sum**2)
  shgradl(:,2) = tmpshgradl(:,2)/shl_sum - &
                 (tmpshl(:)*shgradl_sum(2))/(shl_sum**2)

  shhessl(:,1) = tmpshhessl(:,1)/shl_sum - &
                 tmpshgradl(:,1)*shgradl_sum(1)/(shl_sum**2) - &
                 ((tmpshgradl(:,1)*shgradl_sum(1)+tmpshl(:)*shhessl_sum(1))/ &
                 (shl_sum**2)  &
                 -2d+0*tmpshl(:)*shgradl_sum(1)*shgradl_sum(1)/ &
                 (shl_sum**3))

  shhessl(:,2) = tmpshhessl(:,2)/shl_sum - &
                 tmpshgradl(:,1)*shgradl_sum(2)/(shl_sum**2) - &
                 ((tmpshgradl(:,2)*shgradl_sum(1)+tmpshl(:)*shhessl_sum(2))/ &
                 (shl_sum**2)  &
                 -2d+0*tmpshl(:)*shgradl_sum(1)*shgradl_sum(2)/ &
                 (shl_sum**3))

  shhessl(:,3) = tmpshhessl(:,3)/shl_sum - &
                 tmpshgradl(:,2)*shgradl_sum(2)/(shl_sum**2) - &
                 ((tmpshgradl(:,2)*shgradl_sum(2)+tmpshl(:)*shhessl_sum(3))/ &
                 (shl_sum**2)  &
                 - 2d+0*tmpshl(:)*shgradl_sum(2)*shgradl_sum(2)/ &
                 (shl_sum**3))


  ! Now we calculate the face Jacobian
  xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
  do i = 1, nshl
    do ii = 1, 3
      xu(ii) = xu(ii) + B_NET_U(lIEN(i),ii)*shl(i)
      xd(ii) = xd(ii) + B_NET_D(lIEN(i),ii)*shl(i)
      dxdxi(ii,:)   = dxdxi(ii,:)   + B_NET_D(lIEN(i),ii)*shgradl(i,:)
      ddxddxi(ii,:) = ddxddxi(ii,:) + B_NET_D(lIEN(i),ii)*shhessl(i,:)
    end do
  end do

  ee = dxdxi(1,1)**2+dxdxi(2,1)**2+dxdxi(3,1)**2
  ff = dxdxi(1,1)*dxdxi(1,2)+dxdxi(2,1)*dxdxi(2,2)+dxdxi(3,1)*dxdxi(3,2)
  gg = dxdxi(1,2)**2+dxdxi(2,2)**2+dxdxi(3,2)**2

  DetJb = sqrt(ee*gg-ff**2) ! Jacobian of face mapping

!  nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1)
!  nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1)
!  nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1)

  ! outward normal
  nor(1) = dxdxi(2,1)*dxdxi(3,2) - dxdxi(3,1)*dxdxi(2,2)
  nor(2) = dxdxi(3,1)*dxdxi(1,2) - dxdxi(1,1)*dxdxi(3,2)
  nor(3) = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)

end subroutine eval_SHAPE_bez_sh
  

!================================================================
! Cubic bernstein function
!================================================================
subroutine Bernstein_p3(xi, B, dB, ddB)
  implicit none
  real(8), intent(in)  :: xi
  real(8), intent(out) :: B(4), dB(4), ddB(4)
  real(8) :: x

  x = 0.5d0*(1.0d0+xi)

  B(1) = (1.0d0-x)**3
  B(2) = 3.0d0*x*(1.0d0-x)**2
  B(3) = 3.0d0*(x**2)*(1.0d0-x)
  B(4) = x**3
  
  dB(1) = -3.0d0*(1.0d0-x)**2
  dB(2) = 3.0d0*(1.0d0-x)*(1.0d0-3.0d0*x)
  dB(3) = 3.0d0*x*(2.0d0-3.0d0*x)
  dB(4) = 3.0d0*x**2

  dB = dB*0.5d0

  ddB(1) = 6.0d0*(1.0d0-x)
  ddB(2) = 6.0d0*(3.0d0*x-2.0d0)
  ddB(3) = 6.0d0*(1.0d0-3.0d0*x)
  ddB(4) = 6.0d0*x

  ddB = ddB*0.25d0
end subroutine Bernstein_p3






!================================================================
! subroutine to compute linear triangle shape functions
!================================================================
subroutine eval_SHAPE_tri(gp, shl, shgradl, nor, xu, xd, dxdxi, &
                          nsd, nshl, lIEN, nnode, B_NET_U,   &
                          B_NET_D, DetJb)

  implicit none
      
  integer, intent(in)  :: nshl, nnode, lIEN(nshl), nsd
  real(8), intent(in)  :: gp(2), &
                          B_NET_D(nnode,nsd+1),  &
                          B_NET_U(nnode,nsd+1)
  real(8), intent(out) :: shl(nshl), shgradl(nshl,2), &                          
                          xu(NSD), xd(NSD), dxdxi(NSD,2), &
                          nor(NSD), DetJb

  integer :: i, ii
      
  call lintrishl(gp, NSHL, NSD, shl, shgradl)
      
  ! Now we calculate the face Jacobian
  xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
  do i = 1, nshl
    do ii = 1, 3
      xu(ii) = xu(ii) + B_NET_U(lIEN(i),ii)*shl(i)
      xd(ii) = xd(ii) + B_NET_D(lIEN(i),ii)*shl(i)
      dxdxi(ii,:) = dxdxi(ii,:) + B_NET_D(lIEN(i),ii)*shgradl(i,:)
    end do
  end do

  nor = 0.0d0        
  nor(1) = dxdxi(2,1)*dxdxi(3,2) - dxdxi(3,1)*dxdxi(2,2)
  nor(2) = dxdxi(3,1)*dxdxi(1,2) - dxdxi(1,1)*dxdxi(3,2)
  nor(3) = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)
      
  DetJb  = sqrt(sum(nor*nor))
  nor(:) = nor(:)/DetJb

  if (DetJb < 1.0d-15) then
    write(*,*) "ERROR: Negative Jacobian:", DetJb
    stop
  end if

end subroutine eval_SHAPE_tri



!================================================================
! subroutine to compute linear triangle shape functions
! remove unnecessary computations...
!================================================================
subroutine eval_SHAPE_tri_fast(gp, nsd, nshl, lIEN, nnode, &
                               B_NET, shl, shgradl, xd, dxdxi)
  implicit none
      
  integer, intent(in)  :: nshl, nnode, lIEN(nshl), nsd
  real(8), intent(in)  :: gp(2), B_NET(nnode,nsd+1)
  real(8), intent(out) :: shl(nshl), shgradl(nshl,2), &                          
                          xd(NSD), dxdxi(NSD,2)
  integer :: i, ii
      
  call lintrishl(gp, NSHL, NSD, shl, shgradl)
      
  ! Now we calculate the face Jacobian
  xd = 0.0d0; dxdxi = 0.0d0
  do i = 1, nshl
    do ii = 1, 3
      xd(ii) = xd(ii) + B_NET(lIEN(i),ii)*shl(i)
      dxdxi(ii,:) = dxdxi(ii,:) + B_NET(lIEN(i),ii)*shgradl(i,:)
    end do
  end do
end subroutine eval_SHAPE_tri_fast
