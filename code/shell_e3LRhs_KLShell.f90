subroutine e3LRhs_KLShell(shgradl, shhessl, Dm, Dc, Db, &
                          xKebe, Rhs, nshl, q, NSD,  &
                          B_NET_SH, B_NET_SH_D, lIEN, nnode)

  implicit none  

  integer, intent(in) :: nshl, q, NSD, lIEN(nshl), nnode    
  real(8), intent(in) :: shgradl(NSHL,2), shhessl(NSHL,3), &
                         Dm(3,3), Dc(3,3), Db(3,3), &
                         B_NET_SH(nnode,nsd+1), &
                         B_NET_SH_D(nnode,nsd+1)
  real(8), intent(out):: xKebe(NSD*NSD,NSHL,NSHL), Rhs(NSD,NSHL)

  integer :: ur, kr, dirr, us, ks, dirs, dirt, ddir, i
  real(8) :: dR(NSHL,2), ddR(NSHL,3)
  real(8) :: Gab_r(2,2), Bv_r(3), Tb(3,3), dA
  real(8) :: g(3,2), g3(3), lg3, n(3), gab(2,2), h(3,3), bv(3)
  real(8) :: lg3_3, lg3_5
  real(8) :: E_cu(3), E_ca(3), K_cu(3), K_ca(3)
  real(8) :: dg(3,2,3*NSHL), dE_cu(3), dE_ca(3,3*NSHL), &
             dg3(3,3*NSHL), g3dg3(3*NSHL), g3dg3lg3_3(3*NSHL), &
             dn(3,3*NSHL), dK_cu(3), dK_ca(3,3*NSHL)
  real(8) :: ddE_cu(3), ddE_ca(3,3*NSHL,3*NSHL), &
             ddg3(3), tmp1, tmp2, ddn(3), ddK_cu(3), &
             ddK_ca(3,3*NSHL,3*NSHL), N_ca(3), M_ca(3), &
             dN_ca(3), dM_ca(3)
  real(8) :: kem(3*NSHL,3*NSHL), keb(3*NSHL,3*NSHL), &
             ke(3*NSHL,3*NSHL)
  real(8) :: fiem(3*NSHL), fieb(3*NSHL), fie(3*NSHL)

  real(8) :: tmp33(3,3), tmp32(3,2), tmp3(3), tmp22(2,2)
  
  kem = 0d0
  keb = 0d0
  ke =  0d0

  fiem = 0d0
  fieb = 0d0
  fie  = 0d0

  do i = 1, NSHL
    dR(i,:) = shgradl(i,:)
    ddR(i,1) = shhessl(i,1)
    ddR(i,2) = shhessl(i,3)
    ddR(i,3) = shhessl(i,2)
  end do
  
  ! reference configuration
  call shell_geo(nshl, nsd, nnode, q, dR, ddR, lIEN, B_NET_SH, &
                 tmp32, tmp3, dA, tmp3, Gab_r, tmp22, tmp33, &
                 Bv_r, Tb, tmp33, tmp33)

  ! actual configuration
  call shell_geo(nshl, nsd, nnode, q, dR, ddR, lIEN, B_NET_SH_D, &
                 g, g3, lg3, n, gab, tmp22, h, bv, tmp33, tmp33, &
                 tmp33)

  lg3_3 = lg3**3
  lg3_5 = lg3**5

  ! strain vector [E11,E22,E12] referred to curvilinear coor sys
  E_cu(1) = 0.5d0*(gab(1,1)-Gab_r(1,1))
  E_cu(2) = 0.5d0*(gab(2,2)-Gab_r(2,2))
  E_cu(3) = 0.5d0*(gab(1,2)-Gab_r(1,2))

  ! strain vector [E11,E22,2*E12] referred to cartesian coor sys
  ! E_ca = Tb*E_cu
  E_ca = MATMUL(Tb,E_cu)

  ! curvature vector [K11,K22,K12] referred to curvilinear coor sys
  K_cu = bv - Bv_r

  ! curvature vector [K11,K22,2*K12] referred to cartesian coor sys
  ! K_ca = Tb*K_cu
  K_ca = MATMUL(Tb,K_cu)

  dg = 0.0d0
  ! first variation of strain and curvature w.r.t. dof
  do ur = 1, 3*NSHL
    ! local node number kr and dof direction dirr
    kr = (ur+2)/3
    dirr = ur-3*(kr-1)

    dg(dirr,1,ur) = dR(kr,1)
    dg(dirr,2,ur) = dR(kr,2)
  
    ! strain
    dE_cu(1) = dR(kr,1)*g(dirr,1)
    dE_cu(2) = dR(kr,2)*g(dirr,2)
    dE_cu(3) = 0.5d0*(dR(kr,1)*g(dirr,2) + g(dirr,1)*dR(kr,2))

    dE_ca(1,ur) = Tb(1,1)*dE_cu(1)+Tb(1,2)*dE_cu(2)+Tb(1,3)*dE_cu(3)
    dE_ca(2,ur) = Tb(2,1)*dE_cu(1)+Tb(2,2)*dE_cu(2)+Tb(2,3)*dE_cu(3)
    dE_ca(3,ur) = Tb(3,1)*dE_cu(1)+Tb(3,2)*dE_cu(2)+Tb(3,3)*dE_cu(3)
  
    ! curvature
    dg3(1,ur) = dg(2,1,ur)*g(3,2)-dg(3,1,ur)*g(2,2) &
              + g(2,1)*dg(3,2,ur)-g(3,1)*dg(2,2,ur)
    dg3(2,ur) = dg(3,1,ur)*g(1,2)-dg(1,1,ur)*g(3,2) &
              + g(3,1)*dg(1,2,ur)-g(1,1)*dg(3,2,ur)
    dg3(3,ur) = dg(1,1,ur)*g(2,2)-dg(2,1,ur)*g(1,2) &
              + g(1,1)*dg(2,2,ur)-g(2,1)*dg(1,2,ur)
    
    g3dg3(ur) = g3(1)*dg3(1,ur)+g3(2)*dg3(2,ur)+g3(3)*dg3(3,ur)
    g3dg3lg3_3(ur) = g3dg3(ur)/lg3_3
  
    dn(1,ur) = dg3(1,ur)/lg3 - g3(1)*g3dg3lg3_3(ur)
    dn(2,ur) = dg3(2,ur)/lg3 - g3(2)*g3dg3lg3_3(ur)
    dn(3,ur) = dg3(3,ur)/lg3 - g3(3)*g3dg3lg3_3(ur)
    
    dK_cu(1) = ddR(kr,1)*n(dirr) &
             + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur)+h(3,1)*dn(3,ur)
    dK_cu(2) = ddR(kr,2)*n(dirr) &
             + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur)+h(3,2)*dn(3,ur)
    dK_cu(3) = ddR(kr,3)*n(dirr) &
             + h(1,3)*dn(1,ur)+h(2,3)*dn(2,ur)+h(3,3)*dn(3,ur)
  
    dK_ca(1,ur) = Tb(1,1)*dK_cu(1)+Tb(1,2)*dK_cu(2)+Tb(1,3)*dK_cu(3)
    dK_ca(2,ur) = Tb(2,1)*dK_cu(1)+Tb(2,2)*dK_cu(2)+Tb(2,3)*dK_cu(3)
    dK_ca(3,ur) = Tb(3,1)*dK_cu(1)+Tb(3,2)*dK_cu(2)+Tb(3,3)*dK_cu(3)
  end do

  ! second variation of strain and curvature w.r.t. dofs
  do ur = 1, 3*NSHL
    kr = (ur+2)/3
    dirr = ur-3*(kr-1)
    do us = 1, ur
      ks = (us+2)/3
      dirs = us-3*(ks-1)
      ! strain
      ddE_cu = 0.0d0
      if (dirr == dirs) then
        ddE_cu(1) = dR(kr,1)*dR(ks,1)
        ddE_cu(2) = dR(kr,2)*dR(ks,2)
        ddE_cu(3) = 0.5d0*(dR(kr,1)*dR(ks,2)+dR(kr,2)*dR(ks,1))
      end if
      
      ddE_ca(1,ur,us) = Tb(1,1)*ddE_cu(1)+Tb(1,2)*ddE_cu(2) &
                      + Tb(1,3)*ddE_cu(3)
      ddE_ca(2,ur,us) = Tb(2,1)*ddE_cu(1)+Tb(2,2)*ddE_cu(2) &
                      + Tb(2,3)*ddE_cu(3)
      ddE_ca(3,ur,us) = Tb(3,1)*ddE_cu(1)+Tb(3,2)*ddE_cu(2) &
                      + Tb(3,3)*ddE_cu(3)
      
      ! curvature
      ddg3 = 0.0d0

      dirt = 6-dirr-dirs
      ddir = dirr-dirs      

      if     (ddir==-1) then
        ddg3(dirt) =  dR(kr,1)*dR(ks,2)-dR(ks,1)*dR(kr,2)
      else if (ddir== 2) then
        ddg3(dirt) =  dR(kr,1)*dR(ks,2)-dR(ks,1)*dR(kr,2)
      else if (ddir== 1) then
        ddg3(dirt) = -dR(kr,1)*dR(ks,2)+dR(ks,1)*dR(kr,2)
      else if (ddir==-2) then
        ddg3(dirt) = -dR(kr,1)*dR(ks,2)+dR(ks,1)*dR(kr,2)
      end if
    
      tmp1 = -( ddg3(1)*g3(1) + ddg3(2)*g3(2) + ddg3(3)*g3(3) &
             + dg3(1,ur)*dg3(1,us) + dg3(2,ur)*dg3(2,us) &
             + dg3(3,ur)*dg3(3,us))/lg3_3
      tmp2 = 3.0d0*g3dg3(ur)*g3dg3(us)/lg3_5
    
      ddn(1) = ddg3(1)/lg3 - g3dg3lg3_3(us)*dg3(1,ur) &
             - g3dg3lg3_3(ur)*dg3(1,us) + tmp1*g3(1) + tmp2*g3(1)
      ddn(2) = ddg3(2)/lg3 - g3dg3lg3_3(us)*dg3(2,ur) &
             - g3dg3lg3_3(ur)*dg3(2,us) + tmp1*g3(2) + tmp2*g3(2)
      ddn(3) = ddg3(3)/lg3 - g3dg3lg3_3(us)*dg3(3,ur) &
             - g3dg3lg3_3(ur)*dg3(3,us) + tmp1*g3(3) + tmp2*g3(3)
    
      ddK_cu(1) = ddR(kr,1)*dn(dirr,us) + ddR(ks,1)*dn(dirs,ur) &
                + h(1,1)*ddn(1) + h(2,1)*ddn(2) + h(3,1)*ddn(3)
      ddK_cu(2) = ddR(kr,2)*dn(dirr,us) + ddR(ks,2)*dn(dirs,ur) &
                + h(1,2)*ddn(1) + h(2,2)*ddn(2) + h(3,2)*ddn(3)
      ddK_cu(3) = ddR(kr,3)*dn(dirr,us) + ddR(ks,3)*dn(dirs,ur) &
                + h(1,3)*ddn(1) + h(2,3)*ddn(2) + h(3,3)*ddn(3)
      
      ddK_ca(1,ur,us) = Tb(1,1)*ddK_cu(1) + Tb(1,2)*ddK_cu(2) &
                      + Tb(1,3)*ddK_cu(3)
      ddK_ca(2,ur,us) = Tb(2,1)*ddK_cu(1) + Tb(2,2)*ddK_cu(2) &
                      + Tb(2,3)*ddK_cu(3)
      ddK_ca(3,ur,us) = Tb(3,1)*ddK_cu(1) + Tb(3,2)*ddK_cu(2) &
                      + Tb(3,3)*ddK_cu(3)
    end do
  end do
  
  ! N_ca = Dm*E_ca + Dc*K_ca
  N_ca = MATMUL(Dm,E_ca) + MATMUL(Dc,K_ca)

  ! M_ca = Dc*E_ca + Db*K_ca
  M_ca = MATMUL(Dc,E_ca) + MATMUL(Db,K_ca)
  
  kem = 0.0d0
  keb = 0.0d0
  ! loop over dofs ur and us
  do ur = 1,3*NSHL
    ! dN_ca = Dm*dE_ca(:,ur)
    dN_ca = MATMUL(Dm,dE_ca(:,ur)) + MATMUL(Dc,dK_ca(:,ur))

    ! dM_ca = Db*dK_ca(:,ur)
    dM_ca = MATMUL(Dc,dE_ca(:,ur)) + MATMUL(Db,dK_ca(:,ur))

    do us = 1, ur
      ! membrane stiffness
      kem(ur,us) = dN_ca(1)*dE_ca(1,us) + dN_ca(2)*dE_ca(2,us) &
                 + dN_ca(3)*dE_ca(3,us) &
                 + N_ca(1)*ddE_ca(1,ur,us) + N_ca(2)*ddE_ca(2,ur,us) &
                 + N_ca(3)*ddE_ca(3,ur,us)
      ! bending stiffness
      keb(ur,us) = dM_ca(1)*dK_ca(1,us) + dM_ca(2)*dK_ca(2,us) &
                 + dM_ca(3)*dK_ca(3,us)  &
                 + M_ca(1)*ddK_ca(1,ur,us) + M_ca(2)*ddK_ca(2,ur,us) &
                 + M_ca(3)*ddK_ca(3,ur,us)
    end do
    ! residual
    fiem(ur) = -(N_ca(1)*dE_ca(1,ur) + N_ca(2)*dE_ca(2,ur) &
               + N_ca(3)*dE_ca(3,ur))
    fieb(ur) = -(M_ca(1)*dK_ca(1,ur) + M_ca(2)*dK_ca(2,ur) &
               + M_ca(3)*dK_ca(3,ur))
  end do
  ke  = (kem+keb)*dA
  fie = (fiem+fieb)*dA

  ! transform to Yuri's format  
  do ur = 1, 3*NSHL
    do us = 1, 3*NSHL
      if (us > ur) then
        ke(ur,us) = ke(us,ur)
      end if
    end do
  end do

  do ur = 1,3*NSHL
    kr = (ur+2)/3
    dirr = ur-3*(kr-1)
    do us = 1,3*NSHL
      ks = (us+2)/3
      dirs = us-3*(ks-1)
      i = (dirr-1)*3 + dirs

      xKebe(i,kr,ks) = ke(ur,us)

    end do
    Rhs(dirr,kr) = fie(ur)
  end do
      
end subroutine e3LRhs_KLShell
