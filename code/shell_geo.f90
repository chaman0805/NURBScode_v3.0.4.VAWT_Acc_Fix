subroutine shell_geo(nshl, nsd, nnode, q, dR, ddR, lIEN, B_NET_SH, &
                     g, g3, dA, n, Gab, Gab_con, H, Bv,          &
                     T_Gcon_E, T_E_G, T_G_E)
  implicit none      
  integer, intent(in) :: nshl, nnode, nsd, lIEN(nshl), q
  real(8), intent(in) :: dR(NSHL,2), ddR(NSHL,3), B_NET_SH(nnode,nsd+1)

  real(8), intent(out) :: Bv(3), T_Gcon_E(3,3), dA,    & 
                          T_E_G(3,3), T_G_E(3,3),      &
                          G(3,2), G3(3), N(3), H(3,3), &
                          Gab(2,2), Gab_con(2,2)

  real(8) :: G_con(3,2), E(3,2), invdetGab, lg1, lg_con2, EG(2,2)
  integer :: k, c, b, a, ii
      
  ! covariant base vectors G and hessian H
  G = 0.0d0
  H = 0.0d0
  do ii = 1, nshl
    k = lIEN(ii)
    do a = 1,3
      G(a,1:2) =  dR(ii,1:2)*B_NET_SH(k,a) + G(a,1:2)
      H(a,1:3) = ddR(ii,1:3)*B_NET_SH(k,a) + H(a,1:3)
    end do
  end do
      
      ! basis vector G3
  G3(1) = G(2,1)*G(3,2) - G(3,1)*G(2,2)
  G3(2) = G(3,1)*G(1,2) - G(1,1)*G(3,2)
  G3(3) = G(1,1)*G(2,2) - G(2,1)*G(1,2)

  ! differential area dA = length of G3
  dA = sqrt(sum(G3*G3))

  ! normal vector N
  N = G3/dA
      
  ! curvature coefficients as vector
  Bv(:) = H(1,:)*N(1) + H(2,:)*N(2) + H(3,:)*N(3)

  ! covariant metric Gab
  Gab = matmul(transpose(G),G)

  ! contravariant metric Gab_con and base vectors G_con
  invdetGab = 1.0d0/(Gab(1,1)*Gab(2,2)-Gab(1,2)*Gab(1,2))
  Gab_con(1,1) =  invdetgab*Gab(2,2)
  Gab_con(1,2) = -invdetgab*Gab(1,2)
  Gab_con(2,2) =  invdetgab*Gab(1,1)
  Gab_con(2,1) =  Gab_con(1,2)

  G_con = matmul(G,transpose(Gab_con))

  ! local cartesian coordinates
  lg1 = sqrt(sum(G(:,1)**2))
  E(:,1) = G(:,1)/lg1

  lg_con2 = sqrt(sum(G_con(:,2)**2))
  E(:,2) = G_con(:,2)/lg_con2

  ! quick fix for bending strip
  ! theta_1 need to coincide with the linear direction
  if (q == 1) then
    ! local cartesian coordinates
    lg1 = sqrt(sum(G(:,2)**2))
    E(:,1) = G(:,2)/lg1
    lg_con2 = sqrt(sum(G_con(:,1)**2))
    E(:,2) = G_con(:,1)/lg_con2
  end if
      

  ! Transformation matrix T_Gcon_E from G_con to E
  ! with *2 in last row (for strain in Voigt notation)
  EG = matmul(transpose(E),G_con)

  T_Gcon_E(1,1) = EG(1,1)**2
  T_Gcon_E(1,2) = EG(1,2)**2
  T_Gcon_E(1,3) = 2.0d0*EG(1,1)*EG(1,2)
  T_Gcon_E(2,1) = EG(2,1)**2
  T_Gcon_E(2,2) = EG(2,2)**2
  T_Gcon_E(2,3) = 2.0d0*EG(2,1)*EG(2,2)
  T_Gcon_E(3,1) = 2.0d0*EG(1,1)*EG(2,1)
  T_Gcon_E(3,2) = 2.0d0*EG(1,2)*EG(2,2)
  T_Gcon_E(3,3) = 2.0d0*EG(1,1)*EG(2,2)+EG(1,2)*EG(2,1)

  ! Transformation matrix T_E_G from E to G (for PK2 stress)
  T_E_G(1,1) = EG(1,1)**2
  T_E_G(1,2) = EG(2,1)**2             ! be carefure the difference
  T_E_G(1,3) = 2.0d0*EG(1,1)*EG(2,1)  ! be carefure the difference
  T_E_G(2,1) = EG(1,2)**2
  T_E_G(2,2) = EG(2,2)**2
  T_E_G(2,3) = 2.0d0*EG(1,2)*EG(2,2)
  T_E_G(3,1) = EG(1,1)*EG(1,2)
  T_E_G(3,2) = EG(2,1)*EG(2,2)
  T_E_G(3,3) = EG(1,1)*EG(2,2)+EG(2,1)*EG(1,2)

  ! Transformation matrix T_g_e from g to e (for Cauchy stress)
  EG = matmul(transpose(E),G)

  T_G_E(1,1) = EG(1,1)**2
  T_G_E(1,2) = EG(1,2)**2
  T_G_E(1,3) = 2.0d0*EG(1,1)*EG(1,2)
  T_G_E(2,1) = EG(2,1)**2
  T_G_E(2,2) = EG(2,2)**2
  T_G_E(2,3) = 2.0d0*EG(2,1)*EG(2,2)
  T_G_E(3,1) = EG(1,1)*EG(2,1)
  T_G_E(3,2) = EG(1,2)*EG(2,2)
  T_G_E(3,3) = EG(1,1)*EG(2,2)+EG(1,2)*EG(2,1)

end subroutine shell_geo
