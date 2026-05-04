subroutine IntElmAss_beam(NRB, icnt, col, row, nsd, &
                           NMat, Ec, dens, Rad, C_dp, &
                           g_fact, length,  &
                           RHSG_SH, &
                           LHSK_SH, Delta_t, &
                           alpha_m, alpha_f, gamma_t, beta_t)

  use types_beam
  implicit none

  type(mesh_beam), intent(in) :: NRB

  integer, intent(in) :: icnt, nmat, nsd, &
                         col(NRB%NNODE+1), &
                         row(NRB%NNODE*50*NRB%maxNSHL)
  real(8), intent(in) :: Ec(NMat), dens(NMat), Rad(NMat), C_dp(NMat), g_fact(nsd), &
                         Delta_t,alpha_m, alpha_f, gamma_t, beta_t
  real(8), intent(inout) :: RHSG_SH(NRB%NNODE,NSD),LHSK_SH(NSD*NSD,icnt)
  real(8), intent(out) :: length

  !  Local variables
  integer :: p, nshl, nuk, nvk, ptype, iel, igauss, &
             i, j, ii, jj, kk, ni, psd_switch, psd_u_switch

  real(8) :: pi
  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), gwt, da, VVal, &
             DetJb_SH,&
             tgt(NSD), nor(NSD), binor(NSD), gradnor(NSD), gradbinor(NSD), &
             tgt_u(NSD), nor_u(NSD), binor_u(NSD), &
             xu(NSD), xd(NSD), dxdxi(NSD), ddxddxi(nsd), dddxdddxi(nsd), &
             dxdxi_u(NSD), ddxddxi_u(NSD), &
             ddxddxi_psd(nsd), ddxddxi_u_psd(NSD), &
             ddu_alpham(NSD), du_alphaf(NSD), &
             bvec(NSD), bscale, rtmp1(3), rtmp2

  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:), &
                          shhessl(:), shthdgradl(:), NThk(:), lmass(:), &
                          Rhs(:,:), Rhs_grav(:,:), &
                          xKebe(:,:,:), xKebegp(:,:,:), Rhsgp(:,:), Rhs_ext(:,:), Rhs_extgp(:,:)
  
  
  !write(*,*) "G_fact=", G_fact
  pi = acos(-1.0d0)
  ! get Gaussian points and weights
  gp = 0.0d0; gw = 0.0d0
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 
  
  length  = 0.0d0 

  ! loop over elements
  rtmp1 = 0.0d0; rtmp2 = 0.0d0

  do iel = 1, NRB%NEL
    
    !write(*,*) "iel=", iel

    ! get NURB coordinates
    ni = NRB%INN(iel)
    
    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't

    if (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) then
         
      ! used in calculating quadrature points. The factor of 4.0d0
      ! comes from mapping from the [-1,1] line onto a real segment...
      da = (NRB%U_KNOT(iel,ni+1)-NRB%U_KNOT(iel,ni))/2.0d0
      
      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      ptype = NRB%PTYPE(iel)
                    
      allocate(shl(nshl), shgradl(nshl), &
               shhessl(nshl), shthdgradl(nshl), lIEN(nshl), &
               Rhs(NSD,nshl),  &
               Rhs_grav(NSD,nshl), NThk(nshl), lmass(nshl),  &
               xKebe(NSD*NSD,nshl,nshl),  &
               xKebegp(NSD*NSD,nshl,nshl), &
               Rhsgp(NSD,nshl), Rhs_ext(NSD,nshl), Rhs_extgp(NSD,nshl))

      lIEN = -1
      do i = 1, nshl
        lIEN(i) = NRB%IEN(iel,i)
      end do      

      ! initialization      
      xKebe    = 0.0d0      ! initialize local stiffness matrix
      Rhs      = 0.0d0      ! initialize local load vector
      Rhs_ext  = 0.0d0
      Rhs_grav = 0.0d0
      Nthk     = 0.0d0      ! initialize local nodal thickness
      lmass    = 0.0d0

      ! Loop over integration points (NGAUSS in each direction) 
        do igauss = 1, NRB%NGAUSS
        
          !write(*,*) "igauss=", igauss
          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0; shthdgradl= 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0; dddxdddxi = 0.0d0
          tgt= 0.0d0 ;nor = 0.0d0; binor= 0.0d0
          tgt_u= 0.0d0 ;nor_u = 0.0d0; binor_u= 0.0d0
          dxdxi_u= 0.0d0; ddxddxi_u= 0.0d0; tgt_u= 0.0d0
          ddu_alpham= 0.0d0; du_alphaf= 0.0d0
          ddxddxi_psd=0.0d0 ; ddxddxi_u_psd= 0.0d0
          ! pseudo vector mode is by default open
          psd_switch = 1; psd_u_switch = 1

          call eval_SHAPE_beam(gp(igauss),  &
                                shl, shgradl, shhessl,shthdgradl, &
                                tgt, nor, binor, gradnor, gradbinor, &
                                tgt_u, nor_u, binor_u, &
                                xu, xd, dxdxi, ddxddxi,dddxdddxi, &
                                dxdxi_u, ddxddxi_u, &
                                ddu_alpham, du_alphaf, &
                                ddxddxi_psd, ddxddxi_u_psd, psd_switch, psd_u_switch,&
                                p, nsd, nshl, &
                                lIEN, NRB%NNODE, &
                                NRB%B_NET_U, NRB%B_NET_D_alphaf, &
                                NRB%B_NET_Dt_alphaf, NRB%B_NET_DDt_alpham, DetJb_SH, &
                                ni, nuk, &
                                NRB%U_KNOT(iel,1:nuk))
          !write(*,*) "xu=", xu          
          !write(*,*) "xd=", xd
          !write(*,*) "dxdxi=", dxdxi
          !write(*,*) "dxds=", dxdxi*40d0/3.1416d0
          !write(*,*) "|dxds|=", sqrt(sum(dxdxi**2))*40d0/3.1416d0
          !write(*,*) "ddxddxi=", ddxddxi
          !write(*,*) "dddxdddxi=", dddxdddxi
          !write(*,*) "tgt=", tgt
          !write(*,*) "nor=", nor
          !write(*,*) "binor=", binor
          !write(*,*) "gradnor=", gradnor
          !write(*,*) "gradbinor=", gradbinor
          !write(*,*) "shgradl=", shgradl
          !write(*,*) "shhessl=", shhessl
          gwt = gw(igauss)*da
          length = length + DetJb_SH*gwt

          !dens=8.0d3
          ! Kirchhoff Beam
          xKebegp = 0.0d0
          Rhsgp   = 0.0d0
          Rhs_extgp = 0.0d0
          call e3LRhs_KLBeam(shl, shgradl, shhessl, ptype, Ec(ptype), Rad(ptype), Dens(ptype), C_dp(ptype), &
                              xKebegp, Rhsgp, Rhs_extgp, nshl, nsd, &
                              tgt, nor, binor, &
                              tgt_u, nor_u, binor_u, &
                              xu, xd, dxdxi, ddxddxi, &
                              dxdxi_u, ddxddxi_u, &
                              ddu_alpham, du_alphaf, &
                              ddxddxi_psd, ddxddxi_u_psd, psd_switch, psd_u_switch,&
                              NRB%B_NET, NRB%B_NET_D_alphaf, &
                              NRB%IEN(iel,1:nshl), NRB%NNODE, Delta_t, &
                              alpha_m, alpha_f, gamma_t, beta_t)

          !do ii = 1, nshl
          !  do jj = 1, nshl
          !    do i = 1, nsd*nsd
          !      write(*,*) ii,jj,i,'xKebe=', xKebegp(i,ii,jj)
          !      !write(*,'(/)')
          !    enddo   
          !  enddo
          !enddo


          xKebe = xKebe + xKebegp*gwt*DetJb_SH
          Rhs   = Rhs   + Rhsgp*gwt*DetJb_SH

          if ((ptype .eq. 1) .or. (ptype .eq. 3) ) then

          ! Gravity effect
          bvec = g_fact*dens(ptype)*9.81d0
          do ii = 1, nshl
            Rhs(:,ii) = Rhs(:,ii) + &
                             shl(ii)*bvec(:)*pi*(Rad(ptype)**2.0d0)*DetJb_SH*gwt
          end do
          
          endif


      end do  ! end loop gauss points

      !do ii = 1, nshl
      !  do jj = 1, nshl
      !    do i = 1, nsd*nsd
      !      write(*,*) ii,jj,i,'xKebe=', xKebe(i,ii,jj)
      !      !write(*,'(/)')
      !    enddo   
      !  enddo
      !enddo

      call BCLhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, xKebe)
      call BCRhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, Rhs)

      !do ii = 1, nshl
      !  do jj = 1, nshl
      !    do i = 1, nsd*nsd
      !      write(*,*) ii,jj,i,'xKebe=', xKebe(i,ii,jj)
      !      !write(*,'(/)')
      !    enddo   
      !  enddo
      !enddo

      !write(*,'(///)')

      ! Assemble load vector         
      ! Assemble thickness and lump mass
      ! LocaltoGlobal_3D is removed..
      do ii = 1, nshl
        ! internal 
        RHSG_SH(lIEN(ii),:) = RHSG_SH(lIEN(ii),:) + Rhs(:,ii)
      end do


      call FillSparseMat_3D_shell(nsd, nshl, lIEN, NRB%NNODE, &
                                  NRB%maxNSHL, icnt, col, row, &
                                  xKebe, LHSK_SH)

      deallocate(shl, shgradl, shhessl, shthdgradl, Rhs, Rhs_grav,  &
                 NThk, lmass, xKebe, xKebegp, Rhsgp, Rhs_ext, Rhs_extgp, lIEN)



    end if  ! end if nonzero areas elements

  end do    ! end loop elements
  write(*,*) "length=",length
!  ! output center of mass
!  write(*,*) rtmp1/rtmp2
!  stop

end subroutine IntElmAss_beam
