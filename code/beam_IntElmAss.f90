subroutine IntElmAss_beam(NRB, icnt, col, row, nsd, &
                           NMat, Ec, &
                           g_fact,  &
                           RHSG_SH, RHSG_GRA_SH, &
                           LHSK_SH, tg_SH, mg_SH, T_Flag)

  use types_beam
  implicit none

  type(mesh), intent(in) :: NRB

  integer, intent(in) :: icnt, nmat, nsd, T_Flag, &
                         col(NRB%NNODE+1), &
                         row(NRB%NNODE*50*NRB%maxNBEA)
  real(8), intent(in) :: Ec(NMat), g_fact(nsd)

  real(8), intent(inout) :: RHSG_SH(NRB%NNODE,NSD), &
                            RHSG_GRA_SH(NRB%NNODE,NSD),&
                            LHSK_SH(NSD*NSD,icnt),&
                            tg_SH(NRB%NNODE), mg_SH(NRB%NNODE)

  !  Local variables
  integer :: p, nshl, nuk, nvk, ptype, iel, igauss, &
             i, j, ii, jj, kk, ni

  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), gwt, da, VVal, &
             DetJb_SH, nor(NSD), Dm(3,3), Dc(3,3), Db(3,3),&
             xu(NSD), xd(NSD), dxdxi(NSD), ddxddxi(nsd), &
             dens, bvec(NSD), bscale, rtmp1(3), rtmp2

  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:), &
                          shhessl(:), NThk(:), lmass(:), &
                          Rhs(:,:), Rhs_grav(:,:), &
                          xKebe(:,:,:), xKebegp(:,:,:), Rhsgp(:,:)

  
  ! get Gaussian points and weights
  gp = 0.0d0; gw = 0.0d0
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 

  ! loop over elements
  rtmp1 = 0.0d0; rtmp2 = 0.0d0
  do iel = 1, NRB%NEL
    
    ! get NURB coordinates
    ni = NRB%INN(iel)
    
    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't
    if (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) then
         
      ! used in calculating quadrature points. The factor of 4.0d0
      ! comes from mapping from the [-1,1] line onto a real segment...
      da = (NRB%U_KNOT(iel,ni+1)-NRB%U_KNOT(iel,ni))/2.0d0
      
      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NBEA(iel)
      ptype = NRB%PTYPE(iel)
                    
      allocate(shl(nshl), shgradl(nshl), &
               shhessl(nshl), lIEN(nshl), &
               Rhs(NSD,nshl),  &
               Rhs_grav(NSD,nshl), NThk(nshl), lmass(nshl),  &
               xKebe(NSD*NSD,nshl,nshl),  &
               xKebegp(NSD*NSD,nshl,nshl), &
               Rhsgp(NSD,nshl))

      lIEN = -1
      do i = 1, nshl
        lIEN(i) = NRB%IEN(iel,i)
      end do      

      ! initialization      
      xKebe    = 0.0d0      ! initialize local stiffness matrix
      Rhs      = 0.0d0      ! initialize local load vector
      Rhs_grav = 0.0d0
      Nthk     = 0.0d0      ! initialize local nodal thickness
      lmass    = 0.0d0

      ! Loop over integration points (NGAUSS in each direction) 
        do igauss = 1, NRB%NGAUSS
        
          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
          nor = 0.0d0

          call eval_SHAPE_beam(gp(igauss),  &
                                shl, shgradl, shhessl, nor,  &
                                xu, xd, dxdxi, ddxddxi,  &
                                p, nsd, nshl, &
                                lIEN, NRB%NNODE, &
                                NRB%B_NET_U, NRB%B_NET_D, DetJb_SH, &
                                ni, nuk, &
                                NRB%U_KNOT(iel,1:nuk))

          gwt = gw(igauss)*1.0*da


          ! extensional, coupling and bending material matrices
          Dm = 0.0d0; Dc = 0.0d0; Db = 0.0d0
          if (ptype > 0) then

            ! density
            dens = 8.0d3
 
            ! note: ptype needs to match the material type
            Dm = matA(ptype,:,:)*thi
            Dc = matB(ptype,:,:)*thi**2
            Db = matD(ptype,:,:)*thi**3

          else if (ptype == 0) then    ! bending strips

            dens = 0.0d0
            Db(1,1) = 100.0d0*max(matD(1,1,1),matD(1,2,2))*thi**3
            Db(2,2) = Db(1,1)

          else
            write(*,*) "ERROR: UNDEFINED PTYPE"
            stop
          end if
          != end thickness and stiffness ===============        
         
          ! Kirchhoff-Love Shell
          xKebegp = 0.0d0
          Rhsgp   = 0.0d0
          call e3LRhs_KLShell(shgradl, shhessl, Dm, Dc, Db, &
                              xKebegp, Rhsgp, nshl, q, nsd, nor, &
                              NRB%B_NET, NRB%B_NET_D,  &
                              NRB%IEN(iel,1:nshl), &
                              NRB%NNODE)

          xKebe = xKebe + xKebegp*gwt
          Rhs   = Rhs   + Rhsgp*gwt
          ! end Kirchhoff-Love Shell

          ! project the thickness to the nodes
          if (ptype /= 0) then
            do ii = 1, nshl
              NThk(ii)  = NThk(ii)  + thi*shl(ii)*DetJb_SH*gwt
              lmass(ii) = lmass(ii) +     shl(ii)*DetJb_SH*gwt
            end do
          end if

!          ! for computing center of mass
!          rtmp1 = rtmp1 + xu*dens*thi*DetJb_SH*gwt
!          rtmp2 = rtmp2 +    dens*thi*DetJb_SH*gwt
          
          write(*,*) 'Wrong beam_IntMass'

          ! Gravity effect
          bvec = g_fact*dens*9.81d0
          do ii = 1, nshl
            Rhs_grav(:,ii) = Rhs_grav(:,ii) + &
                             shl(ii)*bvec(:)*thi*DetJb_SH*gwt
          end do

      end do  ! end loop gauss points

      call BCLhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, xKebe)
      call BCRhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, Rhs)
      call BCRhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, Rhs_grav)        

      ! Assemble load vector         
      ! Assemble thickness and lump mass
      ! LocaltoGlobal_3D is removed..
      do ii = 1, NBEA
        ! internal 
        RHSG_SH(lIEN(ii),:) = RHSG_SH(lIEN(ii),:) + Rhs(:,ii)
        ! gravity force
        RHSG_GRA_SH(lIEN(ii),:) = RHSG_GRA_SH(lIEN(ii),:) + Rhs_grav(:,ii)
        ! thickness
        tg_SH(lIEN(ii)) = tg_SH(lIEN(ii)) + NThk(ii)
        ! area
        mg_SH(lIEN(ii)) = mg_SH(lIEN(ii)) + lmass(ii)
      end do

      call FillSparseMat_3D_shell(nsd, nshl, lIEN, NRB%NNODE, &
                                  NRB%maxNBEA, icnt, col, row, &
                                  xKebe, LHSK_SH)
  
      deallocate(shl, shgradl, shhessl, Rhs, Rhs_grav,  &
                 NThk, lmass, xKebe, xKebegp, Rhsgp, lIEN)

    end if  ! end if nonzero areas elements

  end do    ! end loop elements
  
!  ! output center of mass
!  write(*,*) rtmp1/rtmp2
!  stop

end subroutine IntElmAss_beam
