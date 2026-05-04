subroutine IntElmAss_shell(SH, NRB, icnt, col, row, nsd, &
                            g_fact, &
                            RHSG_SH, &
                            LHSK_SH, ashAlpha, &
                            alfi, beti, almi, Delt)
  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(in) :: NRB
  type(shell_bld), intent(in) :: SH

  integer, intent(in) :: icnt, nsd, &
                         col(NRB%NNODE+1), &
                         row(NRB%NNODE*50*NRB%maxNSHL)

  real(8), intent(in) :: g_fact(nsd), &
                         ashAlpha(NRB%NNODE,NSD), &
                         alfi, beti, almi, Delt

  real(8), intent(inout) :: RHSG_SH(NRB%NNODE,NSD), LHSK_SH(NSD*NSD,icnt)

  !  Local variables
  integer :: p, q, nshl, nuk, nvk, ptype, iel, igauss, jgauss, &
             i, j, ii, jj, kk, ni, nj, ct, aa, bb

  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), gwt, da, VVal, &
             DetJb_SH, nor(NSD), thi, Dm(3,3), Dc(3,3), Db(3,3),&
             xu(NSD), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3), &
             dens, bvec(NSD), bscale

  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:), &
                          Rhs(:,:), Rhs_gra(:,:), Rhsgp(:,:), &
                          xMebe(:,:), xKebe(:,:,:), xKebegp(:,:,:), &
                          acl(:,:)

  gp = 0.0d0; gw = 0.0d0  
  DetJb_SH = 0.0d0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 

  ! loop over elements
  do iel = 1, NRB%NEL
  
    ! get NURB coordinates
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)
    
    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't
    if ((NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
        (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1))   ) then
         
      ! used in calculating quadrature points. The factor of 4.0d0
      ! comes from mapping from the [-1,1] line onto a real segment...
      da = (NRB%U_KNOT(iel,ni+1)-NRB%U_KNOT(iel,ni))*  &
           (NRB%V_KNOT(iel,nj+1)-NRB%V_KNOT(iel,nj))/4.0d0

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel); ptype = NRB%PTYPE(iel)

    allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), &
             Rhs(NSD,nshl), Rhs_gra(NSD,nshl), &
             xMebe(nshl,nshl),  &
             xKebe(NSD*NSD,nshl,nshl),  &
             xKebegp(NSD*NSD,nshl,nshl), &
             Rhsgp(NSD,nshl), lIEN(nshl), acl(nshl,nsd))

    lIEN = -1
    do i = 1, nshl
      lIEN(i) = NRB%IEN(iel,i)
    end do

    ! Get local solution arrays
    do i = 1, nshl
      acl(i,:) = ashAlpha(NRB%IEN(iel,i),:)
    end do

    ! initialization   
    xMebe   = 0.0d0     
    xKebe   = 0.0d0      ! initialize local stiffness matrix
    Rhs     = 0.0d0      ! initialize local load vector
    Rhs_gra = 0.0d0

    ! Loop over integration points (NGAUSS in each direction)
    ct = 0 
    do jgauss = 1, NRB%NGAUSS
      do igauss = 1, NRB%NGAUSS

        ct = ct + 1
        
          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
          nor = 0.0d0

          call eval_SHAPE_shell(gp(igauss), gp(jgauss),  &
                                shl, shgradl, shhessl, nor,  &
                                xu, xd, dxdxi, ddxddxi,  &
                                p, q, nsd, nshl, &
                                lIEN, NRB%NNODE, &
                                NRB%B_NET_U, NRB%B_NET_D, DetJb_SH, &
                                ni, nj, nuk, nvk, &
                                NRB%U_KNOT(iel,1:nuk), &
                                NRB%V_KNOT(iel,1:nvk))

          gwt = gw(igauss)*gw(jgauss)*da

          !Total thickness of the section-laminate:
          thi = SH%Thickness(iel,ct)

          ! extensional, coupling and bending material matrices
          Dm = 0.0d0; Dc = 0.0d0; Db = 0.0d0
          if (ptype > 0) then

            dens = SH%Density(iel,ct)

            do i = 1, 3
              do j = 1, 3
                Dm(i,j) = SH%matA(iel,ct,i,j)
                Dc(i,j) = SH%matB(iel,ct,i,j)
                Db(i,j) = SH%matD(iel,ct,i,j)
              end do
            end do

          else if (ptype == 0) then    ! bending strips

            dens = 0.0d0
            Db(1,1) = SH%matD(iel,ct,1,1)/10.0d0
            Db(2,2) = SH%matD(iel,ct,2,2)/10.0d0

          else
            write(*,*) "ERROR: UNDEFINED PTYPE"
            stop
          end if
        != end thickness and stiffness ===============
                
        ! Kirchhoff-Love Shell by J. Kiendl
        xKebegp = 0.0d0
        Rhsgp   = 0.0d0
        call e3LRhs_KLShell(shgradl, shhessl, Dm, Dc, Db, &
                            xKebegp, Rhsgp,  &
                            nshl, q, nsd, &
                            NRB%B_NET, NRB%B_NET_D,  &
                            NRB%IEN(iel,1:nshl), NRB%NNODE)
                            
        xKebe = xKebe + xKebegp*gwt
        Rhs   = Rhs   + Rhsgp*gwt
        ! end Kirchhoff-Love Shell

        ! Build Mass
        do aa = 1, nshl
          do bb = 1, nshl
            xMebe(aa,bb) = xMebe(aa,bb) + &
                           thi*dens*shl(aa)*shl(bb)*DetJb_SH*gwt
          end do
        end do

        ! Gravity effect
        bvec = g_fact*dens*9.81d0
        do aa = 1, nshl
          Rhs(:,aa) = Rhs(:,aa) + &
                          shl(aa)*bvec(:)*thi*DetJb_SH*gwt
        end do

      end do
    end do  ! end loop gauss points

    ! Dynamic part (RHS)
    do aa = 1, nshl
      do bb = 1, nshl
        Rhs(1,aa) = Rhs(1,aa) - xMebe(aa,bb)*acl(bb,1)
        Rhs(2,aa) = Rhs(2,aa) - xMebe(aa,bb)*acl(bb,2)
        Rhs(3,aa) = Rhs(3,aa) - xMebe(aa,bb)*acl(bb,3)
      end do
    end do

    ! Dynamic part (LHS)
    xKebe = alfi*beti*Delt*Delt*xKebe
    xKebe(1,:,:) = xKebe(1,:,:) + almi*xMebe(:,:)
    xKebe(5,:,:) = xKebe(5,:,:) + almi*xMebe(:,:)
    xKebe(9,:,:) = xKebe(9,:,:) + almi*xMebe(:,:)

    call BCLhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, xKebe)
    call BCRhs_3D_shell(nsd, nshl, NRB%NNODE, lIEN, NRB%IBC, Rhs)

    ! Assemble load vector         
    ! Assemble thickness and lump mass
    ! LocaltoGlobal_3D is removed..
    do aa = 1, NSHL
      ! internal 
      RHSG_SH(lIEN(aa),:) = RHSG_SH(lIEN(aa),:) + Rhs(:,aa)

    end do


    call FillSparseMat_3D_shell(nsd, nshl, lIEN, NRB%NNODE, &
                                NRB%maxNSHL, icnt, col, row, &
                                xKebe, LHSK_SH)

    deallocate(shl, shgradl, shhessl, Rhs, Rhs_gra,  &
               xMebe, xKebe, xKebegp, Rhsgp, lIEN, acl)

    end if  ! end if nonzero areas elements

  end do    ! end loop elements
end subroutine IntElmAss_shell




!=======================================================================
! Element loop for T-Spline
!=======================================================================
subroutine IntElmAss_tsp_sh(TSP, BEZ, icnt, col, row, nsd, &
                            NMat, MatA, MatB, MatD, g_fact, &
                            RHSG_SH, RHSG_GRA_SH, &
                            LHSK_SH, T_Flag, ashAlpha, &
                            alfi, beti, almi, Delt, BldRot)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ

  integer, intent(in) :: icnt, nmat, nsd, T_Flag, &
                         col(TSP%NNODE+1), &
                         row(TSP%NNODE*50*TSP%maxNSHL)

  real(8), intent(in) :: MatA(NMat,3,3), MatB(NMat,3,3), &
                         MatD(NMat,3,3), g_fact(nsd), &
                         ashAlpha(TSP%NNODE,NSD), &
                         alfi, beti, almi, Delt, BldRot

  real(8), intent(out) :: RHSG_SH(TSP%NNODE,NSD), &
                          RHSG_GRA_SH(TSP%NNODE,NSD), &
                          LHSK_SH(NSD*NSD,icnt)

  !  Local variables
  integer :: p, q, nshl, nshb, ptype
  integer :: iel, igauss, jgauss, i, j, aa, bb

  real(8) :: gp(TSP%NGAUSS), gw(TSP%NGAUSS), gwt, DetJb_SH
  real(8) :: nor(NSD), xu(NSD), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3) 
  real(8) :: thi, Dm(3,3), Dc(3,3), Db(3,3), dens, bvec(nsd), xb(NSD)

  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:), &
                          Rhs(:,:), Rhs_gra(:,:), Rhsgp(:,:), &
                          xMebe(:,:), xKebe(:,:,:), xKebegp(:,:,:), &
                          acl(:,:)

  gp = 0.0d0; gw = 0.0d0  
  DetJb_SH = 0.0d0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, TSP%NGAUSS) 

  ! loop over elements
  do iel = 1, TSP%NEL
  
    p = BEZ%P(iel); nshl = TSP%NSHL(iel); ptype = TSP%PTYPE(iel)
    q = BEZ%Q(iel); nshb = BEZ%NSHL(iel)

    allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), &
             Rhs(NSD,nshl), Rhs_gra(NSD,nshl), &
             xMebe(nshl,nshl),  &
             xKebe(NSD*NSD,nshl,nshl),  &
             xKebegp(NSD*NSD,nshl,nshl), &
             Rhsgp(NSD,nshl), lIEN(nshl), acl(nshl,nsd))

    lIEN = -1
    do i = 1, nshl
      lIEN(i) = TSP%IEN(iel,i)
    end do

    ! Get local solution arrays
    do i = 1, nshl
      acl(i,:) = ashAlpha(TSP%IEN(iel,i),:)
    end do

    ! initialization   
    xMebe   = 0.0d0     
    xKebe   = 0.0d0      ! initialize local stiffness matrix
    Rhs     = 0.0d0      ! initialize local load vector
    Rhs_gra = 0.0d0

    ! Loop over integration points (NGAUSS in each direction) 
    do jgauss = 1, TSP%NGAUSS
      do igauss = 1, TSP%NGAUSS
        
        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
        
        call eval_SHAPE_bez_sh(gp(igauss), gp(jgauss), &
                               shl, shgradl, shhessl, nor, &
                               xu, xd, dxdxi, ddxddxi,  &
                               p, q, nsd, nshl, nshb,  &
                               TSP%IEN(iel,1:nshl), TSP%NNODE, &
                               TSP%B_NET_U, TSP%B_NET_D, DetJb_SH, &
                               BEZ%Ext(iel,1:nshl,1:nshb))

        ! rotate xu back to straight-up position for evaluating the 
        ! thickness
        xb = 0.0d0
        xb(1) = cos(-BldRot)*xu(1) - sin(-BldRot)*xu(2)
        xb(2) = sin(-BldRot)*xu(1) + cos(-BldRot)*xu(2)
        xb(3) = xu(3)

        ! thickness functions....
        if (T_Flag == 1) then
          if (ptype == 1) then  ! blade skin and spar cap
            call define_thick_breitenberger(xb, thi)
          else if (ptype == 2) then  ! shear webs
            thi = 0.05d0
          else
            write(*,*) "!!!ERROR: Undefined thickness!!!"
            stop
          end if

        else
          call define_thick(xb(2), thi)
        end if

        ! extensional, coupling and bending material matrices
        Dm = 0.0d0; Dc = 0.0d0; Db = 0.0d0
        if (ptype > 0) then
          dens = 2.1d3

          ! note: ptype needs to match the material type
          Dm = matA(ptype,:,:)*thi
          Dc = matB(ptype,:,:)*thi**2
          Db = matD(ptype,:,:)*thi**3

        else
          write(*,*) "ERROR: UNDEFINED PTYPE"
          stop
        end if
        != end thickness and stiffness ===============

        gwt = gw(igauss)*gw(jgauss)
                
        ! Kirchhoff-Love Shell by J. Kiendl
        xKebegp = 0.0d0
        Rhsgp   = 0.0d0
        call e3LRhs_KLShell(shgradl, shhessl, Dm, Dc, Db, &
                            xKebegp, Rhsgp,  &
                            nshl, q, nsd, &
                            TSP%B_NET, TSP%B_NET_D,  &
                            TSP%IEN(iel,1:nshl), TSP%NNODE)
                            
        xKebe = xKebe + xKebegp*gwt
        Rhs   = Rhs   + Rhsgp*gwt
        ! end Kirchhoff-Love Shell

        ! Build Mass
        do aa = 1, nshl
          do bb = 1, nshl
            xMebe(aa,bb) = xMebe(aa,bb) + &
                           thi*dens*shl(aa)*shl(bb)*DetJb_SH*gwt
          end do
        end do

        ! Gravity effect
        bvec = g_fact*dens*9.81d0
        do aa = 1, nshl
          Rhs_gra(:,aa) = Rhs_gra(:,aa) + &
                          shl(aa)*bvec(:)*thi*DetJb_SH*gwt
        end do

      end do
    end do  ! end loop gauss points

    ! Dynamic part (RHS)
    do aa = 1, nshl
      do bb = 1, nshl
        Rhs(1,aa) = Rhs(1,aa) - xMebe(aa,bb)*acl(bb,1)
        Rhs(2,aa) = Rhs(2,aa) - xMebe(aa,bb)*acl(bb,2)
        Rhs(3,aa) = Rhs(3,aa) - xMebe(aa,bb)*acl(bb,3)
      end do
    end do

    ! Dynamic part (LHS)
    xKebe = alfi*beti*Delt*Delt*xKebe
    xKebe(1,:,:) = xKebe(1,:,:) + almi*xMebe(:,:)
    xKebe(5,:,:) = xKebe(5,:,:) + almi*xMebe(:,:)
    xKebe(9,:,:) = xKebe(9,:,:) + almi*xMebe(:,:)

    call BCLhs_3D_shell(nsd, nshl, TSP%NNODE, lIEN, TSP%IBC, xKebe)
    call BCRhs_3D_shell(nsd, nshl, TSP%NNODE, lIEN, TSP%IBC, Rhs)
    call BCRhs_3D_shell(nsd, nshl, TSP%NNODE, lIEN, TSP%IBC, Rhs_gra)

    ! Assemble load vector         
    ! Assemble thickness and lump mass
    ! LocaltoGlobal_3D is removed..
    do aa = 1, NSHL
      ! internal 
      RHSG_SH(lIEN(aa),:) = RHSG_SH(lIEN(aa),:) + Rhs(:,aa)
      ! gravity force
      RHSG_GRA_SH(lIEN(aa),:) = RHSG_GRA_SH(lIEN(aa),:) + Rhs_gra(:,aa)
    end do

    call FillSparseMat_3D_shell(nsd, nshl, lIEN, TSP%NNODE, &
                                TSP%maxNSHL, icnt, col, row, &
                                xKebe, LHSK_SH)

    deallocate(shl, shgradl, shhessl, Rhs, Rhs_gra,  &
               xMebe, xKebe, xKebegp, Rhsgp, lIEN, acl)

  end do    ! end loop elements
end subroutine IntElmAss_tsp_sh
