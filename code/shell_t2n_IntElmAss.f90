!=======================================================================
! L2 Projection: project quantities on T-Spline surface to NURBS Surface
!=======================================================================
subroutine t2n_IntElmAss(NRB, TSP, BEZ, icnt, col, row, nsd, fg, &
                         RHSG_SH, LHSK_SH, mg_SH)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: NRB, TSP, BEZ
  integer,    intent(in)  :: icnt, nsd, col(NRB%NNODE+1), &
                             row(NRB%NNODE*50*NRB%maxNSHL)
  real(8),    intent(in)  :: fg(TSP%NNODE,NSD)
  real(8),    intent(out) :: RHSG_SH(NRB%NNODE,NSD), &
                             LHSK_SH(NSD*NSD,icnt),  &
                             mg_SH(NRB%NNODE)

  !  Local variables
  integer :: p, q, nshl, iel, igauss, jgauss, i, j, aa, bb, &
             ni, nj, ct, nuk, nvk
  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), gwt, DetJb, da, &
             xu(nsd), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3), nor(NSD), fi(NSD)
  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:),  shhessl(:,:), &        
                          Rhs(:,:), lM(:), xMebe(:,:),  xKebe(:,:,:)

  RHSG_SH = 0.0d0; LHSK_SH = 0.0d0; mg_SH = 0.0d0
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 

  ! loop over elements
  do iel = 1, NRB%NEL
  
    ! get NURB coordinates
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)

    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't.
    ! Also, only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
         (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1)) .and. &
         NRB%PTYPE(iel) == 1 ) then

      ! used in calculating quadrature points. The factor of 4.0d0
      ! comes from mapping from the [-1,1] line onto a real segment...
      da = (NRB%U_KNOT(iel,ni+1)-NRB%U_KNOT(iel,ni))* &
           (NRB%V_KNOT(iel,nj+1)-NRB%V_KNOT(iel,nj))/4.0d0

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), lIEN(nshl), &
               Rhs(NSD,nshl), lM(nshl), xKebe(NSD*NSD,nshl,nshl),  &
               xMebe(nshl,nshl))

      lIEN = -1
      do i = 1, nshl
        lIEN(i) = NRB%IEN(iel,i)
      end do      

      ! initialization  
      xMebe = 0.0d0
      xKebe = 0.0d0
      Rhs   = 0.0d0
      lM    = 0.0d0

      ! Loop over integration points (NGAUSS in each direction) 
      ct = 0 
      do jgauss = 1, NRB%NGAUSS
        do igauss = 1, NRB%NGAUSS
          ct = ct + 1
          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
          call eval_SHAPE_shell(gp(igauss), gp(jgauss), shl, shgradl, shhessl, &
                                nor, xu, xd, dxdxi, ddxddxi, &
                                p, q, nsd, nshl, &
                                lIEN, NRB%NNODE, &
                                NRB%B_NET_U, NRB%B_NET_D, DetJb, &
                                ni, nj, nuk, nvk, &
                                NRB%U_KNOT(iel,1:nuk), &
                                NRB%V_KNOT(iel,1:nvk))

          gwt = gw(igauss)*gw(jgauss)*da

          ! get the closest quantities on the other surface
          fi = 0.0d0
          call t2n_int(TSP, BEZ, nsd, NRB%CLE(iel,ct), NRB%CLP(iel,ct,:), &
                       fg, fi)
         
          ! Build Mass and Gravity
          do aa = 1, nshl
            do bb = 1, nshl
              ! consistent mass
              xMebe(aa,bb) = xMebe(aa,bb) + shl(aa)*shl(bb)*DetJb*gwt
            end do 
            ! lumped mass
            lM(aa) = lM(aa) + shl(aa)*DetJb*gwt
            ! right-hand-side: quantities of the other surface
            Rhs(:,aa) = Rhs(:,aa) + fi(:)*shl(aa)*DetJb*gwt
          end do

        end do
      end do  ! end loop gauss points

      ! Assemble vectors
      do aa = 1, NSHL
        ! Right-hand-side
        RHSG_SH(lIEN(aa),:) = RHSG_SH(lIEN(aa),:) + Rhs(:,aa)
 
        ! lumped mass
        mg_SH(lIEN(aa)) = mg_SH(lIEN(aa)) + lM(aa)
      end do

      xKebe(1,:,:) = xMebe
      xKebe(5,:,:) = xMebe
      xKebe(9,:,:) = xMebe
      call FillSparseMat_3D_shell(nsd, nshl, lIEN, NRB%NNODE, &
                                  NRB%maxNSHL, icnt, col, row, &
                                  xKebe, LHSK_SH)

      deallocate(shl, shgradl, shhessl, lIEN, Rhs, &
                 lM, xMebe, xKebe)
    end if
  end do    ! end loop elements
  
end subroutine t2n_IntElmAss



!======================================================================
! Compute values at the closest point
!======================================================================
subroutine t2n_int(TSP, BEZ, nsd, iel, clp, fg, fi)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: TSP, BEZ
  integer,    intent(in)  :: nsd, iel       
  real(8),    intent(in)  :: clp(2), fg(TSP%NNODE,nsd)
  real(8),    intent(out) :: fi(nsd)

  !  Local variables
  integer :: i, j, p, q, nshl, nshb
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), DetJb, nor(NSD)

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:), fl(:,:)
    
  ! Check if the element is on the blade surface (e.g., ptype == 1)
  if (TSP%PTYPE(iel) /= 1 ) then
    write(*,*) "ERROR: T-Spline Element", iel, " is not on the blade surface"
    stop
  end if

  p = BEZ%P(iel); nshl = TSP%NSHL(iel)
  q = BEZ%Q(iel); nshb = BEZ%NSHL(iel)
      
  allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), &
           fl(nshl,nsd))
      
  ! Get local solution arrays            
  do i = 1, nshl
    fl(i,:) = fg(TSP%IEN(iel,i),:)
  end do

  ! Get Element Shape functions and their gradients
  shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
  xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
  call eval_SHAPE_bez_sh(clp(1), clp(2), &
                         shl, shgradl, shhessl, nor, &
                         xu, xd, dxdxi, ddxddxi, p, q, nsd, &
                         nshl, nshb,  &
                         TSP%IEN(iel,1:nshl), TSP%NNODE, &
                         TSP%B_NET_U, TSP%B_NET_D, DetJb, &
                         BEZ%Ext(iel,1:nshl,1:nshb))
  fi = 0.0d0
  do i = 1, nshl
    fi(:) = fi(:) + fl(i,:)*shl(i)
  end do

  deallocate(shl, shgradl, shhessl, fl)
   
end subroutine t2n_int
