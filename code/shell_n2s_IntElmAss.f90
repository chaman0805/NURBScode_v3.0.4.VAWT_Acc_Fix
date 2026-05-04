!=======================================================================
! L2 Projection: get the consistent mass for NURBS surface
!=======================================================================
subroutine n2s_IntElmAss(NRB, icnt, col, row, nsd, LHSK, mg)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: NRB
  integer,    intent(in)  :: icnt, nsd, col(NRB%NNODE+1), &
                             row(NRB%NNODE*50*NRB%maxNSHL)
  real(8),    intent(out) :: LHSK(NSD*NSD,icnt),  &
                             mg(NRB%NNODE)

  !  Local variables
  integer :: p, q, nshl, iel, igauss, jgauss, i, j, aa, bb, ct
  integer :: ni, nj, nuk, nvk
  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), gwt, DetJb, da, &
             xu(NSD), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3), nor(NSD), fi(NSD)
  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:),  shhessl(:,:), &        
                          lM(:), xMebe(:,:),  xKebe(:,:,:)

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
         mod(NRB%PTYPE(iel),2) == 1 ) then

      ! used in calculating quadrature points. The factor of 4.0d0
      ! comes from mapping from the [-1,1] line onto a real segment...
      da = (NRB%U_KNOT(iel,ni+1)-NRB%U_KNOT(iel,ni))* &
           (NRB%V_KNOT(iel,nj+1)-NRB%V_KNOT(iel,nj))/4.0d0

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), lIEN(nshl), &
               lM(nshl), xKebe(NSD*NSD,nshl,nshl), xMebe(nshl,nshl))
               

      lIEN = -1
      do i = 1, nshl
        lIEN(i) = NRB%IEN(iel,i)
      end do      

      ! initialization  
      xMebe = 0.0d0
      xKebe = 0.0d0
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
         
          ! Build Mass and Gravity
          do aa = 1, nshl
            do bb = 1, nshl
              ! consistent mass
              xMebe(aa,bb) = xMebe(aa,bb) + shl(aa)*shl(bb)*DetJb*gwt
            end do 
            ! lumped mass
            lM(aa) = lM(aa) + shl(aa)*DetJb*gwt
          end do

        end do
      end do  ! end loop gauss points

      ! Assemble vectors
      do aa = 1, NSHL 
        ! lumped mass
        mg(lIEN(aa)) = mg(lIEN(aa)) + lM(aa)
      end do

      xKebe(1,:,:) = xMebe
      xKebe(5,:,:) = xMebe
      xKebe(9,:,:) = xMebe
      call FillSparseMat_3D_shell(nsd, nshl, lIEN, NRB%NNODE, &
                                  NRB%maxNSHL, icnt, col, row, &
                                  xKebe, LHSK)

      deallocate(shl, shgradl, shhessl, lIEN, lM, xMebe, xKebe)
    end if
  end do    ! end loop elements
  
end subroutine n2s_IntElmAss
