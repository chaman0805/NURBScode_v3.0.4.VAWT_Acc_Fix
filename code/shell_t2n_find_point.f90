!======================================================================
! First loop through the gauss points and compute the physical
! location on surface 1
!======================================================================
subroutine t2n_find_point(NRB, TSP, BEZ, nsd, NEL_CLOSE, CLOSE_ELM, &
                          s2s_ELM, s2s_CLP)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, NRB
  integer,    intent(in) :: nsd, NEL_CLOSE, &
                            CLOSE_ELM(NRB%NEL,NRB%NGAUSS**2,NEL_CLOSE)

  integer, intent(out) :: s2s_ELM(NRB%NEL,NRB%NGAUSS**2)
  real(8), intent(out) :: s2s_CLP(NRB%NEL,NRB%NGAUSS**2,3)

  !  Local variables
  integer :: iel, tel, igauss, jgauss, i, j, ct, p, q, nshl, tnewt, & 
             P_Flag, maxtel, maxtnewt
  integer :: nuk, nvk, ni, nj
  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), DetJb, nor(NSD), &
             xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), cp(2), tdist

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)
                        
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  s2s_ELM = 0; s2s_CLP = 0.0d0
  maxtel = 0; maxtnewt = 0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 

  ! loop over elements
  do iel = 1, NRB%NEL
!!!      do iel = 599, 599
 
    ! get NURB coordinates
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)

    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't.
    ! Also, only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
         (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1)) .and. &
         NRB%PTYPE(iel) == 1 ) then

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))

      ! Loop over integration points (NGAUSS in each direction)
      ct = 0
      do jgauss = 1, NRB%NGAUSS
        do igauss = 1, NRB%NGAUSS
!      do jgauss = 2, 2
!       do igauss = 2, 2
 
          ct = ct + 1

          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

          call eval_SHAPE_shell(gp(igauss), gp(jgauss), &
                                shl, shgradl, shhessl,  &
                                nor, xu, xd, dxdxi, ddxddxi,  &
                                p, q, nsd, nshl,  &
                                NRB%IEN(iel,1:nshl), NRB%NNODE, &
                                NRB%B_NET_U, NRB%B_NET_D, DetJb, &
                                ni, nj, nuk, nvk, &
                                NRB%U_KNOT(iel,1:nuk), &
                                NRB%V_KNOT(iel,1:nvk))

          write(*,'(60("="))')
          write(*,*) 'NRB: iel, igp, jgp', iel, igauss, jgauss, ct
          write(*,'(60("="))')

          call tsp_find_point(TSP, BEZ, nsd, xd, NEL_CLOSE,  &
                              CLOSE_ELM(iel,ct,:), tel, cp, tdist, &
                              tnewt, P_Flag)

          if (tel   > maxtel)   maxtel   = tel
          if (tnewt > maxtnewt) maxtnewt = tnewt

          if (P_Flag == 0) then
            write(*,'(45("%"))')
            write(*,*) '!!ERROR!! No point was found for', &
                       iel, igauss, jgauss, ct
            write(*,'(45("%"))')
!$$$              stop
          end if

          ! the closest element number and point
          s2s_ELM(iel,ct)     = CLOSE_ELM(iel,ct,tel)
          s2s_CLP(iel,ct,1:2) = cp(:)
          s2s_CLP(iel,ct,3)   = tdist

        end do
      end do

      deallocate(shl, shgradl, shhessl)
    end if
  end do    ! end loop elements
      
  write(*,*) "maxtel, maxtnewt =", maxtel, maxtnewt

end subroutine t2n_find_point
