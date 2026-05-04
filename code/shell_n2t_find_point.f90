!======================================================================
! First loop through the gauss points and compute the physical
! location on surface 1
!======================================================================
subroutine n2t_find_point(TSP, BEZ, NRB, nsd, NEL_CLOSE, CLOSE_ELM, &
                          s2s_ELM, s2s_CLP)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, NRB
  integer,    intent(in) :: nsd, NEL_CLOSE, &
                            CLOSE_ELM(TSP%NEL,TSP%NGAUSS**2,NEL_CLOSE)

  integer, intent(out) :: s2s_ELM(TSP%NEL,TSP%NGAUSS**2)
  real(8), intent(out) :: s2s_CLP(TSP%NEL,TSP%NGAUSS**2,3)

  !  Local variables
  integer :: iel, tel, igauss, jgauss, i, j, ct, p, q, nshl, tnewt, &
             P_Flag, maxtel, maxtnewt
  integer :: nshb
  real(8) :: gp(TSP%NGAUSS), gw(TSP%NGAUSS), DetJb, nor(NSD), &
             xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), cp(2), tdist

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)
                     
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  s2s_ELM = 0; s2s_CLP = 0.0d0
  maxtel = 0; maxtnewt = 0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, TSP%NGAUSS) 

  ! loop over elements
  do iel = 1, TSP%NEL
!$$$      do iel = 3583, 3583

    ! Only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( TSP%PTYPE(iel) == 1 ) then

      p    = BEZ%P(iel)
      q    = BEZ%Q(iel)
      nshl = TSP%NSHL(iel)
      nshb = BEZ%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))

      ! Loop over integration points (NGAUSS in each direction)
      ct = 0
      do jgauss = 1, TSP%NGAUSS
        do igauss = 1, TSP%NGAUSS
!$$$      do jgauss = 3, 3
!$$$        do igauss = 4, 4
 
          ct = ct + 1

          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

          call eval_SHAPE_bez_sh(gp(igauss), gp(jgauss), &
                                 shl, shgradl, shhessl, nor, &
                                 xu, xd, dxdxi, ddxddxi, p, q, nsd, &
                                 nshl, nshb,  &
                                 TSP%IEN(iel,1:nshl), TSP%NNODE, &
                                 TSP%B_NET_U, TSP%B_NET_D, DetJb, &
                                 BEZ%Ext(iel,1:nshl,1:nshb))

          write(*,'(60("="))')
          write(*,*) 'TSP: iel, igp, jgp', iel, igauss, jgauss, ct
          write(*,'(60("="))')

          call nrb_find_point(NRB, nsd, xd, NEL_CLOSE, &
                              CLOSE_ELM(iel,ct,:), tel, cp, tdist, &
                              tnewt, P_Flag)

          if (tel   > maxtel)   maxtel   = tel
          if (tnewt > maxtnewt) maxtnewt = tnewt

          if (P_Flag == 0) then
            write(*,'(45("%"))')
            write(*,*) '!!ERROR!! No point was found for', &
                       iel, igauss, jgauss, ct
            write(*,'(45("%"))')
!$$$          stop
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

end subroutine n2t_find_point
