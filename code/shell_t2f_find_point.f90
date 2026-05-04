!======================================================================
! First loop through the gauss points and compute the physical
! location on surface 1
!======================================================================
subroutine t2f_find_point(FEM, TSP, BEZ, nsd, NEL_CLOSE, CLOSE_ELM, &
                          s2s_ELM, s2s_CLP)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, FEM
  integer,    intent(in) :: nsd, NEL_CLOSE, &
                            CLOSE_ELM(FEM%NEL,FEM%NGAUSS,NEL_CLOSE)

  integer, intent(out) :: s2s_ELM(FEM%NEL,FEM%NGAUSS)
  real(8), intent(out) :: s2s_CLP(FEM%NEL,FEM%NGAUSS,3)

  !  Local variables
  integer :: iel, tel, igauss, i, j, nshl, tnewt, & 
             P_Flag, maxtel, maxtnewt
  real(8) :: gp(FEM%NGAUSS,2), gw(FEM%NGAUSS), DetJb, nor(NSD), &
             xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), cp(2), tdist

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)
                        
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  s2s_ELM = 0; s2s_CLP = 0.0d0
  maxtel = 0; maxtnewt = 0

  ! get Gaussian points and weights    
  call genGPandGW_tri(gp, gw, FEM%NGAUSS)

  ! loop over elements
  do iel = 1, FEM%NEL
!!!      do iel = 599, 599
 
    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))

      ! Loop over integration points (NGAUSS in each direction)
      do igauss = 1, FEM%NGAUSS

        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

        call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
                            xu, xd, dxdxi, nsd, nshl, &
                            FEM%IEN(iel,1:nshl), FEM%NNODE, &
                            FEM%B_NET_U, FEM%B_NET_D, DetJb)

        write(*,'(60("="))')
        write(*,*) 'FEM: iel, igp', iel, igauss
        write(*,'(60("="))')

        call tsp_find_point(TSP, BEZ, nsd, xd, NEL_CLOSE,  &
                            CLOSE_ELM(iel,igauss,:), tel, cp, &
                            tdist, tnewt, P_Flag)

        if (tel   > maxtel)   maxtel   = tel
        if (tnewt > maxtnewt) maxtnewt = tnewt

        if (P_Flag == 0) then
          write(*,'(45("%"))')
          write(*,*) '!!ERROR!! No point was found for', iel, igauss
          write(*,'(45("%"))')
!$$$              stop
        end if

        ! the closest element number and point
        s2s_ELM(iel,igauss)     = CLOSE_ELM(iel,igauss,tel)
        s2s_CLP(iel,igauss,1:2) = cp(:)
        s2s_CLP(iel,igauss,3)   = tdist
      end do

      deallocate(shl, shgradl, shhessl)
    end if
  end do    ! end loop elements
      
  write(*,*) "maxtel, maxtnewt =", maxtel, maxtnewt

end subroutine t2f_find_point
