!======================================================================
! First loop through the gauss points and compute the physical
! location on surface 1
!======================================================================
subroutine f2f_find_point(FEM1, FEM2, nsd, NEL_CLOSE, CLOSE_ELM, &
                          s2s_ELM, s2s_CLP, maxtel)
  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(in) :: FEM1, FEM2
  integer,    intent(in) :: nsd, NEL_CLOSE, &
                            CLOSE_ELM(FEM1%NEL,FEM1%NGAUSS,NEL_CLOSE)

  integer, intent(out) :: s2s_ELM(FEM1%NEL,FEM1%NGAUSS), maxtel
  real(8), intent(out) :: s2s_CLP(FEM1%NEL,FEM1%NGAUSS,3)

  !  Local variables
  integer :: iel, tel, igauss, i, j, nshl, tnewt, & 
             maxtnewt
  real(8) :: gp(FEM1%NGAUSS,2), gw(FEM1%NGAUSS), xd(nsd), dxdxi(nsd,2), &
             cp(2), tdist

  real(8), allocatable :: shl(:), shgradl(:,:)
                        
  gp = 0.0d0; gw = 0.0d0
  s2s_ELM = 0; s2s_CLP = 0.0d0
  maxtel = 0; maxtnewt = 0

  ! get Gaussian points and weights    
  call genGPandGW_tri(gp, gw, FEM1%NGAUSS)

  ! loop over elements
  do iel = 1, FEM1%NEL
 
!!$    if (FEM1%PTYPE(iel) == 1) then

      nshl = FEM1%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2))

      ! Loop over integration points (NGAUSS in each direction)
      do igauss = 1, FEM1%NGAUSS

        ! Get Element Shape functions and their gradients
!!$        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
!!$        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
!!$        call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
!!$                            xu, xd, dxdxi, nsd, nshl, &
!!$                            FEM1%IEN(iel,1:nshl), FEM1%NNODE, &
!!$                            FEM1%B_NET_U, FEM1%B_NET_D, DetJb)

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(gp(igauss,:), nsd, nshl, &
                                 FEM1%IEN(iel,1:nshl), FEM1%NNODE, &
                                 FEM1%B_NET_D, shl, shgradl, xd, dxdxi)

!!$        write(*,'(60("="))')
!!$        write(*,*) 'FEM: iel, igp', iel, igauss
!!$        write(*,'(60("="))')

        call fem_find_point(FEM2, nsd, xd, NEL_CLOSE,  &
                            CLOSE_ELM(iel,igauss,:), tel, cp, &
                            tdist, tnewt)

        if (tel   > maxtel)   maxtel   = tel
        if (tnewt > maxtnewt) maxtnewt = tnewt


        ! the closest element number and point
        s2s_ELM(iel,igauss)     = CLOSE_ELM(iel,igauss,tel)
        s2s_CLP(iel,igauss,1:2) = cp(:)
        s2s_CLP(iel,igauss,3)   = tdist
      end do

      deallocate(shl, shgradl)
!!$    end if
  end do    ! end loop elements
      
  if (ismaster) then
    write(*,*) "maxtel, NEL_CLOSE =", maxtel, NEL_CLOSE
  end if
end subroutine f2f_find_point
