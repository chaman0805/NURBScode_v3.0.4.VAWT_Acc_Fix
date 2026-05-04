!=======================================================================
! L2 Projection: get the consistent mass for FEM surface
!=======================================================================
subroutine f2s_IntElmAss(FEM, icnt, col, row, nsd, LHSK, mg)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM
  integer,    intent(in)  :: icnt, nsd, col(FEM%NNODE+1), &
                             row(FEM%NNODE*50*FEM%maxNSHL)
  real(8),    intent(out) :: LHSK(NSD*NSD,icnt),  &
                             mg(FEM%NNODE)

  !  Local variables
  integer :: nshl, iel, igauss, i, j, aa, bb
  real(8) :: gp(FEM%NGAUSS,2), gw(FEM%NGAUSS), DetJb, &
             xu(NSD), xd(NSD), dxdxi(NSD,2), nor(NSD), fi(NSD)
  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:),  &        
                          lM(:), xMebe(:,:),  xKebe(:,:,:)

  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0

  ! get Gaussian points and weights
  call genGPandGW_tri(gp, gw, FEM%NGAUSS)

  ! loop over elements
  do iel = 1, FEM%NEL

  
    ! only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2), lIEN(nshl), &
               lM(nshl), xKebe(NSD*NSD,nshl,nshl), xMebe(nshl,nshl))
               
      lIEN = -1
      do i = 1, nshl
        lIEN(i) = FEM%IEN(iel,i)
      end do      

      ! initialization  
      xMebe = 0.0d0
      xKebe = 0.0d0
      lM    = 0.0d0

      ! Loop over integration points (NGAUSS in each direction) 
      do igauss = 1, FEM%NGAUSS

        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
                            xu, xd, dxdxi, nsd, nshl, &
                            FEM%IEN(iel,1:nshl), FEM%NNODE, &
                            FEM%B_NET_U, FEM%B_NET_D, DetJb)
         
        ! Build Mass and Gravity
        do aa = 1, nshl
          do bb = 1, nshl
            ! consistent mass
            xMebe(aa,bb) = xMebe(aa,bb) + shl(aa)*shl(bb)*DetJb*gw(igauss)
          end do 
          ! lumped mass
          lM(aa) = lM(aa) + shl(aa)*DetJb*gw(igauss)
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
      call FillSparseMat_3D_shell(nsd, nshl, lIEN, FEM%NNODE, &
                                  FEM%maxNSHL, icnt, col, row, &
                                  xKebe, LHSK)

      deallocate(shl, shgradl, lIEN, lM, xMebe, xKebe)
    end if
  end do    ! end loop elements
  
end subroutine f2s_IntElmAss
