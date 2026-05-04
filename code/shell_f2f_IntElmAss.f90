!=======================================================================
! L2 Projection: project quantities on T-Spline surface to FEM1 Surface
!=======================================================================
subroutine f2f_IntElmAss(FEM1, FEM2, icnt, col, row, nsd, fg, &
                         RHSG_SH, LHSK_SH, mg_SH)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM1, FEM2
  integer,    intent(in)  :: icnt, nsd, col(FEM1%NNODE+1), &
                             row(FEM1%NNODE*50*FEM1%maxNSHL)
  real(8),    intent(in)  :: fg(FEM2%NNODE,NSD)
  real(8),    intent(out) :: RHSG_SH(FEM1%NNODE,NSD), &
                             LHSK_SH(NSD*NSD,icnt),  &
                             mg_SH(FEM1%NNODE)

  !  Local variables
  integer :: p, q, nshl, iel, igauss, jgauss, i, j, aa, bb
  real(8) :: gp(FEM1%NGAUSS,2), gw(FEM1%NGAUSS), DetJb, &
             xu(nsd), xd(NSD), dxdxi(NSD,2), nor(NSD), fi(NSD)
  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:),  &        
                          Rhs(:,:), lM(:), xMebe(:,:),  xKebe(:,:,:)

  RHSG_SH = 0.0d0; LHSK_SH = 0.0d0; mg_SH = 0.0d0
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0

  ! get Gaussian points and weights
  call genGPandGW_tri(gp, gw, FEM1%NGAUSS)

  ! loop over elements
  do iel = 1, FEM1%NEL
  
    ! only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if (FEM1%PTYPE(iel) == 1) then

      nshl = FEM1%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2), lIEN(nshl), &
               Rhs(NSD,nshl), lM(nshl), xKebe(NSD*NSD,nshl,nshl),  &
               xMebe(nshl,nshl))

      lIEN = -1
      do i = 1, nshl
        lIEN(i) = FEM1%IEN(iel,i)
      end do      

      ! initialization  
      xMebe = 0.0d0
      xKebe = 0.0d0
      Rhs   = 0.0d0
      lM    = 0.0d0

      ! Loop over integration points (NGAUSS in each direction) 
      do igauss = 1, FEM1%NGAUSS

        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
                            xu, xd, dxdxi, nsd, nshl, &
                            FEM1%IEN(iel,1:nshl), FEM1%NNODE, &
                            FEM1%B_NET_U, FEM1%B_NET_D, DetJb)

        ! get the closest quantities on the other surface
        fi = 0.0d0
        call f2f_int(FEM2, nsd, FEM1%CLE(iel,igauss), &
                     FEM1%CLP(iel,igauss,:), fg, fi)
         
        ! Build Mass and Gravity
        do aa = 1, nshl
          do bb = 1, nshl
            ! consistent mass
            xMebe(aa,bb) = xMebe(aa,bb) + shl(aa)*shl(bb)*DetJb*gw(igauss)
          end do 
          ! lumped mass
          lM(aa) = lM(aa) + shl(aa)*DetJb*gw(igauss)
          ! right-hand-side: quantities of the other surface
          Rhs(:,aa) = Rhs(:,aa) + fi(:)*shl(aa)*DetJb*gw(igauss)
        end do

      end do ! end loop gauss points

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
      call FillSparseMat_3D_shell(nsd, nshl, lIEN, FEM1%NNODE, &
                                  FEM1%maxNSHL, icnt, col, row, &
                                  xKebe, LHSK_SH)

      deallocate(shl, shgradl, lIEN, Rhs, &
                 lM, xMebe, xKebe)
    end if
  end do    ! end loop elements
  
end subroutine f2f_IntElmAss





!======================================================================
! Compute values at the closest point
!======================================================================
subroutine f2f_int(FEM, nsd, iel, clp, fg, fi)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM
  integer,    intent(in)  :: nsd, iel       
  real(8),    intent(in)  :: clp(2), fg(FEM%NNODE,nsd)
  real(8),    intent(out) :: fi(nsd)

  !  Local variables
  integer :: i, j, nshl
  real(8) :: xd(nsd), dxdxi(nsd,2)

  real(8), allocatable :: shl(:), shgradl(:,:), fl(:,:)

  ! Check if the element is on the blade surface (e.g., ptype == 1)
  if (FEM%PTYPE(iel) /= 1 ) then
    write(*,*) "ERROR: NURBS Element", iel, " is not on the blade surface"
    stop
  end if

  nshl = FEM%NSHL(iel)
      
  allocate(shl(nshl), shgradl(nshl,2), fl(nshl,nsd))
      
  ! Get local solution arrays            
  do i = 1, nshl
    fl(i,:) = fg(FEM%IEN(iel,i),:)
  end do

  ! Get Element Shape functions and their gradients
  shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
  call eval_SHAPE_tri_fast(clp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                           FEM%NNODE, FEM%B_NET_D, &
                           shl, shgradl, xd, dxdxi)

  fi = 0.0d0
  do i = 1, nshl
    fi(:) = fi(:) + fl(i,:)*shl(i)
  end do

  deallocate(shl, shgradl, fl)
   
end subroutine f2f_int
