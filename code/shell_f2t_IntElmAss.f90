!=======================================================================
! L2 Projection: project quantities on FEM surface to T-Spline Surface
!=======================================================================
subroutine f2t_IntElmAss(TSP, BEZ, FEM, icnt, col, row, nsd, fg, &
                         RHSG_SH, LHSK_SH, mg_SH)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: TSP, BEZ, FEM
  integer,    intent(in)  :: icnt, nsd, col(TSP%NNODE+1), &
                             row(TSP%NNODE*50*TSP%maxNSHL)
  real(8),    intent(in)  :: fg(FEM%NNODE,NSD)
  real(8),    intent(out) :: RHSG_SH(TSP%NNODE,NSD), &
                             LHSK_SH(NSD*NSD,icnt),  &
                             mg_SH(TSP%NNODE)

  !  Local variables
  integer :: p, q, nshl, nshb, iel, igauss, jgauss, i, j, aa, bb, ct
  real(8) :: gp(TSP%NGAUSS), gw(TSP%NGAUSS), gwt, DetJb, &
             xu(nsd), xd(NSD), dxdxi(NSD,2), ddxddxi(nsd,3), nor(NSD), fi(NSD)
  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:),  shhessl(:,:), &        
                          Rhs(:,:), lM(:), xMebe(:,:),  xKebe(:,:,:)

  RHSG_SH = 0.0d0; LHSK_SH = 0.0d0; mg_SH = 0.0d0
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, TSP%NGAUSS) 

  ! loop over elements
  do iel = 1, TSP%NEL
  
    if (TSP%PTYPE(iel) == 1) then

      p = BEZ%P(iel); nshl = TSP%NSHL(iel)
      q = BEZ%Q(iel); nshb = BEZ%NSHL(iel)
      
      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), lIEN(nshl), &
               Rhs(NSD,nshl), lM(nshl), xKebe(NSD*NSD,nshl,nshl),  &
               xMebe(nshl,nshl))

      lIEN = -1
      do i = 1, nshl
        lIEN(i) = TSP%IEN(iel,i)
      end do      

      ! initialization  
      xMebe = 0.0d0
      xKebe = 0.0d0
      Rhs   = 0.0d0
      lM    = 0.0d0

      ! Loop over integration points (NGAUSS in each direction) 
      ct = 0 
      do jgauss = 1, TSP%NGAUSS
        do igauss = 1, TSP%NGAUSS
          ct = ct + 1
          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
          call eval_SHAPE_bez_sh(gp(igauss), gp(jgauss), &
                                 shl, shgradl, shhessl, nor, &
                                 xu, xd, dxdxi, ddxddxi, p, q, nsd, &
                                 nshl, nshb,  lIEN, TSP%NNODE, &
                                 TSP%B_NET_U, TSP%B_NET_D, DetJb, &
                                 BEZ%Ext(iel,1:nshl,1:nshb))

          gwt = gw(igauss)*gw(jgauss)

          ! get the closest quantities on the other surface
          fi = 0.0d0
          call f2t_int(FEM, nsd, TSP%CLE(iel,ct), TSP%CLP(iel,ct,:), &
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
      call FillSparseMat_3D_shell(nsd, nshl, lIEN, TSP%NNODE, &
                                  TSP%maxNSHL, icnt, col, row, &
                                  xKebe, LHSK_SH)

      deallocate(shl, shgradl, shhessl, lIEN, Rhs, &
                 lM, xMebe, xKebe)
    end if
  end do    ! end loop elements
  
end subroutine f2t_IntElmAss



!======================================================================
! Compute values at the closest point
!======================================================================
subroutine f2t_int(FEM, nsd, iel, clp, fg, fi)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM
  integer,    intent(in)  :: nsd, iel       
  real(8),    intent(in)  :: clp(2), fg(FEM%NNODE,nsd)
  real(8),    intent(out) :: fi(nsd)

  !  Local variables
  integer :: i, j, nshl
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), DetJb, nor(NSD)

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
  shl = 0.0d0; shgradl = 0.0d0
  xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
      call eval_SHAPE_tri(clp, shl, shgradl, nor, xu, xd, dxdxi, &
                          nsd, nshl, FEM%IEN(iel,1:nshl), FEM%NNODE, &
                          FEM%B_NET_U, FEM%B_NET_D, DetJb)

  fi = 0.0d0
  do i = 1, nshl
    fi(:) = fi(:) + fl(i,:)*shl(i)
  end do

  deallocate(shl, shgradl, fl)
   
end subroutine f2t_int
