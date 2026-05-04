!======================================================================
! FaceAssembly to build local RHS and LHS for non-matching boundaries.
!======================================================================
subroutine FaceAsse_DG_w1t1(dgAlpha, ugAlpha, ugmAlpha, pgAlpha)
  use aAdjKeep
  use mpi
  use commonvars
  implicit none

  real(8), intent(in) :: dgAlpha(NNODE,NSD), ugAlpha(NNODE,NSD), &
                         ugmAlpha(NNODE,NSD), pgAlpha(NNODE)                        
  ! Local variables
  real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:)       
  real(8), allocatable :: dl(:,:), ul(:,:), uml(:,:), pl(:), xl(:,:), &
                          wl(:)
  real(8), allocatable :: xKebe11(:,:,:), xMebe(:,:), xGebe(:,:,:), &
                          xDebe1(:,:,:)
  real(8), allocatable :: Rhsu(:,:), Rhsp(:)

  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, FID
  
  real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)    
  real(8) :: ui(NSD), pri, duidxi(NSD,NSD), Gij(NSD,NSD), &
             Ginv(NSD,NSD), ti(NSD), nor(NSD), umi(NSD), &
             dxidx(NSD,NSD)
  real(8) :: tmp1(NSD), tmp2(NSD,NSD), tauB

  NSHL = -1

  call genGPandGWb(gp, gw, NGAUSSb)

  ! loop over faces
  do b = 1, NBOUND

    ! non-matching boundaries
    FID = bound(b)%Face_ID
    if (FID==11 .or. FID==12) then

      do ifac = 1, bound(b)%NFACE
        
        iel = bound(b)%F2E(ifac)
               
        if (NSHL /= ELMNSHL(iel)) then
          if (NSHL >= 0) then
            deallocate(shlu, shgradgu, shhessgu, dl, ul, uml, pl, xl, &
                       wl, xKebe11, xMebe, xGebe, xDebe1, Rhsu, Rhsp)
          end if

          NSHL = ELMNSHL(iel)

          allocate(shlu(NSHL), shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD), &
                   dl(NSHL,NSD), ul(NSHL,NSD), uml(NSHL,NSD), pl(NSHL), &
                   xl(NSHL,NSD), wl(NSHL), &
                   xKebe11(NSD*NSD,NSHL,NSHL), xMebe(NSHL,NSHL), &
                   xGebe(NSD,NSHL,NSHL), xDebe1(NSD,NSHL,NSHL), &      
                   Rhsu(NSD,NSHL), Rhsp(NSHL))
        end if

        ! get the mapping between face element and volume surface
        call genGPMap(nshl, NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)
           
        ! Get local solution arrays   	       	     
        do i = 1, NSHL
          j = IEN(iel,i)
          xl (i,:) = xg      (j,:)
          dl (i,:) = dgAlpha (j,:)
          ul (i,:) = ugAlpha (j,:)
          uml(i,:) = ugmAlpha(j,:)
          pl (i)   = pgAlpha (j)
          wl (i)   = wg      (j)
        end do
     
        ! Get the density and viscosity. 
        ! Now just set it to water and constant for each element
        rho = rhow
        mu  = muw

        ! initialize local resistance matrix	
        xKebe11 = 0.0d0         
        xGebe   = 0.0d0
        xDebe1  = 0.0d0
        xMebe   = 0.0d0
    
        ! initialize local load vector
        Rhsu = 0.0d0    
        Rhsp = 0.0d0   

        ! Loop over integration points (NGAUSSb in each direction) 
        do igauss = 1, NGAUSSb

          call eval_faceshape(nshl, iel, gp(igauss,:), mgp(igauss,:), &
                              bound(b)%FACE_OR(ifac), &
                              xl, dl, wl, shlu, shgradgu, dxidx, &
                              Gij, Ginv, nor)

          ! Interpolate	    
          pri = sum(pl*shlu)
     
          do i = 1, NSD
            ui (i) = sum(ul (:,i)*shlu)
            umi(i) = sum(uml(:,i)*shlu)
          end do
    
          do j = 1, NSD
            do i = 1, NSD
              duidxi(i,j) = sum(ul(:,i)*shgradgu(:,j))
            end do
          end do

          ! traction
          ti(:) = -pri*nor(:) + mu*(duidxi(:,1)*nor(1) + &
                                    duidxi(:,2)*nor(2) + &
                                    duidxi(:,3)*nor(3) + &
                                    duidxi(1,:)*nor(1) + &
                                    duidxi(2,:)*nor(2) + &
                                    duidxi(3,:)*nor(3))

!!!          ti(:) = -pri*nor(:)

          call e3bSTAB_DG(tauB, ui-umi, nor, dxidx)

          call e3bRHS_DG(nshl, shlu, shgradgu, gw(igauss), ui, &
                         ti, nor, tauB, Rhsu, Rhsp, 1.0d0)

          call e3bLHS_DG(nshl, shlu, shgradgu, ui, umi, duidxi, tauB, &
                         gw(igauss), nor, xKebe11, xGebe, xDebe1)     
        end do

        call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsp)
        call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe) 
        call LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp)

      end do
    end if   
  end do
   
  if (NSHL >= 0) then
    deallocate(shlu, shgradgu, shhessgu, dl, ul, xl, wl, &
               xKebe11, xMebe, xGebe, xDebe1, Rhsu, Rhsp)
  end if
end subroutine FaceAsse_DG_w1t1






!======================================================================
! FaceAssembly to build non-local RHS for Matrix-Free method at the
! non-matching boundaries. For the approximation of R'du
!======================================================================
subroutine FaceAsse_DG_w1t2(dgAlpha, ugAlpha, ugmAlpha, pgAlpha, &
                            RHSGu2, RHSGp2, NM)
  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell
  implicit none
     
  type(shell_nmb), intent(inout) :: NM
       
  real(8), intent(in) :: dgAlpha(NNODE,NSD),  ugAlpha(NNODE,NSD), &
                         ugmAlpha(NNODE,NSD), pgAlpha(NNODE)

  real(8), intent(out) :: RHSGu2(NNODE,NSD), RHSGp2(NNODE)

  ! Local variables
  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, &
             FID
  
  real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)

  real(8), allocatable :: shlu(:), shgradgu(:,:)
  real(8), allocatable :: dl(:,:), ul(:,:), uml(:,:), &
                          pl(:), xl(:,:), wl(:), tl(:,:)
  real(8), allocatable :: Rhsu(:,:), Rhsp(:) 

  real(8) :: ui(NSD), pri, duidxi(NSD,NSD), dxidx(NSD,NSD), &
             Gij(NSD,NSD), Ginv(NSD,NSD), ti(NSD), &
             nor(NSD), umi(NSD)

  real(8) :: tauB, unor, upos, uneg, sgn

  ! lumped mass and traction
  real(8) :: NMLMs(NNODE), NMTra(NNODE,NSD), NMRHS1(NNODE,NSD), &
             NMRHS2(NNODE,NSD), tgLocal(NNODE,NSD), &
             ugLocal(NNODE,NSD), pgLocal(NNODE,NSD), r2tmp(NNODE,NSD)

  ! First compute the traction (or viscous force) then put
  ! on the surface nodes (stored in the volume array)
  NMLMs = 0.0d0; NMTra = 0.0d0
  call FaceAsse_Tra(dgAlpha, ugAlpha, pgAlpha, NMLMs, NMTra)

  ! project the traction, velocity and pressure
  ! NMTra already has zero at slaves
  tgLocal = 0.0d0
  call FaceAsse_Proj(NM%FEM(1), NM%FEM(2), NMTra, tgLocal)
  call FaceAsse_Proj(NM%FEM(2), NM%FEM(1), NMTra, tgLocal)

  ! we want to zero out ug at slaves for the projection
  r2tmp = ugAlpha
  call zeroslaves(r2tmp, NSD)
  ugLocal = 0.0d0
  call FaceAsse_Proj(NM%FEM(1), NM%FEM(2), r2tmp, ugLocal)
  call FaceAsse_Proj(NM%FEM(2), NM%FEM(1), r2tmp, ugLocal)

  ! instead of projecting the traction, project the pressure only,
  ! then use the local normal vector
!!!  r2tmp      = 0.0d0
!!!  r2tmp(:,1) = pgAlpha
!!!  call zeroslaves(r2tmp, NSD)
!!!  pgLocal = 0.0d0
!!!  call FaceAsse_Proj(NM%FEM(1), NM%FEM(2), r2tmp, pgLocal)
!!!  call FaceAsse_Proj(NM%FEM(2), NM%FEM(1), r2tmp, pgLocal)


  NSHL = -1

  call genGPandGWb(gp, gw, NGAUSSb)

  ! loop over faces
  do b = 1, NBOUND

    FID = bound(b)%Face_ID
    if (FID==11 .or. FID==12) then

      do ifac = 1, bound(b)%NFACE
  
        iel = bound(b)%F2E(ifac)
               
        if (NSHL /= ELMNSHL(iel)) then
          if (NSHL >= 0) then
            deallocate(shlu, shgradgu, dl, ul, uml, pl,  &
                       xl, wl, Rhsu, Rhsp)
          end if

          NSHL = ELMNSHL(iel)

          allocate(shlu(NSHL), shgradgu(NSHL,NSD), dl(NSHL,NSD),     &
                   ul(NSHL,NSD), uml(NSHL,NSD), tl(NSHL,NSD),     &
                   pl(NSHL), xl(NSHL,NSD), wl(NSHL), & 
                   Rhsu(NSD,NSHL), Rhsp(NSHL))
        end if

        ! get the mapping between face element and volume surface
        call genGPMap(nshl, NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)

        ! Get local solution arrays   	       	     
        do aa = 1, NSHL
          k = IEN(iel,aa)
          xl (aa,:) = xg      (k,:)
          dl (aa,:) = dgAlpha (k,:)
          wl (aa)   = wg      (k)
          uml(aa,:) = ugmAlpha(k,:)

          ul (aa,:) = ugLocal (k,:)
          tl (aa,:) = tgLocal (k,:)
!!!          pl (aa)   = pgLocal (k,1)
        end do

        ! Get the density and viscosity. 
        ! Now just set it to water and constant for each element
        rho = rhow
        mu  = muw

        ! initialize local load vector
        Rhsu = 0.0d0    
        Rhsp = 0.0d0   

        ! Loop over integration points (NGAUSSb in each direction) 
        do igauss = 1, NGAUSSb

          call eval_faceshape(nshl, iel, gp(igauss,:), mgp(igauss,:), &
                              bound(b)%FACE_OR(ifac), &
                              xl, dl, wl, shlu, shgradgu, dxidx, &
                              Gij, Ginv, nor)
          
          ! Get the quantities from projection
!!!          pri = sum(pl*shlu)
     
          do i = 1, NSD
            ui(i)  = sum(ul(:,i)*shlu)
            umi(i) = sum(uml(:,i)*shlu)
            ti(i)  = sum(tl(:,i)*shlu)
          end do     

          ! traction is computed using the projected pressure and
          ! the local normal vector. (notice the sign of "pri")
!!!          ti(:) = pri*nor(:)

          call e3bSTAB_DG(tauB, ui-umi, nor, dxidx)

          call e3bRHS_DG(nshl, shlu, shgradgu, gw(igauss), ui, &
                         ti, nor, tauB, Rhsu, Rhsp, -1.0d0)
   
        end do

        ! Assemble RHS
        do aa = 1, NSHL
          k = IEN(iel,aa)      
          RHSGu2(k,:) = RHSGu2(k,:) + Rhsu(:,aa)    
          RHSGp2(k)   = RHSGp2(k)   + Rhsp(aa)    
        end do

      end do

    end if

  end do

  if (NSHL >= 0) then
    deallocate(shlu, shgradgu, dl, ul, uml, pl, xl, wl, Rhsu, Rhsp)
  end if
end subroutine FaceAsse_DG_w1t2




!======================================================================
! Projection between non-matching surfaces
!======================================================================
subroutine FaceAsse_Proj(FEM1, FEM2, NMRHS1, NMRHS2)  
  use aAdjKeep
  use commonvars
  use mpi
  use defs_shell
 
  implicit none 

  !type(shell_nmb), intent(inout) :: NM
  type(mesh), intent(inout) :: FEM1, FEM2
  real(8),    intent(in)    :: NMRHS1(NNODE,NSD)
  real(8),    intent(inout) :: NMRHS2(NNODE,NSD)

  real(8), allocatable :: r2tmp1(:,:), r2tmp2(:,:)

  integer :: i, bb

  !--------------------------------------------------------------
  ! assign the local quantities from volume element boundary to
  ! shell mesh
  !--------------------------------------------------------------
  ! loop through each non-matching boundaries
  allocate(r2tmp1(FEM1%NNODE,NSD))
  r2tmp1 = 0.0d0

  bb = FEM1%iBound
  do i = 1, bound(bb)%NNODE
    r2tmp1(bound(bb)%L2SNODE(i),:) = NMRHS1(bound(bb)%BNODES(i),:)
  end do

  FEM1%FORCE = r2tmp1

  ! sum up the local force such that the shell force vector is complete
  if (numnodes > 1) then
    call MPI_ALLREDUCE(r2tmp1, FEM1%FORCE, FEM1%NNODE*NSD, &
                       MPI_DOUBLE_PRECISION,                   &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if

  deallocate(r2tmp1)

  ! project from surface 1 to surface 2
  call f2f_l2project_LM(FEM1, FEM2, nsd, FEM1%FORCE, &
                        FEM2%FORCE)

!  write(*,*) FEM2%FORCE
!  stop
  !--------------------------------------------------------------
  ! map the solutions back to volume boundaries
  !--------------------------------------------------------------
  bb = FEM2%iBound
  do i = 1, bound(bb)%NNODE
    NMRHS2(bound(bb)%BNODES(i),:) = FEM2%FORCE(bound(bb)%L2SNODE(i),:)
  end do

  if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

end subroutine FaceAsse_Proj




!=======================================================================
! L2Projection: project quantities from surface 1 to surface 2
! Element loop through surface 2, get quantities from closest points
! on surface 1
!=======================================================================
subroutine f2f_l2project_LM(FEM1, FEM2, nsd, fg1, fg2)
  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM2, FEM1
  integer,    intent(in)  :: nsd
  real(8),    intent(in)  :: fg1(FEM1%NNODE,NSD)
  real(8),    intent(out) :: fg2(FEM2%NNODE,NSD)
                             
  !  Local variables
  integer :: nshl, iel, igauss, i, j, aa, bb
  real(8) :: gp(FEM2%NGAUSS,2), gw(FEM2%NGAUSS), DetJb, &
             xu(nsd), xd(NSD), dxdxi(NSD,2), nor(NSD), fi(NSD), &
             mg2(FEM2%NNODE)
  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:), Rhs(:,:), lm(:)

  fg2 = 0.0d0; mg2 = 0.0d0
  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0

  ! get Gaussian points and weights
  call genGPandGW_tri(gp, gw, FEM2%NGAUSS)

  ! loop over elements
  do iel = 1, FEM2%NEL

    nshl = FEM2%NSHL(iel)

    allocate(shl(nshl), shgradl(nshl,2), lIEN(nshl), &
             Rhs(NSD,nshl), lM(nshl))

    lIEN = -1
    do i = 1, nshl
      lIEN(i) = FEM2%IEN(iel,i)
    end do      

    ! initialization  
    Rhs   = 0.0d0
    lM    = 0.0d0

    ! Loop over integration points (NGAUSS in each direction) 
    do igauss = 1, FEM2%NGAUSS

      ! Get Element Shape functions and their gradients
      shl = 0.0d0; shgradl = 0.0d0
      xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
      call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
                          xu, xd, dxdxi, nsd, nshl, &
                          FEM2%IEN(iel,1:nshl), FEM2%NNODE, &
                          FEM2%B_NET_U, FEM2%B_NET_D, DetJb)

      ! get the closest quantities on the other surface
      fi = 0.0d0
      call f2f_vint(FEM1, nsd, FEM2%CLE(iel,igauss), &
                    FEM2%CLP(iel,igauss,:), fg1, fi)
         
      do aa = 1, nshl
        ! lumped mass
        lM(aa) = lM(aa) + shl(aa)*DetJb*gw(igauss)
        ! right-hand-side: quantities of the other surface
        Rhs(:,aa) = Rhs(:,aa) + fi(:)*shl(aa)*DetJb*gw(igauss)
      end do

    end do ! end loop gauss points

    ! Assemble vectors
    do aa = 1, NSHL
      ! Right-hand-side
      fg2(lIEN(aa),:) = fg2(lIEN(aa),:) + Rhs(:,aa)
 
      ! lumped mass
      mg2(lIEN(aa)) = mg2(lIEN(aa)) + lM(aa)
    end do


    deallocate(shl, shgradl, lIEN, Rhs, lM)
  end do    ! end loop elements
  
  ! devided by lumped mass to get nodal traction
  forall (i = 1:FEM2%NNODE, mg2(i) > 1.0d-15) 
    fg2(i,:) = fg2(i,:)/mg2(i)
  end forall

end subroutine f2f_l2project_LM



!======================================================================
! Compute values at the closest point on surface 1
!======================================================================
subroutine f2f_vint(FEM, nsd, iel, clp, fg, fi)
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

  nshl = FEM%NSHL(iel)
      
  allocate(shl(nshl), shgradl(nshl,2), fl(nshl,nsd))
      
  ! Get local solution arrays            
  do i = 1, nshl
    fl(i,:) = fg(FEM%IEN(iel,i),:)
  end do

  ! Get Element Shape functions and their gradients
  shl = 0.0d0; shgradl = 0.0d0; xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
  call eval_SHAPE_tri(clp, shl, shgradl, nor, xu, xd, dxdxi, &
                      nsd, nshl, FEM%IEN(iel,1:nshl), FEM%NNODE, &
                      FEM%B_NET_U, FEM%B_NET_D, DetJb)

  fi = 0.0d0
  do i = 1, nshl
    fi(:) = fi(:) + fl(i,:)*shl(i)
  end do

  deallocate(shl, shgradl, fl)
end subroutine f2f_vint



!======================================================================
! Build traction on the surface, then do l2-projection put to the 
! control points
! NMLM:  Lumped mass for non-matching
! NMTra: Traction for non-matching
!======================================================================
subroutine FaceAsse_Tra(dgAlpha, ugAlpha, pgAlpha, NMLMs, NMTra)
  use aAdjKeep
  use mpi
  use commonvars
  implicit none
          
  real(8), intent(in) :: dgAlpha(NNODE,NSD),   ugAlpha(NNODE,NSD), &
                         pgAlpha(NNODE)

  real(8), intent(out) :: NMLMs(NNODE), NMTra(NNODE,NSD)

  ! Local variables
  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, &
             FID
  
  real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)
    
  real(8), allocatable :: shlu(:), shgradgu(:,:)
  real(8), allocatable :: dl(:,:), ul(:,:), pl(:), &
                          xl(:,:), wl(:)
  real(8), allocatable :: rhs1(:), rhs2(:,:)

  real(8) :: pri, duidxi(NSD,NSD), dxidx(NSD,NSD), &
             Gij(NSD,NSD), Ginv(NSD,NSD), ti(NSD), &
             nor(NSD)
		
  NSHL = -1
	 
  call genGPandGWb(gp, gw, NGAUSSb)

  ! loop over faces
  do b = 1, NBOUND

    ! Non-matching boundaries
    FID = bound(b)%Face_ID
    if (FID==11 .or. FID==12) then

      do ifac = 1, bound(b)%NFACE
  
        iel = bound(b)%F2E(ifac)

        if (NSHL /= ELMNSHL(iel)) then
          if (NSHL >= 0) then
            deallocate(shlu, shgradgu, dl, ul, pl, xl, wl, rhs1, rhs2)
          end if

          NSHL = ELMNSHL(iel)
          allocate(shlu(NSHL), shgradgu(NSHL,NSD), dl(NSHL,NSD),  &
                   ul(NSHL,NSD), pl(NSHL), xl(NSHL,NSD), wl(NSHL), &
                   rhs1(NSHL), rhs2(NSD,NSHL))
        end if
	
        ! get the mapping between face element and volume surface
        call genGPMap(nshl, NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)
     	  
        ! Get local solution arrays   	       	     
        do aa = 1, NSHL
          k = IEN(iel,aa)
          xl  (aa,:) = xg       (k,:)
          dl  (aa,:) = dgAlpha  (k,:)
          ul  (aa,:) = ugAlpha  (k,:)
          pl  (aa)   = pgAlpha  (k)
          wl  (aa)   = wg       (k)
        end do

        ! Get the density and viscosity. 
        ! Now just set it to water and constant for each element
        rho = rhow
        mu  = muw

        ! initialize local load vector
        rhs1 = 0.0d0; rhs2 = 0.0d0
	  
        ! Loop over integration points
        do igauss = 1, NGAUSSb

          call eval_faceshape(nshl, iel, gp(igauss,:), mgp(igauss,:), &
                              bound(b)%FACE_OR(ifac), &
                              xl, dl, wl, shlu, shgradgu, dxidx, &
                              Gij, Ginv, nor)
          
          !-----------------------------------------------
          ! Get the quantities at the intergration point
          !-----------------------------------------------

          ! pressure
          pri = sum(pl*shlu)	      
	    
          ! velocity gradient
          do j = 1, NSD
            do i = 1, NSD
              duidxi(i,j) = sum(ul(:,i)*shgradgu(:,j))
            end do
          end do

          ! traction
          ti(:) = -pri*nor(:) !+ mu*(duidxi(:,1)*nor(1) + &
                              !      duidxi(:,2)*nor(2) + &
                              !      duidxi(:,3)*nor(3) + &
                              !      duidxi(1,:)*nor(1) + &
                              !      duidxi(2,:)*nor(2) + &
                              !      duidxi(3,:)*nor(3))

          do aa = 1, NSHL
            ! lumped mass
            rhs1(aa) = rhs1(aa) + shlu(aa)*DetJb*gw(igauss)
            ! traction
            rhs2(:,aa) = rhs2(:,aa) + shlu(aa)*ti(:)*DetJb*gw(igauss)
          end do

        end do

         ! Assemble LM and Tra
        do aa = 1, NSHL
          k = IEN(iel,aa)    
          NMLMs(k)   = NMLMs(k)   + rhs1(aa)  
          NMTra(k,:) = NMTra(k,:) + rhs2(:,aa)
        end do

      end do ! NFACE
    end if
  end do  ! NBOUND

  if (NSHL >= 0) then
    deallocate(shlu, shgradgu, dl, ul, pl, xl, wl, rhs1, rhs2)
  end if

  ! communicate
  if (numnodes > 1) then
    call commu(NMLMs, 1  , 'in ')
    call commu(NMTra, NSD, 'in ')
  end if

  ! devided by lumped mass to get nodal traction
  forall (i = 1:NNODE, NMLMs(i) > 1.0d-15) 
    NMTra(i,:) = NMTra(i,:)/NMLMs(i)
  end forall

end subroutine FaceAsse_Tra
