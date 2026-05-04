!======================================================================
!
!======================================================================
subroutine IntElmAss_3D(inewt, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                        acgmAlpha, pgAlpha, phigAlpha, FFlag)
  use aAdjKeep
  use commonvars 
  use mpi
  implicit none
  
  integer, intent(in) :: inewt, FFlag
  real(8), intent(in) :: dgAlpha(NNODE,NSD),   ugAlpha(NNODE,NSD), &
                         acgAlpha(NNODE,NSD),  ugmAlpha(NNODE,NSD), &
                         acgmAlpha(NNODE,NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE)
  ! Local variables  
  real(8), parameter :: damp = 0.5d0
  
  integer :: iel, igauss, i, j, hess_flag, idx, nshl, NGAUSS
  
  real(8) :: volm, vol_ex, vol_int, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, pri, kappa, k_dc, tauP
 
  real(8), allocatable :: xKebe11(:,:,:), xGebe(:,:,:), xDebe1(:,:,:), &
                          xMebe(:,:), Rhsu(:,:), Rhsm(:,:), Rhsp(:)

  real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:)

  real(8), allocatable :: fl(:,:), dl(:,:), ul(:,:),wl(:), acl(:,:), &
                          uml(:,:), acml(:,:), pl(:), dlold(:,:),    &
                          xl(:,:), dumb(:,:), phil(:), ulold(:,:)

  real(8), allocatable :: gp(:,:), gw(:)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD),    &
             fi(NSD), uadvi(NSD), xi(NSD), ddidxi(NSD,NSD), phi, &
             dphidxi(NSD), duidxi(NSD,NSD),                      &
             duidxixj(NSD,NSD,NSD),duiolddxi(NSD,NSD),           &
             dxidx(NSD,NSD), dpridxi(NSD), rLi(NSD)

  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(NSD,NSD), Stens(NSD,NSD), Ctens(NSD,NSD,NSD,NSD)  

  volm    = 0.0d0
  vol_ex  = 0.0d0
  vol_int = 0.0d0
  cfl     = 0.0d0

  NGAUSS = -1
  NSHL   = -1

  ! loop over elements
  do iel = 1, NELEM
       
    if (NSHL /= ELMNSHL(iel)) then
      
      if (NSHL >= 0) then
        deallocate(shlu, shgradgu, shhessgu, gp, gw, &
                   fl, dl, ul, wl, acl, uml, acml, pl, dlold, &
                   xl, dumb, phil, ulold, & 
                   xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsm, Rhsp)
      end  if
     
      NSHL   = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)

      allocate(xKebe11(NSD*NSD,NSHL,NSHL), xGebe(NSD,NSHL,NSHL),&
               xDebe1(NSD,NSHL,NSHL), xMebe(NSHL,NSHL),&
               Rhsu(NSD,NSHL), Rhsm(NSD,NSHL), Rhsp(NSHL),&
               shlu(NSHL), shgradgu(NSHL,NSD),shhessgu(NSHL,NSD,NSD), &
               fl(NSHL,NSD), dl(NSHL,NSD), ul(NSHL,NSD),wl(NSHL), &
               acl(NSHL,NSD), uml(NSHL,NSD), acml(NSHL,NSD),  &
               pl(NSHL), dlold(NSHL,NSD), xl(NSHL,NSD),  &
               dumb(NSHL,NSD), phil(NSHL),ulold(NSHL,NSD),&
               gp(NGAUSS,NSD), gw(NGAUSS))

      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if

    fl = 0.0d0
!    do i = 1, NSHL
!      fl(i,:) = LD_el(iel,:)
!    end do
    
    ! Get local solution arrays 
    do i = 1, NSHL
      idx = IEN(iel,i)
      xl   (i,:) = xg       (idx,:)
      dl   (i,:) = dgAlpha  (idx,:)
      dlold(i,:) = dgold    (idx,:)
      wl   (i)   = wg       (idx)
      ul   (i,:) = ugAlpha  (idx,:)
      ulold(i,:) = ugold    (idx,:)
      acl  (i,:) = acgAlpha (idx,:)
      uml  (i,:) = ugmAlpha (idx,:)
      acml (i,:) = acgmAlpha(idx,:)
      pl   (i)   = pgAlpha  (idx)
      phil (i)   = phigAlpha(idx)
    end do
         
    dumb = 0.0d0
    
    ! initialize local stiffness matrix
    xKebe11 = 0.0d0     
    xGebe   = 0.0d0
    xDebe1  = 0.0d0
    xMebe   = 0.0d0
    
    ! initialize local load vector 
    Rhsu = 0.0d0    
    Rhsm = 0.0d0 
    Rhsp = 0.0d0
     
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS
     
      ! Get Element Shape functions and their gradients
      ! initialize
      shlu     = 0.0d0
      shgradgu = 0.0d0
      shhessgu = 0.0d0
      hess_flag = NS_hess_flag
      
      ! Fluid
      if (EL_TYP(iel) == 0) then
     
        call eval_shape(nshl, iel, gp(igauss,:), xl, dl, wl, shlu, &
                        shgradgu, shhessgu, dxidx, Gij, Ginv, hess_flag)
        
        vol_int = vol_int + DetJ*gw(igauss)

        call e3int_fluid(nshl, xl,dl, ul, acl, uml, acml, &
                         pl, fl, phil, shlu, shgradgu, &
                         shgradgu, shhessgu,dxidx, Ginv, &
                         di, ui, aci, umi, acmi, pri, fi, &
                         ddidxi, duidxi, duidxixj, dpridxi, & 
                         phi, dphidxi, xi, rLi)

        do i = 1, NSD
          do j = 1, NSD
            duiolddxi(i,j) = sum(ulold(:,i)*shgradgu(:,j))
          end do
        end do


        ! ALE Advective Velocity
        uadvi(:) = ui(:) - umi(:)

        tauM = 0.0d0; tauP = 0.0d0; tauC = 0.0d0; tauBar = 0.0d0

        call e3STAB_3D(Gij, Ginv, uadvi, rLi, &
                       tauM, tauP, tauC, tauBar, uprime, cfl_loc)

        cfl(1)  = max(cfl(1), cfl_loc(1))
!!!        cfl(2)  = max(cfl(2), cfl_loc(2))
   
        k_dc = 0.0d0
!!!        call  e3DC_beta2(uadvi,duiolddxi,Gij,rLi,dxidx, k_dc)

        if (FFlag == 0) then
          call e3LHS_3D_fluid_Old(nshl, ui, umi, aci, pri, duidxi, dpridxi, &
                                  rLi, tauM, tauP, tauC, tauBar, &
                                  k_dc, gw(igauss), shlu, shgradgu, &
                                  shhessgu, xKebe11, xGebe, &
                                  xDebe1, xMebe)

!!!        call e3LHS_3D_fluid(ui, umi, aci, pri, duidxi, &
!!!                            dpridxi, rLi, tauM,  tauC, k_dc, &
!!!                            gw(igauss), shlu, shgradgu,  &
!!!                            shgradgu, shhessgu, &
!!!                            xKebe11, &
!!!                            xGebe, xDebe1, xMebe)
        end if

        call e3RHS_3D_fluid(nshl, ui, aci, umi, acmi, uadvi,&
                            pri, rLi, fi, duidxi, ddidxi, &
                            tauM, tauP, tauC, tauBar, k_dc,&
                            gw(igauss), shlu, shgradgu, uprime, &
                            Rhsu, Rhsp)

!!!    if (isnan( sum(xKebe11*xKebe11))) then 	  	     
!!!      write(*,*)  sum(xKebe11*xKebe11)   
!!!    end if		 

      ! Solid
      else if (EL_TYP(iel) == 1) then
       
        Ftens = 0.0d0; Stens = 0.0d0; Ctens = 0.0d0
       
        call eval_shape(nshl, iel, gp(igauss,:), xl, dumb, wl, shlu, &
                        shgradgu, shhessgu, dxidx, Gij, Ginv, hess_flag)
       
        vol_int = vol_int + DetJ*gw(igauss)

        call e3int_struct(nshl, dl, ul, acl, fl, shlu, shgradgu, &
                          di, ui, aci, fi, ddidxi)
       
       
        call e3tensors_nonlin(NSD, ddidxi, mu, kappa, Ftens, &
                              Stens, Ctens)

        Stens = Stens + FPKS(iel,:,:)
       
        if (FFlag == 0) then
          call e3LHS_3D_struct(nshl, gw(igauss), shlu, shgradgu, Ftens, &
                               Stens, Ctens, xKebe11, damp)
        end if

        call e3Rhs_3D_struct(nshl, aci, ui, fi, Ftens, Stens, gw(igauss), &
                             shlu, shgradgu, Rhsu, damp)
       
      end if
    
    end do    
     
!    volm = 0.0d0
!    call voltet(xl(3,:), xl(2,:), xl(1,:), xl(4,:), volm)
!    vol_ex = vol_ex + volm

    if (FFlag == 0) then      
      call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
                    xMebe, Rhsu, Rhsp)

      ! Assemble global LHS Matrix
      call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe)
    end if

    ! Assemble load RHS vector
    call LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp) 

  end do

  deallocate(xKebe11, xGebe, shlu, shgradgu, shhessgu, &
             fl, dl, ul, wl, acl, uml, acml, pl, dlold, &
             xl, dumb, phil, ulold, &
             Rhsu, Rhsm, Rhsp, xDebe1, xMebe, gp, gw)

  if (FFlag == 0) then
    ! Find largest CFL-number and output to screen
    if (numnodes > 1) then
      gcfl = cfl
      call MPI_ALLReduce(gcfl, cfl, 2, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, MPI_COMM_WORLD,mpi_err) 
    end if
    if (ismaster) then
      write(*,'(40("-"))') 
      write(*,*) "    CFL = ", cfl(1)
    end if
  end if

end subroutine IntElmAss_3D



!======================================================================
!
!======================================================================     
subroutine IntElmMesh_3D(dgAlpha)
  
  use aAdjKeep
  use commonvars

  implicit none
  
  real(8), intent(in) :: dgAlpha(NNODE,NSD)
  
  integer :: nshl
  
  ! Local variables  
  integer :: iel, igauss, i, hess_flag, idx, NGAUSS

  real(8) :: volm
  real(8), allocatable :: gp(:,:), gw(:)
  real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:), &
                          fl(:,:), dl(:,:), ul(:,:),wl(:), acl(:,:), &
                          uml(:,:), acml(:,:), pl(:), dlold(:,:), &
                          xl(:,:), phil(:)

  real(8), allocatable :: xKebe22(:,:,:), Rhsm(:,:)
      
  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), pri, &
             fi(NSD), ddidxi(NSD,NSD), duidxi(NSD,NSD),  &
             duidxixj(NSD,NSD,NSD), &
             dxidx(NSD,NSD), dpridxi(NSD), rLi(NSD)
  
  real(8) :: kappa, phi, dphidxi(NSD), xi(NSD)

  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)

  NGAUSS = -1
  NSHL   = -1 

  ! loop over elements
  do iel = 1, NELEM
     
    if (NSHL /= ELMNSHL(iel)) then
        
      if (NSHL >= 0) then
        deallocate(shlu, shgradgu, shhessgu, fl, dl, ul, wl, acl, &
                   uml, acml, pl, dlold, xl, phil, xKebe22, Rhsm, &
                   gp, gw)
      end if
     
      NSHL   = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
        
      allocate(shlu(NSHL), shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD), &
               fl(NSHL,NSD), dl(NSHL,NSD), ul(NSHL,NSD), wl(NSHL), &
               acl(NSHL,NSD), uml(NSHL,NSD), acml(NSHL,NSD),  &
               pl(NSHL), dlold(NSHL,NSD), xl(NSHL,NSD), phil(NSHL), &
               xKebe22(NSD*NSD,NSHL,NSHL), Rhsm(NSD,NSHL), &
               gp(NGAUSS,NSD), gw(NGAUSS))
      
      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if
  
    rho    = 1.0d0
    mu     = 1.0d0  !3448275.86206897d0
    lambda = 1.0d0  !31034482.7586207d0
    kappa  = 0.0d0

    ! Get local solution arrays  
    fl = 0.0d0
    do i = 1, NSHL
      idx = IEN(iel,i)
      xl   (i,:) = xg     (idx,:)
      dl   (i,:) = dgAlpha(idx,:)
      dlold(i,:) = dgold  (idx,:)
      wl   (i)   = wg     (idx)
      ul   (i,:) = 0.0d0
      acl  (i,:) = 0.0d0
      uml  (i,:) = 0.0d0
      acml (i,:) = 0.0d0
      pl   (i)   = 0.0d0
      phil (i)   = 0.0d0
    end do
    
    ! initialize local stiffness matrix
    xKebe22 = 0.0d0
    
    ! initialize local load vector
    Rhsm = 0.0d0
    
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS
     
      ! Get Element Shape functions and their gradients
      ! initialize
      shlu      = 0.0d0
      shgradgu  = 0.0d0
      shhessgu  = 0.0d0
      hess_flag = 0
         
      ! Fluid
      if (EL_TYP(iel) == 0) then 
                
        call eval_shape(nshl, iel, gp(igauss,:), xl, dlold, wl, shlu, &
                        shgradgu, shhessgu, dxidx, Gij, Ginv, hess_flag)
               
        call e3int_fluid(nshl, xl, dl-dlold, ul, acl, uml, acml, &
                         pl, fl, phil, shlu, shgradgu, &
                         shgradgu, shhessgu, dxidx, Ginv, &
                         di, ui, aci, umi, acmi, pri, fi, &
                         ddidxi, duidxi, duidxixj, dpridxi, phi, &
                         dphidxi,xi, rLi)
    
        ! vval determines Young's modulus     
        call e3LHS_3D_mesh(nshl, gw(igauss), shgradgu, xKebe22)
        
        call e3Rhs_3D_mesh(nshl, ddidxi, gw(igauss), shlu, shgradgu, Rhsm)

      end if     
    end do

    ! Modify local matrix for Dirichlet boundary conditions
    call BCMesh_3D(nshl, iel, xKebe22, Rhsm)     
    call FillSparseMesh_3D(nshl, iel, xKebe22)
    call LocaltoGlobal_3D_mesh(nshl, iel, Rhsm) ! Assemble load vector
     
  end do


  deallocate(shlu, shgradgu, shhessgu, fl, dl, ul, wl, acl, uml, acml, & 
             pl, dlold, xl, phil, xKebe22, Rhsm, gp, gw)

end subroutine IntElmMesh_3D
