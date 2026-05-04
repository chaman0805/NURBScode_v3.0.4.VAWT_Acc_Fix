!======================================================================
!
!======================================================================
subroutine IntElmAss_mass   

  use aAdjKeep
  use commonvars
  use mpi

  implicit none
  integer :: elemCount, nodeCount, itmp
  
  ! Local variables
  integer :: iel, igauss, idx, i, j, k, hess_flag, aa, bb, NGAUSS
  
  integer :: nshl

  real(8), allocatable :: shlu(:), shgradgu(:,:), shhessgu(:,:,:)
  real(8), allocatable :: xMebe(:,:), xl(:,:), dl(:,:), wl(:)
  real(8), allocatable :: gp(:,:), gw(:)
  real(8) :: rtmp, vol_ex, vol_int, volm, vole, nodec(NNODE), dxidx(NSD,NSD)
  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)
  
  LHSmass = 0d0
 


  ! loop over elements to compute isotropic forcing
  elemCount = 0
  vol_ex    = 0
  vol_int   = 0
  
  do iel = 1, NELEM
  
    NSHL   = ELMNSHL(iel)
    NGAUSS = ELMNGAUSS(iel)
    allocate(shlu(NSHL), shgradgu(NSHL,NSD), shhessgu(NSHL,NSD,NSD),&     
             xMebe(NSHL,NSHL), xl(NSHL,NSD),dl(NSHL,NSD),wl(NSHL),&
             gp(NGAUSS,NSD), gw(NGAUSS))
	
    ! get Gaussian points and weights  
    call genGPandGW(gp, gw, NGAUSS)

    do i = 1, NSHL
      idx = IEN(iel,i)
      xl(i,:) = xg(idx,:)
      dl(i,:) = dg(idx,:)
      wl(i)   = wg(idx)
    end do

    xMebe = 0d0
    vole  = 0d0
    elemCount = elemCount + 1

    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS

      ! Get Element Shape functions and their gradients
      shlu     = 0d0   ! initialize
      shgradgu = 0d0
      hess_flag = 0
      call eval_shape(nshl, iel, gp(igauss,:), xl, dl, wl, shlu, shgradgu, &
                 	  shhessgu, dxidx, Gij, Ginv, hess_flag)

      vole = vole + DetJ*gw(igauss)

      ! Calculate element lfull mass
      do aa = 1, NSHL
        do bb = 1, NSHL
          xMebe(aa,bb) = xMebe(aa,bb) + &
     		         shlu(aa)*shlu(bb)*DetJ*gw(igauss)
        end do
      end do	   

    end do

    call FillSparseMat_fm(nshl, iel, xMebe)
    
    vol_int = vol_int + vole

    if (.not.iga) then
      volm = 0.0d0
      call voltet(xl(3,:), xl(2,:), xl(1,:), xl(4,:), volm)
      vol_ex = vol_ex + volm
    end if
	
    deallocate(shlu, shgradgu, shhessgu)
    deallocate(xMebe, xl, dl, wl, gp, gw)
  end do
  
  nodec = 1
  if (numnodes > 1) then
    rtmp  = vol_ex
    call MPI_ALLREDUCE(rtmp, vol_ex, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
     
    rtmp  = vol_int
    call MPI_ALLREDUCE(rtmp, vol_int, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
     
    itmp = elemCount
    call MPI_ALLREDUCE(itmp, elemCount, 1, MPI_INTEGER, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if
       
  nodec = 1.0d0  
  if (numnodes > 1) then
    call commu(nodec, 1,'in ')
  end if

  nodeCount = 0
  do i = 1, NNODE
    if (nodec(i) >= 0.5d0) nodeCount = nodeCount + 1
  end do  
  
  if (numnodes > 1) then
    itmp = nodeCount
    call MPI_ALLREDUCE(itmp, nodeCount, 1, MPI_INTEGER, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  
  end if
    
  hglob = (vol_int/nodeCount)**(1.0d0/3.0d0)

  if (ismaster) then
    if (.not.iga) write(*,"(a,ES12.4)") " Exact  Volume   = ", vol_ex
                  write(*,"(a,ES12.4)") " Integ. Volume   = ", vol_int
                  write(*,"(a, I12)")   " Element count   = ", elemCount
                  write(*,"(a, I12)")   " Node    count   = ", nodeCount
                  write(*,"(a,ES12.4)") " Avg. Elem. size = ", hglob     
  end if

end subroutine IntElmAss_mass


!======================================================================
! Compute Tangent contribution
!======================================================================
subroutine FaceAssembly_area()
  
  use aAdjKeep
  use mpi
  use commonvars
  implicit none
  
  integer :: b, ifac, idx, iel, igauss, i
  
  integer :: nshl

  real(8) :: gp(NGAUSSb,2), gw(NGAUSSb), mgp(NGAUSSb,3)
  
  real(8), allocatable :: xl(:,:), dl(:,:), wl(:)     
  real(8), allocatable :: shlu(:), shgradgu(:,:)
  real(8) :: nor(NSD), dxidx(NSD,NSD)
  real(8) :: Gij(NSD,NSD), Ginv(NSD,NSD)
  real(8) :: area, larea

  ! shb will be the shape function array while shbg will hold the
  ! gradients of the shape functions
  call genGPandGWb(gp, gw, NGAUSSb)
  area = 0.0d0

  ! Loop over Faces 
  do b = 1, NBOUND
    do ifac = 1, bound(b)%NFACE
  
      call genGPMap(NGAUSSb, bound(b)%FACE_OR(ifac), iga, mgp)

      iel = bound(b)%F2E(ifac)
	 
      NSHL = ELMNSHL(iel)
     
      allocate(xl(NSHL,NSD), dl(NSHL,NSD), wl(NSHL), shlu(NSHL),  &
	       shgradgu(NSHL,NSD))
                           
      do i = 1, NSHL
        idx = IEN(iel,i)
        xl(i,:) = xg(idx,:)
        dl(i,:) = dg(idx,:)
        wl(i)   = wg(idx)
      end do
    
      larea = 0.0d0
      ! Loop over integration points    
      do igauss = 1, NGAUSSb
             
        call eval_faceshape(nshl, iel, gp(igauss,:), mgp(igauss,:), &
                            bound(b)%FACE_OR(ifac), &
                            xl, dl, wl, shlu, shgradgu, dxidx, &
                            Gij, Ginv, nor)
     
        larea = larea + DetJb*gw(igauss)
      end do
      area = area + larea
	  
	  deallocate(xl, dl, wl, shlu, shgradgu)
    end do   
  end do  
     
  if (numnodes > 1) then
    larea = area
    call MPI_ALLREDUCE(larea,area, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if
  
  if (ismaster) then
    write(*,"(a,ES12.4)") " Surface area    = ", area  
  end if
    
end subroutine FaceAssembly_area
