!======================================================================
! For every gauss point on surface 1, find the closest point on 
! surface 2. The point need to be between -1 and 1 in the parametric
! domain.
!======================================================================
subroutine fem_find_point(FEM, nsd, x1, NEL_CLOSE, CLOSE_ELM, &
                          tkel, tp, tdist, tnewt)

  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(in) :: FEM
  integer,    intent(in) :: nsd, NEL_CLOSE, CLOSE_ELM(NEL_CLOSE)
  real(8),    intent(in) :: x1(nsd)
  integer,    intent(out) :: tkel, tnewt
  real(8),    intent(out) :: tp(2), tdist

  integer, parameter :: nvp = 3
  real(8), parameter :: eps = 1.0d-13

  !  Local variables
  integer :: iel, kel, tel, i, j, nshl, inewt, nnewt
  real(8) :: gp(2), DetJb, &
             xu(nsd), xd(nsd), dxdxi(nsd,2), nor(NSD), &
             dist, xdist(nsd), vp(nvp,2), &
             tx1(nsd), txi(nsd), &
             xt1(nsd), xt2(nsd), xt3(nsd), &
             matA(2,2), vecB(2), rtmp

  real(8), allocatable :: shl(:), shgradl(:,:)
  real(8) :: lgp(2)

  ! vertices of the parametric domain
  vp(1,1) = 1.0d0; vp(1,2) = 0.0d0
  vp(2,1) = 0.0d0; vp(2,2) = 1.0d0
  vp(3,1) = 0.0d0; vp(3,2) = 0.0d0

  tp    = 0.0d0    ! target point
  tdist = 1.0d6    ! initialize the distance
  tel   = 0
  tkel  = 0
  tnewt = 0
  txi   = 0.0d0
  tx1   = x1

  ! loop over elements
  do kel = 1, NEL_CLOSE
     
    ! get the real element number
    iel = CLOSE_ELM(kel)
    if (iel == -1) exit

    ! Only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( FEM%PTYPE(iel) == 1 ) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2))
     
      !--------------------------------------------------------
      ! Find the perpendicular (closest) points.
      ! Linear problem for triangles
      ! xt are the physical coordinates of the 3 vertices
      !--------------------------------------------------------
      gp   = 0.0d0    ! initial guess
      dist = 0.0d0
 
      ! Get the physical coordinates of the vertices
      ! input vp, output xt
      shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
      call eval_SHAPE_tri_fast(vp(1,:), nsd, nshl, FEM%IEN(iel,1:nshl), &
                               FEM%NNODE, FEM%B_NET_D, &
                               shl, shgradl, xt1, dxdxi)

      shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
      call eval_SHAPE_tri_fast(vp(2,:), nsd, nshl, FEM%IEN(iel,1:nshl), &
                               FEM%NNODE, FEM%B_NET_D, &
                               shl, shgradl, xt2, dxdxi)

      shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
      call eval_SHAPE_tri_fast(vp(3,:), nsd, nshl, FEM%IEN(iel,1:nshl), &
                               FEM%NNODE, FEM%B_NET_D, &
                               shl, shgradl, xt3, dxdxi)

      matA = 0.0d0; vecB = 0.0d0
      matA(1,1) = sum((xt1-xt3)*(xt1-xt3))
      matA(1,2) = sum((xt1-xt3)*(xt2-xt3))
      matA(2,1) = sum((xt1-xt3)*(xt2-xt3))
      matA(2,2) = sum((xt2-xt3)*(xt2-xt3))

      vecB(1) = sum((x1-xt3)*(xt1-xt3))
      vecB(2) = sum((x1-xt3)*(xt2-xt3))

      rtmp = matA(1,1)*matA(2,2)-matA(1,2)*matA(2,1)

      gp(1) = ( matA(2,2)*vecB(1)-matA(1,2)*vecB(2))/rtmp
      gp(2) = (-matA(2,1)*vecB(1)+matA(1,1)*vecB(2))/rtmp

      !---------------------------------------------------------------
      ! This is what we originally want. Find the normal point
      ! that satisfies 0<=gp(1)<=1, 0<=gp(2)<=1, 0<=1-gp(1)-gp(2)<=1 
      !---------------------------------------------------------------
      if ( (gp(1)>=0.0d0 .and. gp(1)<=1.0d0) .and. &
           (gp(2)>=0.0d0 .and. gp(2)<=1.0d0) .and. &
           ( (1.0d0-gp(1)-gp(2))>=0.0d0 .and. &
             (1.0d0-gp(1)-gp(2))<=1.0d0       )      ) then

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(gp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! there could be multiple solution
        ! only recorded the smallest distance
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = gp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if


      !-----------------------------------------------------
      ! if gp is located outside of the criterion, then 
      ! special treatment is needed
      !-----------------------------------------------------
      ! first deal with vertices
      !-----------------------------------------------------
      else if (gp(1)<0.0d0 .and. gp(2)<0.0d0) then

        lgp(1) = 0.0d0
        lgp(2) = 0.0d0

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(lgp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! store the shortest distance at the vertices
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = lgp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if

      else if (gp(1)<0.0d0 .and. gp(2)>1.0d0) then

        lgp(1) = 0.0d0
        lgp(2) = 1.0d0

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(lgp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! store the shortest distance at the vertices
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = lgp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if

      else if (gp(1)>1.0d0 .and. gp(2)<0.0d0) then

        lgp(1) = 1.0d0
        lgp(2) = 0.0d0

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(lgp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! store the shortest distance at the vertices
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = lgp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if

 
      !-----------------------------------------------------
      ! edges
      !-----------------------------------------------------
      else if (gp(1)<0.0d0 .and. gp(2)>=0.0d0 .and. gp(2)<=1.0d0) then

        lgp(1) = 0.0d0

        ! need to solve for lgp(2)
        lgp(2) = sum((x1-xt3)*(xt2-xt3))/sum((xt2-xt3)*(xt2-xt3))
        if (lgp(2) < 0.0d0) lgp(2) = 0.0d0
        if (lgp(2) > 1.0d0) lgp(2) = 1.0d0        

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(lgp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! store the shortest distance to the edges
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = lgp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if

      else if (gp(1)>=0.0d0 .and. gp(1)<=1.0d0 .and. gp(2)<0.0d0) then

        lgp(1) = sum((x1-xt3)*(xt1-xt3))/sum((xt1-xt3)*(xt1-xt3))
        if (lgp(1) < 0.0d0) lgp(1) = 0.0d0
        if (lgp(1) > 1.0d0) lgp(1) = 1.0d0  
        lgp(2) = 0.0d0

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(lgp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! store the shortest distance to the edges
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = lgp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if

      !-----------------------------------------------------
      ! the rest... project to the x1-x2 edge
      !-----------------------------------------------------
      else

        matA = 0.0d0; vecB = 0.0d0
        matA(1,1) = 1.0d0
        matA(1,2) = 1.0d0
        matA(2,1) = sum((xt1-xt2)*(xt1-xt3))
        matA(2,2) = sum((xt1-xt2)*(xt2-xt3))

        vecB(1) = 1.0d0
        vecB(2) = sum((x1-xt3)*(xt1-xt2))

        rtmp = matA(1,1)*matA(2,2)-matA(1,2)*matA(2,1)

        lgp(1) = ( matA(2,2)*vecB(1)-matA(1,2)*vecB(2))/rtmp
        lgp(2) = (-matA(2,1)*vecB(1)+matA(1,1)*vecB(2))/rtmp
        
        if (lgp(1) > 1.0d0 .or. lgp(2) < 0.0d0) then
          lgp(1) = 1.0d0; lgp(2) = 0.0d0
        end if

        if (lgp(1) < 0.0d0 .or. lgp(2) > 1.0d0) then
          lgp(1) = 0.0d0; lgp(2) = 1.0d0
        end if

        !--------------

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(lgp, nsd, nshl, FEM%IEN(iel,1:nshl), &
                                 FEM%NNODE, FEM%B_NET_D, &
                                 shl, shgradl, xd, dxdxi)

        dist = sqrt(sum((xd-x1)**2))

        ! store the shortest distance to the edges
        if (dist < tdist) then 
          ! update the target results
          tdist    = dist
          tp       = lgp
          tel      = iel
          tkel     = kel
          txi      = xd
        end if

      end if

      deallocate(shl, shgradl)

    end if  ! end if ptype

  end do  ! end loop elements

!  if (ismaster) then
!    write(*,*) 'iel, kel =', tel, tkel
!    write(*,*) 'gp  =', tp(1), tp(2)
!    write(*,*) 'loc1=', tx1
!    write(*,*) 'loc2=', txi
!    write(*,*) 'dist=', tdist
!  end if

  if (tp(1) > 1.0d0 .or. tp(1) < 0.0d0 .or. &
      tp(2) > 1.0d0 .or. tp(2) < 0.0d0) stop

  if (tdist >= 5.0d-1 .and. ismaster) then
    write(*,'(45("%"))')
    write(*,*) "!!WARNNING!! Large Distance:", tdist
    write(*,'(45("%"))')
!!!    stop
  end if

end subroutine fem_find_point
