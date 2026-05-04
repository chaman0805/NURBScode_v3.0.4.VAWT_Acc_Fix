!======================================================================
! Subroutine to preprocess the data. Find the list of close elements
! For NURBS or T-Spline, check the vertices and center point
! For Triangles? check the vertices should be enough...
!======================================================================
subroutine fem_find_elm(FEM, nsd, x1, ELM_DIST)
  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM
  integer,    intent(in)  :: nsd           
  real(8),    intent(in)  :: x1(nsd)
  real(8),    intent(out) :: ELM_DIST(FEM%NEL)

  !  Local variables
  integer, parameter :: nvp = 1
  integer :: iel, i, j, nshl
  real(8) :: xd(nsd), dxdxi(nsd,2), vp(nvp,2), dist(nvp)
  real(8), allocatable :: shl(:), shgradl(:,:)  

  ! vertices of the parametric domain
!!$  vp(1,1) = 1.0d0; vp(1,2) = 0.0d0
!!$  vp(2,1) = 0.0d0; vp(2,2) = 1.0d0
!!$  vp(3,1) = 0.0d0; vp(3,2) = 0.0d0
  vp(1,1) = 0.33333333333333d0
  vp(1,2) = 0.33333333333333d0

  ! initialize with a large distance
  ELM_DIST = 1.0d6

  ! loop over elements
  do iel = 1, FEM%NEL

!!$    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2))
    
      do i = 1, nvp
        ! Get Element Shape functions and their gradients
!!$        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
!!$        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
!!$        call eval_SHAPE_tri(vp(i,:), shl, shgradl, nor, &
!!$                            xu, xd, dxdxi, nsd, nshl, &
!!$                            FEM%IEN(iel,1:nshl), FEM%NNODE, &
!!$                            FEM%B_NET_U, FEM%B_NET_D, DetJb)

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(vp(i,:), nsd, nshl, &
                                 FEM%IEN(iel,1:nshl), FEM%NNODE, &
                                 FEM%B_NET_D, shl, shgradl, xd, dxdxi)

        dist(i) = sqrt(sum((xd-x1)**2))
      end do

      ! pick the shortest distance out of vertices or center
      ELM_DIST(iel) = minval(dist)

      deallocate(shl, shgradl)
!!$    end if
  end do    ! end loop elements
   
end subroutine fem_find_elm




!======================================================================
! Same as fem_find_elm, but only loop through the list of radial 
! elements
!======================================================================
subroutine fem_find_elm_rad(rad_elm_num, rad_elm_list, FEM, nsd, x1, &
                            ELM_DIST)
  use defs_shell
  implicit none

  type(mesh), intent(in)  :: FEM
  integer,    intent(in)  :: nsd, rad_elm_num, rad_elm_list(rad_elm_num)         
  real(8),    intent(in)  :: x1(nsd)
  real(8),    intent(out) :: ELM_DIST(FEM%NEL)

  !  Local variables
  integer, parameter :: nvp = 1
  integer :: iel, kel, i, j, nshl
  real(8) :: xd(nsd), dxdxi(nsd,2), vp(nvp,2), dist(nvp)

  real(8), allocatable :: shl(:), shgradl(:,:)

  ! vertices of the parametric domain
  vp = 0.0d0
  vp(1,1) = 0.333333333333333d0
  vp(1,2) = 0.333333333333333d0

  ! initialize with a large distance
  ELM_DIST = 1.0d6

  ! loop over elements
  do kel = 1, rad_elm_num

    iel = rad_elm_list(kel)

!!$    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2))
    
      do i = 1, nvp
        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(vp(i,:), nsd, nshl, &
                                 FEM%IEN(iel,1:nshl), FEM%NNODE, &
                                 FEM%B_NET_D, shl, shgradl, xd, dxdxi)

        dist(i) = sqrt(sum((xd-x1)**2))
      end do

      ! pick the shortest distance out of vertices or center
      ELM_DIST(iel) = minval(dist)

      deallocate(shl, shgradl)
!!$    end if
  end do    ! end loop elements
end subroutine fem_find_elm_rad




!======================================================================
! compute the element size (h) and record the largest one (ref config.)
!======================================================================
subroutine fem_find_elm_size(FEM, nsd)
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: FEM
  integer,    intent(in)    :: nsd

  !  Local variables
  integer :: iel, i, j, k, nshl, igauss
  real(8) :: xd(nsd), dxdxi(nsd,2), &
             gp(FEM%NGAUSS,2), gw(FEM%NGAUSS), Ginv(2,2)

  real(8), allocatable :: shl(:), shgradl(:,:)
  real(8) :: Elm_Size_Avg

  FEM%Elm_Size = 0.0d0

  gp = 0.0d0; gw = 0.0d0
  call genGPandGW_tri(gp, gw, FEM%NGAUSS)

  ! loop over elements
  do iel = 1, FEM%NEL

!!$    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2))

      do igauss = 1, FEM%NGAUSS

        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(gp(igauss,:), nsd, nshl, &
                                 FEM%IEN(iel,1:nshl), FEM%NNODE, &
                                 FEM%B_NET_U, shl, shgradl, xd, dxdxi)
        Ginv = 0.0d0
        do j = 1, 2
          do i = 1, 2
            Ginv(i,j)  = sum(dxdxi(:,i)*dxdxi(:,j))
          end do
        end do

        Elm_Size_Avg = Elm_Size_Avg + maxval(sqrt(Ginv))

        if (maxval(sqrt(Ginv)) > FEM%Elm_Size) then
          FEM%Elm_Size = maxval(sqrt(Ginv))
        end if

      end do

      deallocate(shl, shgradl)

!!$    end if
  end do    ! end loop elements
end subroutine fem_find_elm_size




!======================================================================
! Compute the physical location of the gauss points (ref config.)
! then one can build an element list based on the radial and z location
! this list will not change during the location because the radial
! location doesn't change
!======================================================================
subroutine fem_find_elm_loc(FEM, nsd)
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: FEM
  integer,    intent(in)    :: nsd

  real(8), allocatable :: shl(:), shgradl(:,:)
  integer :: iel, i, j, k, nshl, igauss
  real(8) :: xd(nsd), dxdxi(nsd,2), gp(FEM%NGAUSS,2), gw(FEM%NGAUSS), &
             Ginv(2,2)

  gp = 0.0d0; gw = 0.0d0
  call genGPandGW_tri(gp, gw, FEM%NGAUSS)

  do iel = 1, FEM%NEL

!!$    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)
      allocate(shl(nshl), shgradl(nshl,2))

      do igauss = 1, FEM%NGAUSS
        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(gp(igauss,:), nsd, nshl, &
                                 FEM%IEN(iel,1:nshl), FEM%NNODE, &
                                 FEM%B_NET_U, shl, shgradl, xd, dxdxi)

        FEM%Elm_Loc(iel,igauss,1:3) = xd

        Ginv = 0.0d0
        do j = 1, 2
          do i = 1, 2
            Ginv(i,j)  = sum(dxdxi(:,i)*dxdxi(:,j))
          end do
        end do

        FEM%Elm_Loc(iel,igauss,4) = maxval(sqrt(Ginv))

      end do

      deallocate(shl, shgradl)

!!$    end if
  end do    ! end loop elements  
end subroutine fem_find_elm_loc
!======================================================================
! Compute the normal on sliding cylinders
!======================================================================

subroutine fem_find_elm_normal(FEM, nsd, id, istep, norm)
  use defs_shell
  use mpi
  implicit none

  type(mesh), intent(inout) :: FEM
  integer,    intent(in)    :: nsd, id, istep

  real(8), intent(out) :: norm(3)

  real(8), allocatable :: shl(:), shgradl(:,:)
  integer :: iel, i, j, k, nshl, igauss
  real(8) :: xd(nsd), dxdxi(nsd,2), gp(FEM%NGAUSS,2), gw(FEM%NGAUSS), &
             Ginv(2,2), vec1(3), vec2(3)

  character(len=30) :: iname, idname

  gp = 0.0d0; gw = 0.0d0
  call genGPandGW_tri(gp, gw, FEM%NGAUSS)

  if(ismaster) then
  write(iname,*) istep
  write(idname,*) id
  open(11,file='Norm.'//trim(adjustl(iname))//'.'//trim(adjustl(idname))//'.dat',status='replace')
  end if
  norm = 0.0d0
  do iel = 1, FEM%NEL

!!$    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)
      allocate(shl(nshl), shgradl(nshl,2))

      do igauss = 1, FEM%NGAUSS
        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(gp(igauss,:), nsd, nshl, &
                                 FEM%IEN(iel,1:nshl), FEM%NNODE, &
                                 FEM%B_NET_U, shl, shgradl, xd, dxdxi)

  !Compute normal for each element:
  if ( (abs(xd(1)) <= 5.0d-2) .and. (abs(xd(2)) <= 5.0d-2)  .and. (xd(3) >= 5.0d0) .and. ismaster) then

  vec1(:) = FEM%B_NET_D(FEM%IEN(iel,2),:) - FEM%B_NET_D(FEM%IEN(iel,1),:)
  vec2(:) = FEM%B_NET_D(FEM%IEN(iel,3),:) - FEM%B_NET_D(FEM%IEN(iel,1),:)

  norm(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  norm(2) = -vec1(1)*vec2(3) + vec1(3)*vec2(1)
  norm(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

  norm(:) = norm(:)/sqrt(sum(norm(:)**2.0d0)) !Unit normal

  write(11,*) norm(:), xd(1:3)

  norm(:) = abs(norm(:)) !Output only in one direction

  end if

      end do

      deallocate(shl, shgradl)

!!$    end if
  end do    ! end loop elements 
  if(ismaster) then
  close(11)
  end if
end subroutine fem_find_elm_normal
!======================================================================
! Compute the rotation matrix to align two cylinders
!======================================================================

subroutine fem_find_rot_mat(norm, Rfix)

  use commonpars
  use mpi
  implicit none

  real(8), intent(in)  :: norm(2,3)
  real(8), intent(out) :: Rfix(3,3)
  real(8) :: theta, axis(3)

  integer :: iel, i, j, k, nshl, igauss
  real(8) :: vec1(3), vec2(3), Identity(3,3), W(3,3)

  character(len=30) :: iname, idname

  Identity = 0.0d0
  W = 0.0d0
  Rfix = 0.0d0
  do i = 1, 3
    Identity(i,i) = 1.0d0
  end do

  theta = sum(norm(1,:)*norm(2,:))/(sqrt(sum(norm(1,:)**2.0d0))*sqrt(sum(norm(2,:)**2.0d0)))
  theta = acos(theta)!*pi/180.0d0

  if(ismaster) write(*,*) 'Angle between two normals:', theta*180.0d0/pi,norm(1,:), norm(2,:)

  theta = theta!*180.0d0/pi !Minus because rotation is clockwise

  ! Define axis of rotation

  axis(1) = norm(1,2)*norm(2,3) - norm(1,3)*norm(2,2)
  axis(2) = -norm(1,1)*norm(2,3) + norm(1,3)*norm(2,1)
  axis(3) = norm(1,1)*norm(2,2) - norm(1,2)*norm(2,1)

  axis(:) = axis(:)/sqrt(sum(axis(:)**2.0d0))

  W(1,2) = -axis(3)
  W(1,3) = axis(2)
  W(2,1) = axis(3)
  W(2,3) = -axis(1)
  W(3,1) = -axis(2)
  W(3,2) = axis(1)

  ! Rodrigues formula for rotation matrix:
  if(ismaster) then
    write(*,*) 'Axis of rotation:', axis
    write(*,*) 'Rotation matrix:'
  end if
  do i = 1, 3
    do j = 1, 3
      Rfix(i,j) = Identity(i,j)*cos(theta)+(1-cos(theta))*axis(i)*axis(j)+W(i,j)*sin(theta)
    end do
  if(ismaster) write(*,*) Rfix(i,:)
  end do

end subroutine fem_find_rot_mat
