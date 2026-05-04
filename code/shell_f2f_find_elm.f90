!======================================================================
! We first check and save the n% closest elements for the gauss points
! on f1. This list is then used for finding closest points. 
! Check the distance between f1 gauss point and f2 element vertices.
!======================================================================
subroutine f2f_find_elm(FEM1, FEM2, nsd, NEL_CLOSE, CLOSE_ELM)
  use mpi
  use defs_shell
  implicit none

  type(mesh), intent(in) :: FEM1, FEM2
  integer,    intent(in) :: nsd, NEL_CLOSE

  ! Only sort for the elements that have distance smaller than 
  ! "elm_fact" times the element size (h)
  real(8), parameter :: elm_fact = 1.2d0

  ! local variables
  integer :: iel, i, j, k, nshl, igauss, rad_elm_num
  real(8) :: xd(nsd), dxdxi(nsd,2), Ginv(2,2), &
             gp(FEM1%NGAUSS,2), gw(FEM1%NGAUSS)

!!$  real(8) :: xu(nsd), DetJb, nor(NSD)

  integer :: ELM_NUM(FEM2%NEL), CLOSE_ELM(FEM1%NEL,FEM1%NGAUSS,NEL_CLOSE)
  real(8) :: ELM_DIST(FEM2%NEL), ELM_DIST_TMP(FEM2%NEL), ELM_DIST_CRI

  real(8), allocatable :: shl(:), shgradl(:,:)

  ELM_NUM = -1

  ELM_DIST     = 1.0d6
  ELM_DIST_TMP = 1.0d6

  gp = 0.0d0; gw = 0.0d0
  call genGPandGW_tri(gp, gw, FEM1%NGAUSS)

  ! loop over elements
  do iel = 1, FEM1%NEL

!!$    if (FEM1%PTYPE(iel) == 1) then

      nshl = FEM1%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2))

      do igauss = 1, FEM1%NGAUSS

        ! Get Element Shape functions and their gradients
!!$        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0; Ginv = 0.0d0
!!$        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
!!$        call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
!!$                            xu, xd, dxdxi, nsd, nshl, &
!!$                            FEM1%IEN(iel,1:nshl), FEM1%NNODE, &
!!$                            FEM1%B_NET_U, FEM1%B_NET_D, DetJb)

        shl = 0.0d0; shgradl = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0
        call eval_SHAPE_tri_fast(gp(igauss,:), nsd, nshl, &
                                 FEM1%IEN(iel,1:nshl), FEM1%NNODE, &
                                 FEM1%B_NET_D, shl, shgradl, xd, dxdxi)

        ! find the distance to vertices or center on the other surface
!!!        call fem_find_elm(FEM2, nsd, xd, ELM_DIST)

        rad_elm_num = FEM1%RAD_ELM_NUM(iel,igauss)
        call fem_find_elm_rad(rad_elm_num,                                 &
                              FEM1%RAD_ELM_LIST(iel,igauss,1:rad_elm_num), &
                              FEM2, nsd, xd, ELM_DIST)

        ! Get the element size and set ELM_DIST_CRI for sorting 
        ! because you don't want to sort all the elements. 
        ! e.g., 10.0*maxval(sqrt(Ginv)) means only sort for the elements
        ! that have distance smaller than 10 times of the element size
!        Ginv = 0.0d0
!        do j = 1, 2
!          do i = 1, 2
!            Ginv(i,j)  = sum(dxdxi(:,i)*dxdxi(:,j))
!          end do
!        end do
!        ELM_DIST_CRI = elm_fact*maxval(sqrt(Ginv))

        ELM_DIST_CRI = elm_fact*FEM1%Elm_Size

        ! Only sort elements with distance smaller than criterion
        k = 0
        do i = 1, FEM2%NEL
          if (ELM_DIST(i) < ELM_DIST_CRI) then
            k = k + 1
            ELM_DIST_TMP(k) = ELM_DIST(i)
            ELM_NUM(k)  = i
          end if
        end do

!!$        ! creat the element list (unsorted)
!!$        ELM_NUM = (/ (i, i=1,FEM2%NEL) /)

        ! sort the distance
!!$        call in_sort(ELM_DIST, ELM_NUM, FEM2%NEL)
        call in_sort(ELM_DIST_TMP(1:k), ELM_NUM(1:k), k)

        ! store the NEL_CLOSEth shortest distance elements
        CLOSE_ELM(iel,igauss,1:NEL_CLOSE) = ELM_NUM(1:NEL_CLOSE)
!!!        CLOSE_ELM(iel,igauss,1:k) = ELM_NUM(1:k)
      end do
      deallocate(shl, shgradl)

!!$    end if
!!$    if (mod(iel,1000)==0 .and. ismaster) write(*,*) "  finish FEM1 element", iel, FEM1%NEL
  end do    ! end loop elements
  
end subroutine f2f_find_elm




!======================================================================
! Build a list of elements based on z and radius.
! Then one can build the element list based on this list.
! The key is this list does not change with rotation, so only need 
! to compute once.
!======================================================================
subroutine f2f_find_elm_rad(FEM1, FEM2, nsd, Center_Rot)
  use mpi
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: FEM1
  type(mesh), intent(in)    :: FEM2
  integer,    intent(in)    :: nsd
  real(8), intent(in)       :: Center_rot(3)

  real(8), parameter :: ring_fact = 2.0d0
  integer :: iel, jel, igauss, jgauss
  integer :: num_rad_elm, Elm_Flag
  real(8) :: z1, r1, z2, r2, he

  ! loop over elements on surface 1
  do iel = 1, FEM1%NEL

!!$    if (FEM1%PTYPE(iel) == 1) then

      ! for every gauss point on surface 1
      do igauss = 1, FEM1%NGAUSS      

        num_rad_elm = 0
        ! z and radius
        z1 = FEM1%Elm_Loc(iel,igauss,3)
        r1 = sqrt((FEM1%Elm_Loc(iel,igauss,1) - Center_Rot(1))**2 + &
                  (FEM1%Elm_Loc(iel,igauss,2) - Center_Rot(2))**2) 

        ! loop over elements on surface 2
        ! to find the radial elements
        do jel = 1, FEM2%NEL

          Elm_Flag = 0

!!$          if (FEM2%PTYPE(jel) == 1) then

            do jgauss = 1, FEM2%NGAUSS

              z2 = FEM2%Elm_Loc(jel,jgauss,3)
              r2 = sqrt((FEM2%Elm_Loc(jel,jgauss,1) - Center_Rot(1))**2 + &
                        (FEM2%Elm_Loc(jel,jgauss,2) - Center_Rot(2))**2)

              he = FEM2%Elm_Loc(jel,jgauss,4)

!              if (abs(z1-z2) < ring_fact*FEM2%Elm_Size .and. &
!                  abs(r1-r2) < ring_fact*FEM2%Elm_Size) then
              if (abs(z1-z2) < ring_fact*he .and. &
                  abs(r1-r2) < ring_fact*he) then

                Elm_Flag = 1
              end if
            end do
!!$          end if

          if (Elm_Flag == 1) then
            num_rad_elm = num_rad_elm + 1
            FEM1%RAD_ELM_LIST(iel,igauss,num_rad_elm) = jel
          end if
        end do

        FEM1%RAD_ELM_NUM(iel,igauss) = num_rad_elm
      end do

!!$    end if

    if (mod(iel,1000)==0 .and. ismaster) write(*,*) "  finish FEM1 element", iel, FEM1%NEL
  end do
end subroutine f2f_find_elm_rad
