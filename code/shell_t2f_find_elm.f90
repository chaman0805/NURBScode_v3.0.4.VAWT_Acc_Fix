!======================================================================
! We first check and save the 10% closest elements for finding closest
! points. Check the distance between element centers (gp=0).
!======================================================================
subroutine t2f_find_elm(FEM, TSP, BEZ, nsd, NEL_CLOSE, CLOSE_ELM)

  use defs_shell
  implicit none

  type(mesh), intent(in) :: FEM, TSP, BEZ
  integer,    intent(in) :: nsd, NEL_CLOSE

  !  Local variables
  integer :: iel, i, j, k, nshl, igauss
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), DetJb, nor(NSD), &
             gp(FEM%NGAUSS,2), gw(FEM%NGAUSS), Ginv(2,2)

  integer :: ELM_NUM(TSP%NEL), CLOSE_ELM(FEM%NEL,FEM%NGAUSS,NEL_CLOSE)
  real(8) :: ELM_DIST(TSP%NEL), ELM_DIST_TMP(TSP%NEL), ELM_DIST_CRI

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)

  ELM_NUM = -1

  ELM_DIST     = 1.0d6
  ELM_DIST_TMP = 1.0d6

  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  call genGPandGW_tri(gp, gw, FEM%NGAUSS)

  ! loop over elements
  do iel = 1, FEM%NEL

    if (FEM%PTYPE(iel) == 1) then

      nshl = FEM%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))

      do igauss = 1, FEM%NGAUSS

        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
        call eval_SHAPE_tri(gp(igauss,:), shl, shgradl, nor, &
                            xu, xd, dxdxi, nsd, nshl, &
                            FEM%IEN(iel,1:nshl), FEM%NNODE, &
                            FEM%B_NET_U, FEM%B_NET_D, DetJb)

        ! find the distance to vertices or center on the other surface
        call tsp_find_elm(TSP, BEZ, nsd, xd, ELM_DIST)

        ! Get the element size and set ELM_DIST_CRI for sorting 
        ! because you don't want to sort all the elements. 
        ! e.g., 10.0*maxval(sqrt(Ginv)) means only sort for the elements
        ! that have distance smaller than 10 times of the element size
        do j = 1, 2
          do i = 1, 2
            Ginv(i,j)  = sum(dxdxi(:,i)*dxdxi(:,j))
          end do
        end do
        ELM_DIST_CRI = 200.0d0*maxval(sqrt(Ginv))

        ! Only sort elements with distance smaller than criterion
        k = 0
        do i = 1, TSP%NEL
          if (ELM_DIST(i) < ELM_DIST_CRI) then
            k = k + 1
            ELM_DIST_TMP(k) = ELM_DIST(i)
            ELM_NUM(k)  = i
          end if
        end do

!!!        ! creat the element list (unsorted)
!!!        ELM_NUM = (/ (i, i=1,TSP%NEL) /)

        ! sort the distance
!!!        call in_sort(ELM_DIST, ELM_NUM, TSP%NEL)
        call in_sort(ELM_DIST_TMP(1:k), ELM_NUM(1:k), k)

        ! store the NEL_CLOSEth shortest distance elements
        CLOSE_ELM(iel,igauss,:) = ELM_NUM(1:NEL_CLOSE)
      end do
      deallocate(shl, shgradl, shhessl)
      if (mod(iel,10)==0) write(*,*) "  finish FEM element", iel, FEM%NEL
    end if
  end do    ! end loop elements
  
end subroutine t2f_find_elm
