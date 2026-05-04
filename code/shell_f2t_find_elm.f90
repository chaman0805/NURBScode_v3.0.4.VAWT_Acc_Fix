!======================================================================
! We first check and save the 10% closest elements for finding closest
! points. Check the distance between element centers (gp=0).
!======================================================================
subroutine f2t_find_elm(TSP, BEZ, FEM, nsd, NEL_CLOSE, CLOSE_ELM)

  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, FEM
  integer,    intent(in) :: nsd, NEL_CLOSE

  !  Local variables
  integer :: i, j, k, jel, iel, p, q, nshl, nshb, igauss, jgauss, ct
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), DetJb, nor(NSD), &
             gp(TSP%NGAUSS), gw(TSP%NGAUSS), Ginv(2,2)

  integer :: ELM_NUM(FEM%NEL), CLOSE_ELM(TSP%NEL,TSP%NGAUSS**2,NEL_CLOSE)
  real(8) :: ELM_DIST(FEM%NEL), ELM_DIST_TMP(FEM%NEL), ELM_DIST_CRI

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)

  ELM_NUM     = -1

  ELM_DIST     = 1.0d6
  ELM_DIST_TMP = 1.0d6

  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  call genGPandGW_shell(gp, gw, TSP%NGAUSS) 

  ! loop over all elements
  do iel = 1, TSP%NEL

    ! Only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( TSP%PTYPE(iel) == 1 ) then

      p = BEZ%P(iel); nshl = TSP%NSHL(iel)
      q = BEZ%Q(iel); nshb = BEZ%NSHL(iel)
      
      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))
 
      ! Loop over integration points (NGAUSS in each direction)
      ct = 0
      do jgauss = 1, TSP%NGAUSS
        do igauss = 1, TSP%NGAUSS
 
          ct = ct + 1

          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

          call eval_SHAPE_bez_sh(gp(igauss), gp(jgauss), &
                                 shl, shgradl, shhessl, &
                                 nor, xu, xd, dxdxi, ddxddxi, &
                                 p, q, nsd, nshl, nshb, &
                                 TSP%IEN(iel,1:nshl), TSP%NNODE, &
                                 TSP%B_NET_U, TSP%B_NET_D, DetJb, &
                                 BEZ%Ext(iel,1:nshl,1:nshb))

          ! find the distance to the other element vertices or centers
          call fem_find_elm(FEM, nsd, xd, ELM_DIST)

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

          ! ELM_DIST stores distance between S1 and S2 elements for
          ! all elements. However, we don't need to waste time sorting 
          ! all of them. We only want to sort those with distance 
          ! smaller than certain value.
          k = 0
          do i = 1, FEM%NEL
            if (ELM_DIST(i) < ELM_DIST_CRI) then
              k = k + 1
              ELM_DIST_TMP(k) = ELM_DIST(i)
              ELM_NUM(k) = i
            end if
          end do

!          if (k < NEL_CLOSE) then
!            write(*,*) "ERROR: k can't be smaller than NEL_CLOSE"
!            stop
!          end if

!!!          ! creat the element list (unsorted)
!!!          ELM_NUM = (/ (i, i=1,FEM%NEL) /)

          ! sort the distance
!!!          call in_sort(ELM_DIST, ELM_NUM, FEM%NEL)
          call in_sort(ELM_DIST_TMP(1:k), ELM_NUM(1:k), k)
        
          ! store the NEL_CLOSEth shortest distance elements
          ! If k < NEL_CLOSE, CLOSE_ELM(iel,ct,k+1:NEL_CLOSE) = -1
          CLOSE_ELM(iel,ct,:) = ELM_NUM(1:NEL_CLOSE)
        end do
      end do

      deallocate(shl, shgradl, shhessl)
      if (mod(iel,10)==0) write(*,*) "  finish TSP element", &
                                     iel, TSP%NEL
    end if
  end do    ! end loop elements
  
end subroutine f2t_find_elm
