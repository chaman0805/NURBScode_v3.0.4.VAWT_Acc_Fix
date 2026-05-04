!======================================================================
! We first check and save the 10% closest elements for finding closest
! points. Check the distance between element centers (gp=0).
!======================================================================
subroutine t2n_find_elm(NRB, TSP, BEZ, nsd, NEL_CLOSE, CLOSE_ELM)

  use defs_shell
  implicit none

  type(mesh), intent(in) :: NRB, TSP, BEZ
  integer,    intent(in) :: nsd, NEL_CLOSE

  !  Local variables
  integer :: iel, i, j, k, ni, nj, p, q, nshl, nuk, nvk, igauss, jgauss, ct
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), DetJb, &
             nor(NSD), gp(NRB%NGAUSS), gw(NRB%NGAUSS)

  integer :: ELM_NUM(TSP%NEL), CLOSE_ELM(NRB%NEL,NRB%NGAUSS**2,NEL_CLOSE)
  real(8) :: ELM_DIST(TSP%NEL), ELM_DIST_TMP(TSP%NEL)

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)

  ! Assume that the distance between two elements of interest is 
  ! smaller than ELM_DIST_CRI
  ! For the wind turbine case, 0.2 did not give the right solution...
  real(8), parameter :: ELM_DIST_CRI = 0.5d0

  ELM_NUM     = -1

  ELM_DIST     = 1.0d6
  ELM_DIST_TMP = 1.0d6

  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 

  ! loop over elements
  do iel = 1, NRB%NEL

    ! get NURB coordinates
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)

    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't.
    ! Also, only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
         (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1)) .and. &
         NRB%PTYPE(iel) == 1 ) then

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))

      ct = 0
      do jgauss = 1, NRB%NGAUSS
        do igauss = 1, NRB%NGAUSS
          ct = ct + 1

          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

          ! Get Element Shape functions and their gradients
          call eval_SHAPE_shell(gp(igauss), gp(jgauss), &
                                shl, shgradl, shhessl, &
                                nor, xu, xd, dxdxi, ddxddxi, &
                                p, q, nsd, nshl, &
                                NRB%IEN(iel,1:nshl), NRB%NNODE,&
                                NRB%B_NET_U, NRB%B_NET_D, DetJb,&
                                ni, nj, nuk, nvk,&
                                NRB%U_KNOT(iel,1:nuk),&
                                NRB%V_KNOT(iel,1:nvk))

          ! find the distance to vertices or center on the other surface
          call tsp_find_elm(TSP, BEZ, nsd, xd, ELM_DIST)

          ! ELM_DIST stores distance between S1 and S2 elements for
          ! all elements. However, we don't need to waste time sorting 
          ! all of them. We only want to sort those with distance 
          ! smaller than certain value. Than find the minimum.
          k = 0
          do i = 1, TSP%NEL
            if (ELM_DIST(i) < ELM_DIST_CRI) then
              k = k + 1
              ELM_DIST_TMP(k) = ELM_DIST(i)
              ELM_NUM(k) = i
            end if
          end do

!          ! k is the number of elements that has distances smaller than 
!          ! the criterion. This number needs to be bigger than 
!          ! NEL_CLOSE, which is n% of the total number of elements.
!          if (k < NEL_CLOSE) then
!            write(*,*) "ERROR: k can't be smaller than NEL_CLOSE", k, NEL_CLOSE
!            write(*,*) "       need to reset ELM_DIST_CRI >", ELM_DIST_CRI
!            stop
!          end if

          ! The above procedure will save time for sorting
  
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
      if (mod(iel,10)==0) write(*,*) "  finish NRB element", iel, NRB%NEL
    end if
  end do    ! end loop elements
  
end subroutine t2n_find_elm
