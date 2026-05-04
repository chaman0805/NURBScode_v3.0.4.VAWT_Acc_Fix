!======================================================================
! We first check and save the 10% closest elements for finding closest
! points. Check the distance between element centers (gp=0).
!======================================================================
subroutine n2t_find_elm(TSP, BEZ, NRB, nsd, NEL_CLOSE, CLOSE_ELM)

  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, NRB
  integer,    intent(in) :: nsd, NEL_CLOSE

  !  Local variables
  integer :: i, j, k, jel, iel, p, q, nshl, nshb, igauss, jgauss, ct
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), DetJb, nor(NSD), &
             gp(TSP%NGAUSS), gw(TSP%NGAUSS)

  integer :: ELM_NUM(NRB%NEL), CLOSE_ELM(TSP%NEL,TSP%NGAUSS**2,NEL_CLOSE)
  real(8) :: ELM_DIST(NRB%NEL), ELM_DIST_TMP(NRB%NEL)

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)

  real(8), parameter :: ELM_DIST_CRI = 0.5d0

  ELM_NUM = -1

  ELM_DIST     = 1.0d6
  ELM_DIST_TMP = 1.0d6

  gp = 0.0d0; gw = 0.0d0; DetJb = 0.0d0
  call genGPandGW_shell(gp, gw, TSP%NGAUSS) 

  ! loop over elements
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
          call nrb_find_elm(NRB, nsd, xd, ELM_DIST)

          ! Only sort elements with distance smaller than criterion
          k = 0
          do i = 1, NRB%NEL
            if (ELM_DIST(i) < ELM_DIST_CRI) then
              k = k + 1
              ELM_DIST_TMP(k) = ELM_DIST(i)
              ELM_NUM(k)  = i
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

!!!          ! creat the element list (unsorted)
!!!          ELM_NUM = (/ (i, i=1,NRB%NEL) /)

          ! sort the distance
!!!          call in_sort(ELM_DIST, ELM_NUM, NRB%NEL)
          call in_sort(ELM_DIST_TMP(1:k), ELM_NUM(1:k), k)

          ! store the NEL_CLOSEth shortest distance elements
          CLOSE_ELM(iel,ct,:) = ELM_NUM(1:NEL_CLOSE)
        end do
      end do

      deallocate(shl, shgradl, shhessl)
      if (mod(iel,10)==0) write(*,*) "  finish TSP element", iel, TSP%NEL
    end if
  end do    ! end loop elements
  
end subroutine n2t_find_elm
