subroutine tsp_find_elm(TSP, BEZ, nsd, x1, ELM_DIST)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: TSP, BEZ
  integer,    intent(in)  :: nsd           
  real(8),    intent(in)  :: x1(nsd)
  real(8),    intent(out) :: ELM_DIST(TSP%NEL)

  !  Local variables
  integer, parameter :: nvp = 5
  integer :: i, j, jel, iel, p, q, nshl, nshb
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), DetJb, nor(NSD), &
             vp(nvp,2), dist(nvp)

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)     

  ! vertices of the parametric domain
  vp(1,1) = -1.0d0; vp(1,2) = -1.0d0
  vp(2,1) =  1.0d0; vp(2,2) = -1.0d0
  vp(3,1) =  1.0d0; vp(3,2) =  1.0d0
  vp(4,1) = -1.0d0; vp(4,2) =  1.0d0
  vp(5,1) =  0.0d0; vp(5,2) =  0.0d0

  ! initialize with a large distance
  ELM_DIST = 1.0d6

  ! loop over elements
  do iel = 1, TSP%NEL
     
    ! Only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( TSP%PTYPE(iel) == 1 ) then

      p = BEZ%P(iel); nshl = TSP%NSHL(iel)
      q = BEZ%Q(iel); nshb = BEZ%NSHL(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))
  
      do i = 1, nvp
        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

        call eval_SHAPE_bez_sh(vp(i,1), vp(i,2), shl, shgradl, shhessl, &
                               nor, xu, xd, dxdxi, ddxddxi, &
                               p, q, nsd, nshl, nshb, &
                               TSP%IEN(iel,1:nshl), TSP%NNODE,&
                               TSP%B_NET_U, TSP%B_NET_D, DetJb,&
                               BEZ%Ext(iel,1:nshl,1:nshb))

        dist(i) = sqrt(sum((xd-x1)**2))
      end do

      ! distance between elements center
      ELM_DIST(iel) = minval(dist)

      deallocate(shl, shgradl, shhessl)
    end if
  end do    ! end loop elements
   
end subroutine tsp_find_elm
