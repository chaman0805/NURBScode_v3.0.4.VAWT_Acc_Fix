subroutine nrb_find_elm(NRB, nsd, x1, ELM_DIST)

  use defs_shell
  implicit none

  type(mesh), intent(in)  :: NRB
  integer,    intent(in)  :: nsd           
  real(8),    intent(in)  :: x1(nsd)
  real(8),    intent(out) :: ELM_DIST(NRB%NEL)

  !  Local variables
  integer, parameter :: nvp = 5
  integer :: iel, i, j, ni, nj, p, q, nshl, nuk, nvk
  real(8) :: xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), &
             DetJb, nor(NSD), vp(nvp,2), dist(nvp)

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
  do iel = 1, NRB%NEL
 
    ! get NURB coordinates    
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)

    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't.
    ! Also, only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
         (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1)) .and. &
         (NRB%PTYPE(iel) /= 0  ) .or. (NRB%PTYPE(iel) /= 2  )) then

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel)
      
      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))
    
      do i = 1, nvp
        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0

        call eval_SHAPE_shell(vp(i,1), vp(i,2), shl, shgradl, shhessl, &
                              nor, xu, xd, dxdxi, ddxddxi, &
                              p, q, nsd, nshl, &
                              NRB%IEN(iel,1:nshl), NRB%NNODE,&
                              NRB%B_NET_U, NRB%B_NET_D, DetJb,&
                              ni, nj, nuk, nvk,&
                              NRB%U_KNOT(iel,1:nuk),&
                              NRB%V_KNOT(iel,1:nvk))

        dist(i) = sqrt(sum((xd-x1)**2))
      end do

      ! pick the shortest distance out of vertices or center
      ELM_DIST(iel) = minval(dist)

      deallocate(shl, shgradl, shhessl)

    end if
  end do    ! end loop elements
   
end subroutine nrb_find_elm
