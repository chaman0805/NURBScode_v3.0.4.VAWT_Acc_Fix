!======================================================================
! For every gauss point on surface 1, find the closest point on 
! surface 2. The point need to be between -1 and 1 in the parametric
! domain.
!======================================================================
subroutine nrb_find_point(NRB, nsd, x1, NEL_CLOSE, CLOSE_ELM, &
                          tkel, tp, tdist, tnewt, P_Flag)

  use defs_shell
  implicit none

  type(mesh), intent(in) :: NRB
  integer,    intent(in) :: nsd, NEL_CLOSE, CLOSE_ELM(NEL_CLOSE)
  real(8),    intent(in) :: x1(nsd)
  integer,    intent(out) :: tkel, tnewt, P_Flag
  real(8),    intent(out) :: tp(2), tdist

  integer, parameter :: nvp = 4
  real(8), parameter :: eps = 1.0d-13

  !  Local variables
  integer :: iel, kel, tel, i, j, p, q, nshl, inewt, nnewt
  integer :: nuk, nvk, ni, nj
  real(8) :: gp(2), dgp(2), DetJb, RhsNorm, &
             xu(nsd), xd(nsd), dxdxi(nsd,2), ddxddxi(nsd,3), nor(NSD), &
             Rhs(2), Lhs(2,2), iLHS(2,2), tmp,  &
             dist, xdist(nsd), vp(nvp,2), &
             tRhsNorm, tx1(nsd), txi(nsd)

  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)

  integer :: linewt, lnnewt
  real(8) :: lgp(2), ldgp, lRhs, lLhs, lxdist(nsd)

!  ! vertices of the parametric domain
!  vp(1,1) = -1.0d0; vp(1,2) = -1.0d0
!  vp(2,1) =  1.0d0; vp(2,2) = -1.0d0
!  vp(3,1) =  1.0d0; vp(3,2) =  1.0d0
!  vp(4,1) = -1.0d0; vp(4,2) =  1.0d0

  nnewt = 150      ! minimum newton (50 was not enough...)
                   ! there was a case took 82 iter to converge..
  lnnewt = 20

  tp    = 0.0d0    ! target point
  tdist = 1.0d6    ! initialize the distance
  tel   = 0
  tkel  = 0
  tnewt = 0
  txi   = 0.0d0
  tx1   = x1
  tRhsNorm = 0.0d0

  P_Flag = 0       ! if find the point, set to 1

!$$$      write(*,*) CLOSE_ELM
!$$$      read(*,*)

  ! loop over elements
  do kel = 1, NEL_CLOSE
!$$$      do kel = 11, 11
!$$$      do iel = 1, NRB%NEL
!$$$      do iel = 2998, 2998
     
    ! get the real element number
    iel = CLOSE_ELM(kel)
    if (iel == -1) exit

    ! get NURB coordinates
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)

    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't.
    ! Also, only loop through the elements on the blade surface
    ! (e.g., ptype == 1)
    if ( (NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
         (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1)) .and. &
         (NRB%PTYPE(iel) /= 0  ) .or. (NRB%PTYPE(iel) /= 2  ) ) then

      p     = NRB%P(iel)
      q     = NRB%Q(iel)
      nshl  = NRB%NSHL(iel)
      nuk   = NRB%NUK(iel)
      nvk   = NRB%NVK(iel)

      allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3))
     
      !--------------------------------------------------------
      ! Use Newton's method to find the perpendicular points.
      ! Only consider those that are located between -1 and 1.
      ! To make sure they are the shortest distance, also need
      ! to check with vertex points.
      !--------------------------------------------------------
      gp   = 0.0d0    ! initial guess
      dist = 0.0d0
      do inewt = 1, nnewt
          
        dgp = 0.0d0
      
        ! Get Element Shape functions and their gradients
        shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
        xu = 0.0d0; xd = 0.0d0; dxdxi= 0.0d0; ddxddxi= 0.0d0 

        call eval_SHAPE_shell(gp(1), gp(2), &
                              shl, shgradl, shhessl, &
                              nor, xu, xd, dxdxi, ddxddxi, &
                              p, q, nsd, nshl, &
                              NRB%IEN(iel,1:nshl), NRB%NNODE,&
                              NRB%B_NET_U, NRB%B_NET_D, DetJb,&
                              ni, nj, nuk, nvk,&
                              NRB%U_KNOT(iel,1:nuk), NRB%V_KNOT(iel,1:nvk))
                              

        !-Newton------------------------------------
        Rhs = 0.0d0; Lhs = 0.0d0; iLhs = 0.0d0

        xdist = xd-x1

        Rhs(1) = sum(xdist*dxdxi(:,1))
        Rhs(2) = sum(xdist*dxdxi(:,2))

        RhsNorm = sqrt(sum(Rhs**2))
        if (isnan(RhsNorm) .eqv. .true.) then
          exit
 
        else if (RhsNorm >= eps) then

          ! does not work without second derivatives
          Lhs(1,1) = sum(dxdxi(:,1)*dxdxi(:,1)) + sum(xdist*ddxddxi(:,1))
          Lhs(1,2) = sum(dxdxi(:,1)*dxdxi(:,2)) + sum(xdist*ddxddxi(:,2))
          Lhs(2,1) = sum(dxdxi(:,2)*dxdxi(:,1)) + sum(xdist*ddxddxi(:,2))
          Lhs(2,2) = sum(dxdxi(:,2)*dxdxi(:,2)) + sum(xdist*ddxddxi(:,3))

          tmp = Lhs(1,1)*Lhs(2,2)-Lhs(1,2)*Lhs(2,1)

          iLhs(1,1) =  Lhs(2,2)/tmp
          iLhs(1,2) = -Lhs(1,2)/tmp
          iLhs(2,1) = -Lhs(2,1)/tmp
          iLhs(2,2) =  Lhs(1,1)/tmp

          dgp = matmul(iLHS,-Rhs)
          gp  = gp + dgp
          !--------------------------------------------

!$$$!!!            ! write out the convergence history
!$$$!!!            if (abs(gp(1))<=1.0d0 .and. abs(gp(2))<=1.0d0) then
!$$$              write(*,*) 'iel, inewt, Rhs=', iel, inewt, RhsNorm
!$$$              write(*,*) gp, sqrt(sum(xdist*xdist))
!$$$!!!            end if

        ! If converged:
        else if (RhsNorm < eps) then

          !-----------------------------------------------------
          ! if gp is located outside of -1 and 1, then find
          ! the shortest distance to the vertices of the 
          ! parametric domain. [-1,-1] -> [1,1]
          !-----------------------------------------------------
          if (abs(gp(1))>1.0d0 .and. abs(gp(2))>1.0d0) then

            lgp(1) = sign(1.0d0,gp(1))
            lgp(2) = sign(1.0d0,gp(2))

!            do i = 1, nvp

              shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
              xu = 0.0d0; xd = 0.0d0; dxdxi= 0.0d0; ddxddxi= 0.0d0 

!              call eval_SHAPE_shell(vp(i,1), vp(i,2), &
              call eval_SHAPE_shell(lgp(1), lgp(2), &
                                    shl, shgradl, shhessl,  &
                                    nor, xu, xd, dxdxi, ddxddxi, &
                                    p, q, nsd, nshl, &
                                    NRB%IEN(iel,1:nshl), NRB%NNODE, &
                                    NRB%B_NET_U, NRB%B_NET_D, DetJb, &
                                    ni, nj, nuk, nvk, &
                                    NRB%U_KNOT(iel,1:nuk), &
                                    NRB%V_KNOT(iel,1:nvk))

              dist = sqrt(sum((xd-x1)**2))

              ! store the shortest distance at the vertices
              if (dist < tdist) then 
                ! record the target results
                tdist = dist   
!                tp    = vp(i,:)
                tp    = lgp
                tel   = iel
                tkel  = kel
                txi   = xd

                P_Flag = 1
              end if
!            end do

            exit    ! exit the newton loop

          !--------------------------------------------------------
          ! Add this block: if we can't find the solution between
          ! -1 and 1, the shortest distance is the corner in 1D.
          ! But in 2D, this is not true. Maybe at the edges, but
          ! definitely not the vertices.
          ! So what we do here is check the converged point, and 
          ! if one of them is between -1 and 1, then bring the 
          ! other one (outside of -1 and 1) to -1 or 1 depends on 
          ! the sign. Then do newton to find the perpendicular 
          ! point subject to the constrain
          !--------------------------------------------------------
          else if (abs(gp(1))> 1.0d0 .and. abs(gp(2))<=1.0d0) then

            ! initial guess
            lgp(1) = sign(1.0d0,gp(1))
            lgp(2) = gp(2)
            do linewt = 1, lnnewt
          
              ldgp = 0.0d0
     
              shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
              xu = 0.0d0; xd = 0.0d0; dxdxi= 0.0d0; ddxddxi= 0.0d0 

              call eval_SHAPE_shell(lgp(1), lgp(2), &
                              shl, shgradl, shhessl, &
                              nor, xu, xd, dxdxi, ddxddxi,& 
                              p, q, nsd, nshl, &
                              NRB%IEN(iel,1:nshl), NRB%NNODE,&
                              NRB%B_NET_U, NRB%B_NET_D, DetJb,&
                              ni, nj, nuk, nvk,&
                              NRB%U_KNOT(iel,1:nuk),&
                              NRB%V_KNOT(iel,1:nvk))

              lxdist = xd-x1

              lRhs = 0.0d0; lLhs = 0.0d0

              lRhs = sum(lxdist*dxdxi(:,2))

              if (abs(lRhs) >= eps) then

                lLhs = sum(dxdxi(:,2)*dxdxi(:,2)) + sum(lxdist*ddxddxi(:,3))

                ldgp = -lRhs/lLhs
                lgp(2) = lgp(2) + ldgp

              else !if (abs(lRhs) < eps) then
                ! check if it goes out of bound again...
                if (abs(lgp(2))<=1.0d0) then
                  dist = sqrt(sum(lxdist**2))

                  if (dist < tdist) then 
                    ! update the target results
                    tdist    = dist
                    tp       = lgp
                    tel      = iel
                    tkel     = kel
                    tnewt    = inewt
                    tRhsNorm = abs(lRhs)
                    txi      = xd

                    P_Flag = 1
                  end if
                  exit   ! exit the inner newton loop
                else
                  exit
                end if
              end if
            end do
            exit     ! exit the outer newton loop

          else if (abs(gp(1))<=1.0d0 .and. abs(gp(2))>1.0d0) then

            ! initial guess
            lgp(1) = gp(1)
            lgp(2) = sign(1.0d0,gp(2))
            do linewt = 1, lnnewt
      
              ldgp = 0.0d0
 
              shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
              xu = 0.0d0; xd = 0.0d0; dxdxi= 0.0d0; ddxddxi= 0.0d0 

              call eval_SHAPE_shell(lgp(1), lgp(2), &
                                    shl, shgradl, shhessl, &
                                    nor, xu, xd, dxdxi, ddxddxi, &
                                    p, q, nsd, nshl, & 
                                    NRB%IEN(iel,1:nshl), NRB%NNODE, &
                                    NRB%B_NET_U, NRB%B_NET_D, DetJb, &
                                    ni, nj, nuk, nvk, &
                                    NRB%U_KNOT(iel,1:nuk), &
                                    NRB%V_KNOT(iel,1:nvk))

              lxdist = xd-x1

              lRhs = 0.0d0; lLhs = 0.0d0

              lRhs = sum(lxdist*dxdxi(:,1))

              if (abs(lRhs) >= eps) then

                lLhs = sum(dxdxi(:,1)*dxdxi(:,1)) + sum(lxdist*ddxddxi(:,1))

                ldgp = -lRhs/lLhs
                lgp(1) = lgp(1) + ldgp

              else !if (abs(lRhs) < eps) then
                ! check if it goes out of bound again...
                if (abs(lgp(1)) <= 1.0d0) then
                  dist = sqrt(sum(lxdist**2))

                  if (dist < tdist) then 
                    ! update the target results
                    tdist    = dist
                    tp       = lgp
                    tel      = iel
                    tkel     = kel
                    tnewt    = inewt
                    tRhsNorm = abs(lRhs)
                    txi      = xd

                    P_Flag = 1
                  end if
                  exit   ! exit the inner newton loop
                else
                  exit
                end if
              end if

            end do   ! inner newton loop
            exit     ! exit the outer newton loop

          !--------------------------------------------------------
          ! This is what we originally want. Find the normal point
          ! that is between -1 and 1.
          !--------------------------------------------------------
          else if (abs(gp(1))<=1.0d0 .and. abs(gp(2))<=1.0d0) then

!$$$                ! converge better
!$$$                if (RhsNorm < 1.0d-12) then

              dist = sqrt(sum(xdist**2))

              ! there could be multiple solution
              ! only recorded the smallest distance
              if (dist < tdist) then 
                ! update the target results
                tdist    = dist
                tp       = gp
                tel      = iel
                tkel     = kel
                tnewt    = inewt
                tRhsNorm = RhsNorm
                txi      = xd

                P_Flag   = 1
              end if
              exit   ! exit the newton loop
!$$$                end if

          else
            write(*,*) "ERROR: Undefined point range"
            stop
          end if

        end if

!$$$           if (inewt == nnewt) then
!$$$             write(*,*) "something is wrong...", kel, iel
!$$$             stop
!$$$           end if

      end do  ! end newton loop

      deallocate(shl, shgradl, shhessl)

    end if
  end do

  write(*,*) 'iel, kel, inewt=', tel, tkel, tnewt
  write(*,*) 'Rhs =', tRhsNorm
  write(*,*) 'gp  =', tp(1), tp(2)
  write(*,*) 'loc1=', tx1
  write(*,*) 'loc2=', txi
  write(*,*) 'dist=', tdist

  if (tdist >= 1.0d-1) then
    write(*,'(45("%"))')
    write(*,*) "!!WARNNING!! Large Distance:", tdist
    write(*,'(45("%"))')
  end if

end subroutine nrb_find_point
