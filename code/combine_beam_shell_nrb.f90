subroutine combine_beam_shell_nrb(NSD, NRB_BEA, NRB_SHL)

  use types_beam
  use defs_shell
  implicit none

  type(mesh_beam),    intent(inout)    :: NRB_BEA
  type(mesh),    intent(inout)    :: NRB_SHL
  real(8), allocatable :: B_NET_TMP(:,:)
  integer,     intent(in)    :: NSD

  integer :: ier, i, ip, p, mcp, nnode, nel, nshl, eloc, eglob
  
  allocate(B_NET_TMP(NRB_SHL%NNODE,NSD+1))
  B_NET_TMP = 0.0d0
  B_NET_TMP = NRB_SHL%B_NET

  do i = 1, NRB_SHL%NNODE
    if((abs(NRB_SHL%B_NET(i,1)) <= 1.0d-12).and.(abs(NRB_SHL%B_NET(i,2)) <= 1.0d-12) &
      .and.(abs(NRB_SHL%B_NET(i,3)) <= 1.0d-12)) then
      B_NET_TMP(i,:) = NRB_BEA%B_NET(i,:)    
    end if
  end do

  NRB_BEA%B_NET = B_NET_TMP
  NRB_SHL%B_NET = B_NET_TMP

  !write (*,*) 'NRB_SHL%NNODE=', NRB_SHL%NNODE
!A  do i = 1, NRB_SHL%NNODE
!A    if (abs(NRB_SHL%B_NET(i,1)-NRB_BEA%B_NET(i,1))<=1.0d-12 .and. abs(NRB_SHL%B_NET(i,2)-NRB_BEA%B_NET(i,2))<=1.0d-12) then
!A      if (abs(NRB_SHL%B_NET(i,3)-NRB_BEA%B_NET(i,3))<=1.0d-12) then
        !write (*,*) 'case 1'
!A      endif
!A    elseif (abs(NRB_SHL%B_NET(i,1)-0.0d0)<=1.0d-12 .and. abs(NRB_SHL%B_NET(i,2)-0.0d0)<=1.0d-12) then
!A      if (abs(NRB_SHL%B_NET(i,3)-0.0d0)<=1.0d-12) then
        !write (*,*) 'case 2'
!A        NRB_SHL%B_NET(i,:)=NRB_BEA%B_NET(i,:)
!A      endif
!A    else
      !write (*,*) 'case 3'
!A        NRB_BEA%B_NET(i,:)=NRB_SHL%B_NET(i,:)
!A      end if
!A    endif
!A  end do

!A  NRB_BEA%B_NET_U = NRB_BEA%B_NET
!A  NRB_BEA%B_NET_D = NRB_BEA%B_NET
!A  NRB_BEA%B_NET_D_old = NRB_BEA%B_NET

!A  NRB_SHL%B_NET_U = NRB_SHL%B_NET
!A  NRB_SHL%B_NET_D = NRB_SHL%B_NET
!A  NRB_SHL%B_NET_D_old = NRB_SHL%B_NET

  ! use the coordinates to setup IBC_SH because one also
  ! need to take the bending strips into consideration
  write(*,*) NRB_BEA%NNODE, NRB_SHL%NNODE
  do i = 1, NRB_BEA%NNODE
    if(NRB_BEA%B_NET(i,1) - NRB_SHL%B_NET(i,1) .ne. 0.0d0) write(*,*) 'BOOMX', NRB_BEA%B_NET(i,1) - NRB_SHL%B_NET(i,1), &
    NRB_BEA%B_NET(i,1:3), NRB_SHL%B_NET(i,1:3)
    if(NRB_BEA%B_NET(i,2) - NRB_SHL%B_NET(i,2) .ne. 0.0d0) write(*,*) 'BOOMY', NRB_BEA%B_NET(i,2) - NRB_SHL%B_NET(i,2), &
    NRB_BEA%B_NET(i,1:3), NRB_SHL%B_NET(i,1:3)
    if(NRB_BEA%B_NET(i,3) - NRB_SHL%B_NET(i,3) .ne. 0.0d0) write(*,*) 'BOOMZ', NRB_BEA%B_NET(i,3) - NRB_SHL%B_NET(i,3), &
    NRB_BEA%B_NET(i,1:3), NRB_SHL%B_NET(i,1:3)
    ! root
!A    write(*,*) NRB_SHL%B_NET(i,1:3)
    if ( (abs(NRB_BEA%B_NET(i,3)-0.0d0)<=1.0d-12).or.(abs(NRB_BEA%B_NET(i,3)-0.2d0)<=1.0d-12)) then 
    !original radius is 3.25
!A      write(*,*) NRB_BEA%B_NET(i,1:3)
      NRB_SHL%IBC(i,:) = 1
      NRB_BEA%IBC(i,:) = 1

      !write (*,*) 'root node =', i
    end if
  end do

  deallocate(B_NET_TMP)



end subroutine combine_beam_shell_nrb
