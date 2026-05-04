!-----------------------------------------------------------------------
! This subroutine modifies the local stiffness matrix and load vector
! to account for constrained nodes. 
!-----------------------------------------------------------------------
subroutine BCLhs_3D_shell(nsd, nshl, nnode, IEN, IBC, xKebe)
  implicit none  
  integer, intent(in)    :: nsd, nshl, nnode, IEN(nshl), IBC(nnode,nsd)
  real(8), intent(inout) :: xKebe(NSD*NSD,NSHL,NSHL)
  integer :: aa, bb, cc
  
  ! loop through local nodes 
  do aa = 1, NSHL
    cc = IEN(aa)    ! global node number of local node aa
    ! Check status of constraint in x direction
    if (IBC(cc,1) == 1) then
      ! update matrix to account for constraint
      do bb = 1, NSHL   ! project away row (X-dir)     
        xKebe(1,aa,bb) = 0.0d0
        xKebe(2,aa,bb) = 0.0d0
        xKebe(3,aa,bb) = 0.0d0
      end do 
      do bb = 1, NSHL   !project away column (X-dir)
        xKebe(1,bb,aa) = 0.0d0
        xKebe(4,bb,aa) = 0.0d0
        xKebe(7,bb,aa) = 0.0d0 
      end do
      xKebe(1,aa,aa) = 1.0d0
    end if
      
    ! Check for constraint in y direction
    if (IBC(cc,2)==1) then
      do bb = 1, NSHL   ! project away row (Y-dir)     
        xKebe(4,aa,bb) = 0.0d0
        xKebe(5,aa,bb) = 0.0d0
        xKebe(6,aa,bb) = 0.0d0
      end do
      do bb = 1, NSHL   !project away column (Y-dir)
        xKebe(2,bb,aa) = 0.0d0
        xKebe(5,bb,aa) = 0.0d0
        xKebe(8,bb,aa) = 0.0d0
      end do
      xKebe(5,aa,aa) = 1.0d0       
    end if
    
    ! Check for constraint in z direction
    if (IBC(cc,3)==1) then
      do bb = 1, NSHL   ! project away row (Y-dir)     
        xKebe(7,aa,bb) = 0.0d0
        xKebe(8,aa,bb) = 0.0d0
        xKebe(9,aa,bb) = 0.0d0
      end do
      do bb = 1, NSHL   !project away column (Y-dir)
        xKebe(3,bb,aa) = 0.0d0
        xKebe(6,bb,aa) = 0.0d0
        xKebe(9,bb,aa) = 0.0d0
      end do
      xKebe(9,aa,aa) = 1.0d0       
    end if     
  end do
end subroutine BCLhs_3D_shell


!---------------------------------------------------
! This subroutine modifies the local load vector
! to account for constrained nodes. 
!---------------------------------------------------
subroutine BCRhs_3D_shell(nsd, nshl, nnode, IEN, IBC, Rhs)  
  implicit none  
  integer, intent(in)    :: nsd, nshl, nnode, IEN(nshl), IBC(nnode,nsd)
  real(8), intent(inout) :: Rhs(NSD,NSHL)
  integer :: aa, cc, ii
  
  ! loop through local nodes 
  do aa = 1, nshl
    ! global node number of local node aa
    cc = IEN(aa)
    ! Check status of constraint in x,y,z direction
    do ii = 1, nsd
      if (IBC(cc,ii) == 1) then 
    ! update load vector to account for constraint
        Rhs(ii,aa) = 0.0d0
      end if
    end do    
  end do
end subroutine BCRhs_3D_shell
