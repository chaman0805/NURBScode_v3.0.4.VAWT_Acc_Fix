!======================================================================
!
!======================================================================
subroutine setBCs_CFD()
  use commonvars
  use aAdjkeep
  use mpi
  
  implicit none
  
  integer :: b, i, j, k, n, d, dir, ptmp, iface
  real(8) :: u(NSD), phi

  IBC  = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0

  do b = 1, NBOUND

    do d = 1, NSD
      if (BCugType(bound(b)%FACE_ID,d) == 1) then
        do j = 1, bound(b)%NNODE
          k = bound(b)%BNODES(j)
          IBC(k,d) = 1
          ug (k,d) = BCugValu(b,d)
          acg(k,d) = ( (ug(k,d)-ugold(k,d))*Dtgl  &
                        + (gami-1.0d0)*acgold(k,d)  )/gami
        end do

      else if (BCugType(bound(b)%FACE_ID,d) == 3) then
        do j = 1, bound(b)%NNODE
          k = bound(b)%BNODES(j)
          IBC(k,d) = 1
          ug (k,d) = ugm(k,d)
          acg(k,d) = ( (ug(k,d)-ugold(k,d))*Dtgl  &
                      + (gami-1.0d0)*acgold(k,d)  )/gami
        end do
      end if
    end do
    
  end do  ! end loop: all faces

!  ! Stationary Case
!  IBC(:,4:6) = 1        

!  ! 2D problem
!  IBC (:,3) = 1
!  IBC (:,6) = 1
!  ug  (:,3) = 0.0d0
!  acg (:,3) = 0.0d0
!  ugm (:,3) = 0.0d0
!  acgm(:,3) = 0.0d0
!  dg  (:,3) = 0.0d0

end subroutine setBCs_CFD




!======================================================================
! Get the rotation matrices
!======================================================================
subroutine get_Rmat(theta, thetd, thedd, Rmat, Rdot, Rddt)
  implicit none

  real(8), intent(in ) :: theta, thetd, thedd
  real(8), intent(out) :: Rmat(3,3), Rdot(3,3), Rddt(3,3)

  Rmat      = 0.0d0
  Rmat(2,2) = cos(theta) - 1.0d0
  Rmat(2,3) =-sin(theta)
  Rmat(3,2) = sin(theta)
  Rmat(3,3) = cos(theta) - 1.0d0

  Rdot      = 0.0d0
  Rdot(2,2) =-sin(theta)*thetd
  Rdot(2,3) =-cos(theta)*thetd
  Rdot(3,2) = cos(theta)*thetd
  Rdot(3,3) =-sin(theta)*thetd

  Rddt      = 0.0d0
  Rddt(2,2) =-cos(theta)*thetd**2 - sin(theta)*thedd
  Rddt(2,3) = sin(theta)*thetd**2 - cos(theta)*thedd
  Rddt(3,2) =-sin(theta)*thetd**2 + cos(theta)*thedd
  Rddt(3,3) =-cos(theta)*thetd**2 - sin(theta)*thedd
end subroutine get_Rmat



!======================================================================
! Get the rotation matrices
!======================================================================
subroutine get_Rmat_Z(theta, thetd, thedd, Rmat, Rdot, Rddt)
  implicit none

  real(8), intent(in ) :: theta, thetd, thedd
  real(8), intent(out) :: Rmat(3,3), Rdot(3,3), Rddt(3,3)

  Rmat      = 0.0d0
  Rmat(1,1) = cos(theta) - 1.0d0
  Rmat(1,2) =-sin(theta)
  Rmat(2,1) = sin(theta)
  Rmat(2,2) = cos(theta) - 1.0d0

  Rdot      = 0.0d0
  Rdot(1,1) =-sin(theta)*thetd
  Rdot(1,2) =-cos(theta)*thetd
  Rdot(2,1) = cos(theta)*thetd
  Rdot(2,2) =-sin(theta)*thetd

  Rddt      = 0.0d0
  Rddt(1,1) =-cos(theta)*thetd**2 - sin(theta)*thedd
  Rddt(1,2) = sin(theta)*thetd**2 - cos(theta)*thedd
  Rddt(2,1) =-sin(theta)*thetd**2 + cos(theta)*thedd
  Rddt(2,2) =-cos(theta)*thetd**2 - sin(theta)*thedd
end subroutine get_Rmat_Z




!======================================================================
! Set the BC for the case of rotating HAWT
!======================================================================
subroutine setMeshBCs_hawt
  use commonvars
  use aAdjkeep
  implicit none
  
  integer :: b, i, j, k, dd
  real(8) :: rot(NSD)

  do b = 1, NBOUND
    ! apply the same mesh BC for inlet/side/outlet/rotor
    do i = 1, bound(b)%NNODE
      j = bound(b)%BNODES(i)     
      IBC(j,4:6) = 1
!A-------------------------------------S 05/11/13 Apply rotation to
!A the surface of rotating cylinder corresponding to the structural
!A rotation computed from Rt
        do dd = 1, NSD
          dg(j,dd) = sum((Rmat(dd,:) - Identity(dd,:))*(xg(j,1:3) - &
                     xrotOld(1:3))) + (xrot(dd) - xrotOld(dd))

          acgm(j,dd) = (((dg(j,dd)-dgold(j,dd)-Delt*ugmold(j,dd)))* &
                       2.0d0/(Delt*Delt) - &
                       (1.0d0-2.0d0*beti)*acgmold(j,dd))/(2.0d0*beti)

          ugm(j,dd) = ugmold(j,dd) + &
                       Delt*((1.0d0-gami)*acgmold(j,dd)+gami*acg(j,dd))
        end do
!A-------------------------------------E 
    end do      
  end do  
end subroutine setMeshBCs_hawt




!======================================================================
! Set the BC for the case of rotating 5MW with tower
!======================================================================
subroutine setMeshBCs_tower
  use commonvars
  use aAdjkeep
  use defs_shell
  implicit none
  
  integer :: b, i, j, k, dd

  ! for the fixed outer box domain (including the tower)
!A  do i = 1, NNODE
!A    if (NodeID(i) == 31) then   
!A      IBC(i,4:6) = 1 
!A    end if
!A  end do

  ! for the disk and the blade
  do b = 1, NBOUND

    if ( (bound(b)%FACE_ID == 11) ) then

      do i = 1, bound(b)%NNODE
        j = bound(b)%BNODES(i)     
        IBC(j,4:6) = 1
!A-------------------------------------S 05/11/13 Apply rotation to
!A the surface of rotating cylinder corresponding to the structural
!A rotation computed from Rt
        do dd = 1, NSD
          dg(j,dd) = sum((Rmat(dd,:) - Identity(dd,:))*(xg(j,1:3) - &
                     xrotOld(1:3))) + (xrot(dd) - xrotOld(dd))

          acgm(j,dd) = (((dg(j,dd)-dgold(j,dd)-Delt*ugmold(j,dd)))* &
                       2.0d0/(Delt*Delt) - &
                       (1.0d0-2.0d0*beti)*acgmold(j,dd))/(2.0d0*beti)

          ugm(j,dd) = ugmold(j,dd) + &
                       Delt*((1.0d0-gami)*acgmold(j,dd)+gami*acg(j,dd))

        end do
!A-------------------------------------E          
      end do
!A fixed the boundary of outer domain
    elseif  ( (bound(b)%FACE_ID == 1) .or. (bound(b)%FACE_ID == 2) .or. (bound(b)%FACE_ID == 3) ) then  
      do i = 1, bound(b)%NNODE
        j = bound(b)%BNODES(i)     !A Points to volumetric node from surface node
        IBC(j,4:6) = 1   
          dg(j,1:3) = 0.0d0
          acgm(j,1:3) = 0.0d0
          ugm(j,1:3) = 0.0d0
      end do

    elseif  ( (bound(b)%FACE_ID == 21) ) then  
      do i = 1, bound(b)%NNODE
        j = bound(b)%BNODES(i)     !mesh BC on a rotor
        IBC(j,4:6) = 1   
      end do


    elseif ( (bound(b)%FACE_ID == 12) ) then

      do i = 1, bound(b)%NNODE
        j = bound(b)%BNODES(i)     
        IBC(j,4:6) = 1
        do dd = 1, NSD
          dg(j,dd) = sum((Rtang(dd,:) - Identity(dd,:))*(xg(j,1:3) - &
                     xrotOld(1:3))) + (xrot(dd) - xrotOld(dd))

          acgm(j,dd) = (((dg(j,dd)-dgold(j,dd)-Delt*ugmold(j,dd)))* &
                       2.0d0/(Delt*Delt) - &
                       (1.0d0-2.0d0*beti)*acgmold(j,dd))/(2.0d0*beti)

          ugm(j,dd) = ugmold(j,dd) + &
                       Delt*((1.0d0-gami)*acgmold(j,dd)+gami*acg(j,dd))
        end do
      end do
    end if
  end do     
end subroutine setMeshBCs_tower


!======================================================================
! Set the BC for the case of rotating 5MW with tower
!======================================================================
subroutine setNMBCs(NM)
  use commonvars
  use mpi
  use aAdjkeep
  use defs_shell
  implicit none
  
  type(shell_nmb), intent(inout) :: NM
  integer :: b, i, j, k, dd

  ! for the fixed outer box domain (including the tower)
!A  do i = 1, NNODE
!A    if (NodeID(i) == 31) then   
!A      IBC(i,4:6) = 1 
!A    end if
!A  end do

  ! for the disk and the blade




      do i = 1, NM%FEM(1)%NNODE
!        j = bound(b)%BNODES(i)     
!A-------------------------------------S 05/11/13 Apply rotation to
!A the surface of rotating cylinder corresponding to the structural
!A rotation computed from Rt
        do dd = 1, NSD

          NM%FEM(1)%dsh(i,dd) = sum((Rmat(dd,:) - Identity(dd,:))*(NM%FEM(1)%B_NET(i,1:3) - &
                     xrotOld(1:3))) + (xrot(dd) - xrotOld(dd))

          NM%FEM(1)%ash(i,dd) = (((NM%FEM(1)%dsh(i,dd)-NM%FEM(1)%dshold(i,dd)-Delt*NM%FEM(1)%ushold(i,dd)))* &
                       2.0d0/(Delt*Delt) - &
                       (1.0d0-2.0d0*beti)*NM%FEM(1)%ashold(i,dd))/(2.0d0*beti)

          NM%FEM(1)%ush(i,dd) = NM%FEM(1)%ushold(i,dd) + &
                       Delt*((1.0d0-gami)*NM%FEM(1)%ashold(i,dd)+gami*NM%FEM(1)%ash(i,dd))
        end do
!A-------------------------------------E          
      end do
!A fixed the boundary of outer domain
      do i = 1, NM%FEM(2)%NNODE
        do dd = 1, NSD
          NM%FEM(2)%dsh(i,dd)   = sum((Rtang(dd,:) - Identity(dd,:))*(NM%FEM(2)%B_NET(i,1:3) - &
                     xrotOld(1:3))) + (xrot(dd) - xrotOld(dd))

        
          NM%FEM(2)%ash(i,dd) = (((NM%FEM(2)%dsh(i,dd)-NM%FEM(2)%dshold(i,dd)-Delt*NM%FEM(2)%ushold(i,dd)))* &
                       2.0d0/(Delt*Delt) - &
                       (1.0d0-2.0d0*beti)*NM%FEM(2)%ashold(i,dd))/(2.0d0*beti)

          NM%FEM(2)%ush(i,dd) = NM%FEM(2)%ushold(i,dd) + &
                       Delt*((1.0d0-gami)*NM%FEM(2)%ashold(i,dd)+gami*NM%FEM(2)%ash(i,dd))

        end do
      end do
   
end subroutine setNMBCs




!======================================================================
!
!======================================================================
subroutine getWave(u, phi, x, t) 
  use commonvars
  use aAdjkeep
  use mpi
  implicit none
  
  real(8) :: u(NSD), phi, x(NSD), t
  real(8) :: kxwt, kzh, amp, xy
  
  amp = wave_periode*(wave_amp/sinh(wave_length*water_depth))  
  xy = dcos(wave_angle)*x(1) + dsin(wave_angle)*x(2)

  kxwt = wave_length*(xy - domain_left - Uin*t) - wave_periode*time
  kzh  = wave_length*min(x(3)-domain_bottom,water_depth)

  u(1) = amp*dcosh(kzh)*dcos(kxwt)*dcos(wave_angle)  + Uin
  u(2) = amp*dcosh(kzh)*dcos(kxwt)*dsin(wave_angle) 
  u(3) = amp*dsinh(kzh)*dsin(kxwt)  
  phi  = wave_amp*dcos(kxwt) + water_level - x(3) 
end subroutine getWave




!======================================================================
!
!======================================================================
subroutine setBCs()
  use commonvars
  use aAdjkeep
  use mpi
  
  implicit none
  
  integer :: b, i, j, k, n, dir, ptmp
  real(8) :: u(NSD), phi

  IBC  = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0
  do b = 1, NBOUND
    if ((bound(b)%FACE_ID >= 7).and.(BCType(7) == 1)) then ! Hull 
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)
	ug    (i,:) = ugm    (i,:)
	ugold (i,:) = ugmold (i,:)
	acg   (i,:) = acgm   (i,:)
	acgold(i,:) = acgmold(i,:)
	   
        IBC(i,1:3) = 1
      end do
    end if 
      
    ! Sides  	  
    if ( ((bound(b)%FACE_ID == 2).and.(BCType(2) == 1)).or. &
         ((bound(b)%FACE_ID == 4).and.(BCType(4) == 1)) ) then 
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)   
	call getWave(u, phi, xg(i,:),time)
		 
	if ( ((bound(b)%FACE_ID == 2).and.(u(2) > 0.0d0)).or. & 
             ((bound(b)%FACE_ID == 4).and.(u(2) < 0.0d0))  ) then
	  ug  (i,:) = u  
	  phig(i)   = phi  
	  
	  acg(i,:) = ( (ug(i,:)-ugold(i,:))*Dtgl &
        	       + (gami-1.0d0)*acgold(i,:)  )/gami
		  
	  rphig(i) = ( (phig(i) - phigold(i))*Dtgl & 
                       + (gami-1.0d0)*rphigold(i)  )/gami
	 
	  IBC(i,1:3) = 1
	  IBC(i,8)   = 1
        end if
      end do
    end if
	
    ! inflow 
    if ((bound(b)%FACE_ID == 5).and.(BCType(5) == 1)) then 
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)     
	call getWave(ug(i,:), phig(i), xg(i,:), time)
	  	  
	acg(i,:)= ( (ug(i,:)-ugold(i,:))*Dtgl &
                    + (gami-1.0d0)*acgold(i,:)  )/gami
		  
        rphig(i) = ( (phig(i) - phigold(i))*Dtgl &
                     + (gami-1.0d0)*rphigold(i)  )/gami
	 
        IBC(i,1:3) = 1
        IBC(i,8)   = 1
      end do 
    endif
	
    ! Bottom 
    if ((bound(b)%FACE_ID == 1).and.(BCType(1) == 1) ) then  
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j) 	   
	ug    (i,3) = 0.0d0
	ugold (i,3) = 0.0d0
	acg   (i,3) = 0.0d0
	acgold(i,3) = 0.0d0
	  
	IBC   (i,3) = 1
      end do
    endif

    ! Top	 	
    if ((bound(b)%FACE_ID == 6).and.(BCType(6) == 1) ) then  
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)    
        ug    (i,2) = 0.0d0
        ugold (i,2) = 0.0d0
        acg   (i,2) = 0.0d0
        acgold(i,2) = 0.0d0
	  
        IBC   (i,2) = 1
      end do  
    end if
    
  end do

end subroutine setBCs
     



!======================================================================
!
!======================================================================
subroutine setMeshBCs()
  use commonvars
  use aAdjkeep

  implicit none
  
  integer :: b, i, j, k
  real(8) :: rot(NSD)

  do b = 1, NBOUND
 
    if (bound(b)%FACE_ID >= 7) then ! Object
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j) 
  
        rot = matmul(Rn1, (xg(i,:)-xcg))
	dg(i,:) = dbn1 + rot - (xg(i,:)-xcg)
	    
        acgm(i,:) = (1.0d0/beti)*(-(ugmold(i,:)/Delt) &
                  + ((dg(i,:)-dgold(i,:))/(Delt*Delt)) &
               	  + (beti-0.5d0)*acgmold(i,:))
     
        ugm (i,:) = ugmold(i,:) &  
      	          + Delt*((1.0d0-gami)*acgmold(i,:)+ gami*acgm (i,:))
     
        IBC(i,4:6) = 1
      end do 
    end if
	  
    ! Outflow/Inflow
    if ((bound(b)%FACE_ID == 3).or.(bound(b)%FACE_ID == 5)) then
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)   
        dg     (i,1:3) = 0.0d0
        ugm    (i,1:3) = 0.0d0
        acgm   (i,1:3) = 0.0d0
        dgold  (i,1:3) = 0.0d0
        ugmold (i,1:3) = 0.0d0
        acgmold(i,1:3) = 0.0d0     
    
        IBC    (i,4:6) = 1 
      end do 
    end if

    ! Sides		    
    if ((bound(b)%FACE_ID == 2).or.(bound(b)%FACE_ID == 4)) then 
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j) 
        dg     (i,1:3) = 0.0d0
        ugm    (i,1:3) = 0.0d0
        acgm   (i,1:3) = 0.0d0
        dgold  (i,1:3) = 0.0d0
        ugmold (i,1:3) = 0.0d0
        acgmold(i,1:3) = 0.0d0     
    
        IBC    (i,4:6) = 1
      end do  
    end if

    ! Top/bottom	
    if ((bound(b)%FACE_ID == 1).or.(bound(b)%FACE_ID == 6)) then  
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j) 	   
        dg     (i,1:3) = 0.0d0
        ugm    (i,1:3) = 0.0d0
        acgm   (i,1:3) = 0.0d0
        dgold  (i,1:3) = 0.0d0
        ugmold (i,1:3) = 0.0d0
        acgmold(i,1:3) = 0.0d0     
    
        IBC    (i,4:6) = 1
      end do 
    end if
		     
  end do  
   
end subroutine setMeshBCs
