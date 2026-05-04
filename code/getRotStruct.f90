!======================================================================
! preserves the propper rotation
!======================================================================
  subroutine integrateBodyRotation(Qb1,Qb0,Rfix,wb,Dtgl)

    implicit none
  
    real(8), intent(out) :: Qb1(3,3)
    real(8), intent(in)  :: Qb0(3,3),wb(3),Dtgl, Rfix(3,3)
  
    real(8) :: RHS(3,3),LHS(3,3),Ainv(3,3),det, Rtemp(3,3)
  
    integer :: i,j,k
    real(8) :: Eijk

    RHS = Qb0*Dtgl
    LHS = 0d0   
    do i = 1, 3  
      LHS(i,i) = Dtgl
      do j = 1, 3   
        do k = 1, 3
          RHS(i,:) = RHS(i,:) + 0.5d0*Eijk(i,j,k)*wb(j)*Qb0(k,:)
          LHS(i,k) = LHS(i,k) - 0.5d0*Eijk(i,j,k)*wb(j)
        enddo
      enddo     
    enddo  
          
    call get_inverse_3x3(LHS,Ainv,det)
  
    do i = 1, 3
      do j = 1, 3
        Qb1(i,j) = sum(Ainv(i,:)*RHS(:,j))
      enddo
    enddo
        do i = 1, 3
      do j = 1, 3
        Rtemp(i,j) = sum(Rfix(i,:)*Qb1(:,j))
      end do
    end do

    Qb1 = Rtemp
  
  end subroutine integrateBodyRotation

!---------------------------------------------------------------------- 
! subroutine to get the omega for rotation
!---------------------------------------------------------------------- 
subroutine get_omega(NRB, Density, omega_full,omega_tang)
  use commonvars
  use aAdjkeep
  use defs_shell
  implicit none

  type(mesh), intent(in) :: NRB 
  real(8), intent(in) :: Density(NRB%NEL,NRB%NGAUSS**2)
!  real(8), intent(out) :: Omega_tensor(3,3), Omega_tang_tensor(3,3)
  real(8), intent(out) :: omega_full(3), omega_tang(3)
  real(8) :: M(3), Inert(3,3), Iinv(3,3), det
  integer :: b, i, j, k, ii, jj, iel, igauss, jgauss, ct

  integer :: p, q, nshl, nuk, nvk, ptype, &
             kk, ni, nj, aa, bb

  real(8) :: gp(NRB%NGAUSS), gw(NRB%NGAUSS), gwt, da, VVal, &
             DetJb_SH, nor(NSD), dens, &
             xu(NSD), xd(NSD), vel(NSD), dxdxi(NSD,2), ddxddxi(nsd,3)

  integer, allocatable :: lIEN(:)
  real(8), allocatable :: shl(:), shgradl(:,:), shhessl(:,:)

  real(8), allocatable :: mprod(:), iprod(:,:)	

  do i = 1, NRB%NNODE
    if((abs(NRB%B_NET(i,3)- 3.0d0) <=1.0d-4).and.(abs(NRB%B_NET(i,1)- &
      0.0d0) <=1.0d-4).and.(abs(NRB%B_NET(i,2)-0.0d0) <=1.0d-4)) then
      xrot(:) = NRB%B_NET_D(i,1:3) !Center of rotation
    elseif((abs(NRB%B_NET(i,3)- 9.0d0) <=1.0d-4).and.(abs(NRB%B_NET(i,1)- &
      0.0d0) <=1.0d-4).and.(abs(NRB%B_NET(i,2)-0.0d0) <=1.0d-4)) then
      xh_E(:) = NRB%B_NET_D(i,1:3) !End of the hub
    end if
  end do


  normal = (xh_E - xrot)
  normal = normal/sqrt(sum(normal(:)**2)) !Unit normal

  write(*,*) 'Hub:',xrot
  write(*,*) xh_E


  gp = 0.0d0; gw = 0.0d0  
  DetJb_SH = 0.0d0

  ! get Gaussian points and weights     
  call genGPandGW_shell(gp, gw, NRB%NGAUSS) 
  allocate(mprod(3), iprod(3,3))
  mprod = 0.0d0 ; iprod = 0.0d0 !Local Angular momentum and angular moment of inetria
  ! loop over elements
  do iel = 1, NRB%NEL
  
    ! get NURB coordinates
    ni = NRB%INN(iel,1); nj = NRB%INN(iel,2)
    
    ! Check to see if current element has nonzero area, 
    ! skip if it doesn't
    if ((NRB%U_KNOT(iel,ni) /= NRB%U_KNOT(iel,ni+1)) .and. &
        (NRB%V_KNOT(iel,nj) /= NRB%V_KNOT(iel,nj+1))   ) then
         
      ! used in calculating quadrature points. The factor of 4.0d0
      ! comes from mapping from the [-1,1] line onto a real segment...
      da = (NRB%U_KNOT(iel,ni+1)-NRB%U_KNOT(iel,ni))*  &
           (NRB%V_KNOT(iel,nj+1)-NRB%V_KNOT(iel,nj))/4.0d0

      p = NRB%P(iel); nuk = NRB%NUK(iel); nshl = NRB%NSHL(iel)
      q = NRB%Q(iel); nvk = NRB%NVK(iel); ptype = NRB%PTYPE(iel)

    allocate(shl(nshl), shgradl(nshl,2), shhessl(nshl,3), lIEN(nshl))




    lIEN = -1
    do i = 1, nshl
      lIEN(i) = NRB%IEN(iel,i)
    end do

    ! Loop over integration points (NGAUSS in each direction)
    ct = 0 
    do jgauss = 1, NRB%NGAUSS
      do igauss = 1, NRB%NGAUSS

        ct = ct + 1
          
          ! Get Element Shape functions and their gradients
          shl = 0.0d0; shgradl = 0.0d0; shhessl = 0.0d0
          xu = 0.0d0; xd = 0.0d0; dxdxi = 0.0d0; ddxddxi = 0.0d0
          nor = 0.0d0; vel = 0.0d0

          call eval_SHAPE_shell(gp(igauss), gp(jgauss),  &
                                shl, shgradl, shhessl, nor,  &
                                xu, xd, dxdxi, ddxddxi,  &
                                p, q, nsd, nshl, &
                                lIEN, NRB%NNODE, &
                                NRB%B_NET_U, NRB%B_NET_D, DetJb_SH, &
                                ni, nj, nuk, nvk, &
                                NRB%U_KNOT(iel,1:nuk), &
                                NRB%V_KNOT(iel,1:nvk))

          do i = 1, nshl
            do ii = 1, 3
              vel(ii) = vel(ii) + NRB%ush(lIEN(i),ii)*shl(i)
            end do
          end do


          gwt = gw(igauss)*gw(jgauss)*da

          dens = Density(iel, ct)

          if ( ptype .eq. 7 ) then !Only integrate over the hub

              mprod(1) = mprod(1) + ((xd(2)-xrot(2))* &
                              dens*vel(3)- &
                              (xd(3)-xrot(3))* &
                              dens*vel(2))*DetJb_SH*gwt
              mprod(2) = mprod(2) - ((xd(1)-xrot(1))* &
                              dens*vel(3)- &
                              (xd(3)-xrot(3))* &
                              dens*vel(1))*DetJb_SH*gwt
	      mprod(3) = mprod(3) + ((xd(1)-xrot(1))* &
                              dens*vel(2)- &
                              (xd(2)-xrot(2))* &
                              dens*vel(1))*DetJb_SH*gwt
		
              do ii = 1, 3
	        do jj = 1, 3
		  iprod(ii,jj) = iprod(ii,jj) + (dens*sum((xd(:)-xrot(:))* &
                                    (xd(:)-xrot(:)))*Identity(ii,jj)- &
		                     dens*(xd(ii)-xrot(ii))* &
                                     (xd(jj)-xrot(jj)))*DetJb_SH*gwt
		end do
              end do      


          end if
        end do
      end do  ! end loop gauss points


    deallocate(shl, shgradl, shhessl, lIEN)

    end if  ! end if nonzero areas elements

  end do    ! end loop elements

  call get_inverse_3x3(iprod,Iinv,det)

  do i = 1, 3
    omega_full(i) = sum(Iinv(i,:)*mprod(:))
  end do

  omega_tang(:) = omega_full(:) - sum(omega_full(:)*normal(:))*normal(:) !Tangential component

!Construct angular tensors

!  Omega_tensor = 0.0d0
!  Omega_tensor(1,2) = -omega_full(3)
!  Omega_tensor(1,3) =  omega_full(2)
!  Omega_tensor(2,1) =  omega_full(3)
!  Omega_tensor(2,3) = -omega_full(1)
! Omega_tensor(3,1) = -omega_full(2)
!  Omega_tensor(3,2) =  omega_full(1)

!  Omega_tang_tensor = 0.0d0
!  Omega_tang_tensor(1,2) = -omega_tang(3)
!  Omega_tang_tensor(1,3) =  omega_tang(2)
!  Omega_tang_tensor(2,1) =  omega_tang(3)
!  Omega_tang_tensor(2,3) = -omega_tang(1)
!  Omega_tang_tensor(3,1) = -omega_tang(2)
!  Omega_tang_tensor(3,2) =  omega_tang(1)

  deallocate(mprod, iprod)

end subroutine get_omega

!----------------------------------------------------------------------
!  Function for levi-civita symbol ==> for outerproduct
!---------------------------------------------------------------------- 
function Eijk(i,j,k)
  implicit none
  real(8)  Eijk
  integer :: i,j,k
  Eijk = real((j-i)*(k-i)*(k-j)/2,8)
end function Eijk

