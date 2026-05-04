!=============================================================== 
! define the thickness distribution for the root of blade
!===============================================================
! subroutine define_thick_root(zval, thi)

subroutine define_thick_root(zval, thi)


  real(8), intent(in) :: zval
  real(8), intent(out):: thi
!  real(8):: zval, thi

  real(8) :: Zone(14), Thick(14)

  integer :: i, num

  num = 14

!  zval=0.900  
  if (zval>1.0d0 .or. zval<0.0d0) then
    write (*,*) 'Wrong input thickness'
    Stop
  endif

  thi = 0.0d0
  Zone = (/0.000, 0.005, 0.007, 0.009, 0.011, 0.013, 0.024, &
           0.026, 0.047, 0.068, 0.089, 0.114, 0.146, 1.000 /)
  Thick = (/160.0, 140.0, 120.0, 100.0, 80.0,  70.0,  63.0, &
            55.0,  40.0,  25.0,  15.0,   5.0,  0.0,   0.0 /)

  if (abs (zval-1.0d0) < 1.0d-12) then
    thi = Thick (num)
  else  
    do i = 1, num    
      if (Zone(i)>zval) then
        thi= Thick(i)+ (zval-Zone(i))*(Thick(i-1)-Thick(i))/(Zone(i-1)-Zone(i))
        exit
      end if
    end do 
  endif
  thi = thi*1.0d-3/1.64d0
!  write (*,*) 'zval = ', zval, 'Thi =', thi

end subroutine define_thick_root

!=============================================================== 
! define the thickness distribution for the spar cap of blade
!===============================================================
subroutine define_thick_sparcap(zval, thi)

!subroutine define_thick_sparcap

!  implicit none

  real(8), intent(in) :: zval
  real(8), intent(out):: thi
!  real(8):: zval, thi

  real(8) :: Zone(34), Thick(34)

  integer :: i, num

  num = 34

!  zval=0.900  
  if (zval>1.0d0 .or. zval<0.0d0) then
    write (*,*) 'Wrong input thickness'
    Stop
  endif

  thi = 0.0d0
  Zone = (/0.000, 0.005, 0.007, 0.009, 0.011, 0.013, 0.024, &
           0.026, 0.047, 0.068, 0.089, 0.114, 0.146, 0.163, &
           0.179, 0.195, 0.222, 0.249, 0.277, 0.358, 0.439, &
           0.521, 0.602, 0.667, 0.683, 0.732, 0.765, 0.846, &
           0.895, 0.944, 0.957, 0.972, 0.986, 1.000/)
  Thick = (/0.0,    1.0,   2.0,   3.0,   4.0,  10.0,  13.0, &
            13.0,  20.0,  30.0,  51.0,  68.0,  94.0, 111.0, &
           119.0, 136.0, 136.0, 136.0, 128.0, 119.0, 111.0, &
           102.0,  85.0,  68.0,  64.0,  47.0,  34.0,  17.0, &
             9.0,   5.0,   5.0,   5.0,   5.0,   0.0/)

  if (abs (zval-1.0d0) < 1.0d-12) then
    thi = Thick (num)
  else  
    do i = 1, num    
      if (Zone(i)>zval) then
        thi= Thick(i)+ (zval-Zone(i))*(Thick(i-1)-Thick(i))/(Zone(i-1)-Zone(i))
        exit
      end if
    end do 
  endif
  thi = thi*1.0d-3/1.64d0
!  write (*,*) 'zval = ', zval, 'Thi =', thi

end subroutine define_thick_sparcap


!============================================================================= 
! define the thickness distribution for the trailing edge (TE) of blade (ELT)
!=============================================================================
subroutine define_thick_TE_ELT(zval, thi)

!subroutine define_thick_TE_ELT

!  implicit none

  real(8), intent(in) :: zval
  real(8), intent(out):: thi
!  real(8):: zval, thi

  real(8) :: Zone(34), Thick(34)

  integer :: i, num

  num = 34

!  zval=0.900  
  if (zval>1.0d0 .or. zval<0.0d0) then
    write (*,*) 'Wrong input thickness'
    Stop
  endif

  thi = 0.0d0
  Zone = (/0.000, 0.005, 0.007, 0.009, 0.011, 0.013, 0.024, &
           0.026, 0.047, 0.068, 0.089, 0.114, 0.146, 0.163, &
           0.179, 0.195, 0.222, 0.249, 0.277, 0.358, 0.439, &
           0.521, 0.602, 0.667, 0.683, 0.732, 0.765, 0.846, &
           0.895, 0.944, 0.957, 0.972, 0.986, 1.000/)
  Thick = (/0.0,    1.0,   2.0,   3.0,   5.0,   7.0,   8.0, &
            9.0,   13.0,  18.0,  25.0,  33.0,  40.0,  50.0, &
           60.0,   60.0,  60.0,  60.0,  30.0,  30.0,  15.0, &
            8.0,    4.0,   4.0,   4.0,   4.0,   4.0,   4.0, &
            4.0,    4.0,   4.0,   4.0,   4.0,   0.0/)

  if (abs (zval-1.0d0) < 1.0d-12) then
    thi = Thick (num)
  else  
    do i = 1, num    
      if (Zone(i)>zval) then
        thi= Thick(i)+ (zval-Zone(i))*(Thick(i-1)-Thick(i))/(Zone(i-1)-Zone(i))
        exit
      end if
    end do 
  endif
  thi = thi*1.0d-3/1.64d0
!  write (*,*) 'zval = ', zval, 'Thi =', thi

end subroutine define_thick_TE_ELT

!============================================================================= 
! define the thickness distribution for the trailing edge (TE) of blade (foam)
!=============================================================================
subroutine define_thick_TE_foam(zval, thi)

!subroutine define_thick_TE_foam

!  implicit none

  real(8), intent(in) :: zval
  real(8), intent(out):: thi
!  real(8):: zval, thi

  real(8) :: Zone(34), Thick(34)

  integer :: i, num

  num = 34

!  zval=0.900  
  if (zval>1.0d0 .or. zval<0.0d0) then
    write (*,*) 'Wrong input thickness'
    Stop
  endif

  thi = 0.0d0
  Zone = (/0.000, 0.005, 0.007, 0.009, 0.011, 0.013, 0.024, &
           0.026, 0.047, 0.068, 0.089, 0.114, 0.146, 0.163, &
           0.179, 0.195, 0.222, 0.249, 0.277, 0.358, 0.439, &
           0.521, 0.602, 0.667, 0.683, 0.732, 0.765, 0.846, &
           0.895, 0.944, 0.957, 0.972, 0.986, 1.000/)
  Thick = (/0.0,    0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
            0.0,    0.0,   0.0,  60.0,  60.0,  60.0,  60.0, &
           60.0,   60.0,  60.0,  60.0,  40.0,  40.0,  20.0, &
           10.0,   10.0,  10.0,  10.0,  10.0,  10.0,  10.0, &
           10.0,   10.0,  10.0,  10.0,  10.0,   0.0/)

  if (abs (zval-1.0d0) < 1.0d-12) then
    thi = Thick (num)
  else  
    do i = 1, num    
      if (Zone(i)>zval) then
        thi= Thick(i)+ (zval-Zone(i))*(Thick(i-1)-Thick(i))/(Zone(i-1)-Zone(i))
        exit
      end if
    end do 
  endif
  thi = thi*1.0d-3/1.64d0
!  write (*,*) 'zval = ', zval, 'Thi =', thi

end subroutine define_thick_TE_foam

!============================================================================= 
! define the thickness distribution for the leading edge panel of blade
!=============================================================================
subroutine define_thick_LE(zval, thi)

!subroutine define_thick_LE

!  implicit none

  real(8), intent(in) :: zval
  real(8), intent(out):: thi
!  real(8):: zval, thi

  real(8) :: Zone(34), Thick(34)

  integer :: i, num

  num = 34

!  zval=0.900  
  if (zval>1.0d0 .or. zval<0.0d0) then
    write (*,*) 'Wrong input thickness'
    Stop
  endif

  thi = 0.0d0
  Zone = (/0.000, 0.005, 0.007, 0.009, 0.011, 0.013, 0.024, &
           0.026, 0.047, 0.068, 0.089, 0.114, 0.146, 0.163, &
           0.179, 0.195, 0.222, 0.249, 0.277, 0.358, 0.439, &
           0.521, 0.602, 0.667, 0.683, 0.732, 0.765, 0.846, &
           0.895, 0.944, 0.957, 0.972, 0.986, 1.000/)
  Thick = (/0.0,    0.0,   0.0,   0.0,   0.0,   1.0,   3.5, &
           13.0,  100.0, 100.0, 100.0, 100.0, 100.0,  60.0, &
           60.0,   60.0,  60.0,  60.0,  60.0,  60.0,  60.0, &
           60.0,   60.0,  60.0,  55.0,  45.0,  30.0,  15.0, &
           10.0,    0.0,   0.0,   0.0,   0.0,   0.0/)

  if (abs (zval-1.0d0) < 1.0d-12) then
    thi = Thick (num)
  else  
    do i = 1, num    
      if (Zone(i)>zval) then
        thi= Thick(i)+ (zval-Zone(i))*(Thick(i-1)-Thick(i))/(Zone(i-1)-Zone(i))
        exit
      end if
    end do 
  endif
  thi = thi*1.0d-3/1.64d0
!  write (*,*) 'zval = ', zval, 'Thi =', thi

end subroutine define_thick_LE

!============================================================================= 
! define the thickness distribution for the aft panel of blade
!=============================================================================
subroutine define_thick_Aft(zval, thi)

!subroutine define_thick_Aft

!  implicit none

  real(8), intent(in) :: zval
  real(8), intent(out):: thi
!  real(8):: zval, thi

  real(8) :: Zone(34), Thick(34)

  integer :: i, num

  num = 34

!  zval=0.900  
  if (zval>1.0d0 .or. zval<0.0d0) then
    write (*,*) 'Wrong input thickness'
    Stop
  endif

  thi = 0.0d0
  Zone = (/0.000, 0.005, 0.007, 0.009, 0.011, 0.013, 0.024, &
           0.026, 0.047, 0.068, 0.089, 0.114, 0.146, 0.163, &
           0.179, 0.195, 0.222, 0.249, 0.277, 0.358, 0.439, &
           0.521, 0.602, 0.667, 0.683, 0.732, 0.765, 0.846, &
           0.895, 0.944, 0.957, 0.972, 0.986, 1.000/)
  Thick = (/0.0,    0.0,   0.0,   0.0,   0.0,   1.0,   3.5, &
           13.0,   30.0,  50.0,  60.0,  60.0,  60.0,  60.0, &
           60.0,   60.0,  60.0,  60.0,  60.0,  60.0,  60.0, &
           60.0,   60.0,  60.0,  55.0,  45.0,  30.0,  15.0, &
           10.0,    5.0,   5.0,   5.0,   5.0,   0.0/)

  if (abs (zval-1.0d0) < 1.0d-12) then
    thi = Thick (num)
  else  
    do i = 1, num    
      if (Zone(i)>zval) then
        thi= Thick(i)+ (zval-Zone(i))*(Thick(i-1)-Thick(i))/(Zone(i-1)-Zone(i))
        exit
      end if
    end do 
  endif
  thi = thi*1.0d-3/1.64d0
!  write (*,*) 'zval = ', zval, 'Thi =', thi

end subroutine define_thick_Aft


!-----------------------------------------------------------------------------------
!---------------------------Ming-Chen Thickness-------------------------------------
!-----------------------------------------------------------------------------------
!=============================================================== 
! define the thickness distribution for the blade
!===============================================================
subroutine define_thick_sandia(yval, thi)

  implicit none

  real(8), intent(in) :: yval
  real(8), intent(out):: thi

  real(8) :: TZone(3), Thick(3), thick_lin

  TZone(1) = 0.7d0
  Thick(1) = 0.0114285714d0

  TZone(2) = 1.8d0
  Thick(2) = 0.00857142857d0

  Tzone(3) = 9.0d0
  Thick(3) = 0.00285714286d0

  if (yval <= TZone(1)) then
    thi = Thick(1)

  else if (yval <= Tzone(2)) then
    thi = thick_lin(yval,Thick(1),Thick(2),Tzone(1),Tzone(2))

  else if (yval <= Tzone(3)) then
    thi = thick_lin(yval,Thick(2),Thick(3),Tzone(2),Tzone(3))

  else
!   write(*,*) "ERROR: Thickness profile undefined!!!"
!   stop
    thi = Thick(3)
  end if
end subroutine define_thick_sandia




!=============================================================== 
! define the thickness distribution for the blade
!===============================================================
subroutine define_thick(yval, thi)

  implicit none

  real(8), intent(in) :: yval
  real(8), intent(out):: thi

  real(8) :: TZone(7), Thick(7), thick_lin

  TZone(1) = 5.5d0
  Thick(1) = 0.08d0

  TZone(2) = 10.0d0
  Thick(2) = 0.065d0

  Tzone(3) = 15.8d0
  Thick(3) = 0.05d0
      
  Tzone(4) = 54.5d0
  Thick(4) = 0.029d0

  Tzone(5) = 59.9d0
  Thick(5) = 0.025d0

  Tzone(6) = 62.1d0
  Thick(6) = 0.0225d0

  Tzone(7) = 63.0d0
  Thick(7) = 0.02d0

  if (yval <= TZone(1)) then
    thi = Thick(1)

  else if (yval <= Tzone(2)) then
    thi = thick_lin(yval,Thick(1),Thick(2),Tzone(1),Tzone(2))

  else if (yval <= Tzone(3)) then
    thi = thick_lin(yval,Thick(2),Thick(3),Tzone(2),Tzone(3))

  else if (yval <= Tzone(4)) then
    thi = thick_lin(yval,Thick(3),Thick(4),Tzone(3),Tzone(4))

  else if (yval <= Tzone(5)) then
    thi = thick_lin(yval,Thick(4),Thick(5),Tzone(4),Tzone(5))

  else if (yval <= Tzone(6)) then
    thi = thick_lin(yval,Thick(5),Thick(6),Tzone(5),Tzone(6))

  else if (yval <= Tzone(7)) then
    thi = thick_lin(yval,Thick(6),Thick(7),Tzone(6),Tzone(7))

  else
!   write(*,*) "ERROR: Thickness profile undefined!!!"
!   stop
    thi = Thick(7)
  end if

!$$$      !----------------------------------------------------
!$$$      ! evaluate the thickness, linear distribution so far
!$$$      ! root:0.1m -> tip:0.02m (was 0.3-0.04)
!$$$      ! vknot:0.0 -> 1.8 or 1.766666667
!$$$      !----------------------------------------------------
!$$$      TZone(1) = 0.15d0   ! Root zone, 8cm
!$$$      Thick(1) = 0.08d0
!$$$
!$$$      TZone(2) = 0.45d0   ! Trasition zone, 8cm->2cm
!$$$!      Thick(2) = 0.02d0
!$$$      Thick(2) = 0.05d0
!$$$
!$$$      Tzone(3) = 1.8d0     ! Aerodynamic zone, 2cm->0.5cm
!$$$!      Thick(3) = 0.005d0
!$$$      Thick(3) = 0.02d0
      
!$$$      !----------------------------------------------------
!$$$      ! evaluate the thickness, linear distribution so far
!$$$      !----------------------------------------------------
!$$$
!$$$      if (VVal <= TZone(1)) then
!$$$        thi = Thick(1) 
!$$$      else if (VVal <= Tzone(2)) then
!$$$        thi = (Thick(2)-Thick(1))*(VVal-Tzone(1))/
!$$$     &        (TZone(2)-Tzone(1)) + Thick(1)
!$$$      else if (VVal <= Tzone(3)) then
!$$$        thi = (Thick(3)-Thick(2))*(VVal-Tzone(2))/
!$$$     &        (TZone(3)-Tzone(2)) + Thick(2)
!$$$      else
!$$$        write(*,*) "ERROR: Thickness profile undefined!!!"
!$$$        stop
!$$$      end if

end subroutine define_thick

!-------------------------------------
! linear function for thickness
!-------------------------------------
function thick_lin(yval, thick1, thick2, tzone1, tzone2)
  implicit none
  real(8) :: thick_lin
  real(8), intent(in) :: yval, thick1, thick2, tzone1, tzone2

  thick_lin = (thick2-thick1)*(yval-tzone1)/(tzone2-tzone1)+thick1

end function thick_lin


!===============================================================
! define the thickness distribution for the blade
!===============================================================
subroutine define_thick_breitenberger(xi, thi)
  implicit none
  real(8), intent(in) :: xi(3)
  real(8), intent(out):: thi
  integer :: k
  real(8) :: thifac, tmpfac
  real(8) :: root, rout(0:7), rin(0:7)
  real(8) :: thiout(0:7), thiin(0:7)
  real(8) :: TZone(0:7), Thick(0:7), thick_lin
  real(8) :: get_mid, get_rad

  thifac = 0.342482d0   ! same volume as the old design

  ! root
  TZone(0) = 5.5d0
  Thick(0) = 0.06d0
  root     = 4.9d0

  ! < 5.5
  TZone(1)  = 5.5d0
  Thick(1)  = 0.06d0
  rout(1)   = 0.85d0
  rin(1)    = 0.8d0
  thiin(1)  = 4.9d0
  thiout(1) = 1.3d0

  ! < 10.0
  TZone(2)  = 10.0d0
  Thick(2)  = 0.06d0
  rout(2)   = 0.85d0
  rin(2)    = 0.8d0
  thiin(2)  = 4.9d0
  thiout(2) = 1.3d0

  ! < 15.8
  Tzone(3)  = 15.8d0
  Thick(3)  = 0.05d0
  rout(3)   = 0.85d0
  rin(3)    = 0.8d0
  thiin(3)  = 4.9d0
  thiout(3) = 1.3d0

  ! < 54.8
  Tzone(4)  = 54.5d0
  Thick(4)  = 0.029d0
  rout(4)   = 0.5d0
  rin(4)    = 0.4d0
  thiin(4)  = 4.9d0
  thiout(4) = 0.9d0

  ! < 59.9
  Tzone(5)  = 59.9d0
  Thick(5)  = 0.025d0
  rout(5)   = 0.4d0
  rin(5)    = 0.3d0
  thiin(5)  = 4.9d0
  thiout(5) = 0.7d0

  ! < 62.1
  Tzone(6)  = 62.1d0
  Thick(6)  = 0.0225d0
  rout(6)   = 0.3d0
  rin(6)    = 0.2d0
  thiin(6)  = 4.9d0
  thiout(6) = 0.7d0

  ! < 63.0
  Tzone(7)  = 63.0d0
  Thick(7)  = 0.02d0
  rout(7)   = 0.4d0
  rin(7)    = 0.1d0
  thiin(7)  = 0.7d0
  thiout(7) = 0.7d0

  if (xi(2) <= Tzone(0)) then
    thi = (Thick(1)*thifac*root)

  else if (xi(2) <= TZone(1)) then
    k = 1
    ! send in thick(k-1) = thick(k) so that the thickness 
    ! function will be "thi = Thick(1)*thifac" only.
    ! because thick(2)-thick(1) = 0 in this case
    call get_thick(xi, (/thick(k-1),thick(k-1)/), tzone(k-1:k), &
                   rout(k-1:k), rin(k-1:k), thiin(k),  &
                   thiout(k), thifac, thi)

  else if (xi(2) <= TZone(2)) then        
    k = 2
    call get_thick(xi, thick(k-1:k), tzone(k-1:k), rout(k-1:k), &
                   rin(k-1:k), thiin(k), thiout(k), thifac, thi)

  else if (xi(2) <= TZone(3)) then        
    k = 3
    call get_thick(xi, thick(k-1:k), tzone(k-1:k), rout(k-1:k), &
                   rin(k-1:k), thiin(k), thiout(k), thifac, thi)

  else if (xi(2) <= TZone(4)) then        
    k = 4
    call get_thick(xi, thick(k-1:k), tzone(k-1:k), rout(k-1:k), &
                   rin(k-1:k), thiin(k), thiout(k), thifac, thi)

  else if (xi(2) <= TZone(5)) then        
    k = 5
    call get_thick(xi, thick(k-1:k), tzone(k-1:k), rout(k-1:k), &
                   rin(k-1:k), thiin(k), thiout(k), thifac, thi)

  else if (xi(2) <= TZone(6)) then
    k = 6
    call get_thick(xi, thick(k-1:k), tzone(k-1:k), rout(k-1:k), &
                   rin(k-1:k), thiin(k), thiout(k), thifac, thi)

  else
    k = 7
    call get_thick(xi, thick(k-1:k), tzone(k-1:k), rout(k-1:k), &
                   rin(k-1:k), thiin(k), thiout(k), thifac, thi)
  end if

end subroutine define_thick_breitenberger


!-------------------------------------
! get mid point for spar cap
!-------------------------------------
function get_mid(xi)
  implicit none
   
  real(8), intent(in) :: xi(3)
  real(8):: root, tip, get_mid
  root = 0.15d0
  tip  = 0.05d0

  get_mid = (tip-root)*(xi(2)-2.0d0)/(63.0d0-2.0d0)+root

end function get_mid

!-------------------------------------
! get radius for spar cap
!-------------------------------------
function get_rad(xi, rad, tzone)
  implicit none
   
  real(8), intent(in) :: xi(3), rad(2), tzone(2)
  real(8):: root, tip, get_rad

  get_rad = rad(1)-(rad(1)-rad(2))*(xi(2)-tzone(1))/(tzone(2)-tzone(1))

end function get_rad

!-------------------------------------
! get the thickness
!-------------------------------------
subroutine get_thick(xi, thick, tzone, radout, radin, &
                     thiin, thiout, thifac, thi)

  implicit none

  real(8), intent(in) :: xi(3), thick(2), tzone(2), radout(2), &
                         radin(2), thiin, thiout, thifac
  real(8), intent(out):: thi
  real(8) :: get_mid, get_rad, dist(2), tmpfac, rout, rin, mid

  thi = Thick(1)*thifac-(xi(2)-TZone(1))/(TZone(2)-TZone(1))* &
        (Thick(1)-Thick(2))*thifac

  ! Spar cab
  rout = get_rad(xi,radout,tzone)
  rin  = get_rad(xi,radin, tzone)

  mid = get_mid(xi)

  dist(1) = mid-rout
  dist(2) = mid+rout
  if (xi(1)>dist(1) .and. xi(1)<dist(2)) then
    dist(1) = mid-rin
    dist(2) = mid+rin
    if (xi(1)>dist(1) .and. xi(1)<dist(2)) then
      thi = thiin*thi
    else
      tmpfac = thiin-(abs(xi(1)-mid)-rin)/(rout-rin)*(thiin-thiout)
      thi = thi*tmpfac
    end if
  else
    thi =  thi*thiout
  end if
end subroutine get_thick
