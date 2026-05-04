!========================================================================
! Main routine to call all the subroutines   
! Find closest points between FEM and T-Spline
! f2f means for a Gauss point on f1, we want to find the closest point
! on f2
!========================================================================
subroutine find_close_point_FEM_FEM(FEM1, FEM2, NSD, istep, Rstep, ifq)
  use mpi
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: FEM1, FEM2
  integer,    intent(in)    :: NSD, istep, Rstep, ifq

  integer :: i, j, k, ier, nf1, nf2, nf3, itmp1, itmp2

  integer, allocatable :: f2f_ELM(:,:)
  real(8), allocatable :: f2f_CLP(:,:,:)

  !-------------------------------------------------------
  ! build close element list
  !-------------------------------------------------------
  if (mod(istep,ifq)==0 .or. istep==(Rstep+1)) then
    if (ismaster) then
      write(*,*) "*f2f: Build closest-elements list."
    end if
    FEM1%ELM_Close = -1
    call f2f_find_elm(FEM1, FEM2, nsd, FEM1%NEL_Close, FEM1%ELM_Close)
  end if

  !-------------------------------------------------------
  ! find close points 
  !-------------------------------------------------------
  allocate(f2f_ELM(FEM1%NEL,FEM1%NGAUSS), &
           f2f_CLP(FEM1%NEL,FEM1%NGAUSS,3))
  f2f_ELM = 0; f2f_CLP = 0.0d0

  call f2f_find_point(FEM1, FEM2, nsd, FEM1%NEL_Close, &
                      FEM1%ELM_Close(:,:,1:FEM1%NEL_Close), &
                      f2f_ELM, f2f_CLP, FEM1%NEL_Close_t)

  !-------------------------------------------------------
  ! if the close points locate outside 80% of the element list
  !-------------------------------------------------------
  if (FEM1%NEL_Close_t >= FEM1%NEL_Close*4/5) then
    if (ismaster) then
      write(*,*) "*f2f: Build closest-elements list again!!!"!!!, &
!!!                 FEM1%NEL_Close_t, FEM1%NEL_Close
    end if
    FEM1%ELM_Close = -1
    call f2f_find_elm(FEM1, FEM2, nsd, FEM1%NEL_Close, FEM1%ELM_Close)

    f2f_ELM = 0; f2f_CLP = 0.0d0
    call f2f_find_point(FEM1, FEM2, nsd, FEM1%NEL_Close, &
                        FEM1%ELM_Close(:,:,1:FEM1%NEL_Close), &
                        f2f_ELM, f2f_CLP, FEM1%NEL_Close_t)
  end if

  !-------------------------------------------------------
  ! check if the closest point is between 0 and 1 (for triangles)
  !-------------------------------------------------------
  if (maxval(f2f_CLP(:,:,1:2)) > 1.0d0 .or. &
      minval(f2f_CLP(:,:,1:2)) < 0.0d0) then
    write(*,*) "ERROR: Point located outside of 0 and 1"
    stop
  end if

  FEM1%CLE = f2f_ELM
  FEM1%CLP = f2f_CLP(:,:,1:2)

  deallocate(f2f_ELM, f2f_CLP)

end subroutine find_close_point_FEM_FEM
