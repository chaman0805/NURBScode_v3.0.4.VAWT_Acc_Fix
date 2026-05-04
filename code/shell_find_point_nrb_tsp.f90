!========================================================================
! Main routine to call all the subroutines        
! Find closest points between NURBS and T-Spline
!========================================================================
subroutine find_close_point_NRB_TSP(NRB, TSP, BEZ, NSD)

  use mpi
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: NRB, TSP, BEZ
  integer,    intent(in) :: NSD

  integer :: i, j, k, ier, nf1, nf2, nf3, itmp1, itmp2

  integer, allocatable :: CLOSE_ELM(:,:,:)
  integer :: NEL_CLOSE
    
  integer, allocatable :: n2t_ELM(:,:),   t2n_ELM(:,:)
  real(8), allocatable :: n2t_CLP(:,:,:), t2n_CLP(:,:,:)

  real(8) :: clock1, clock2

  character(len=30) :: fname, ch1, ch2, fmat

  nf1 = 11; nf2 = 21; nf3 = 31

  !************************************************************
  ! Begin finding the closest points on NURBS surface
  !************************************************************
  ! first check if the file exist. If exist, then read the 
  ! closest point from the files instead of building them
  open(nf1, file='n2t_close_point.dat', status='old', iostat=ier)
  if (ier == 0) then

    if (ismaster) then
      write(*,*) "n2t: Reading the closest-points list"
    end if

    read(nf1,*) itmp1, itmp2
    if (itmp1/=TSP%NEL .or. itmp2/=TSP%NGAUSS**2) then
      write(*,*) "ERROR: TSP%NEL or TSP%NGAUSS does not match the"
      write(*,*) "       numbers in n2t_close_point.dat"
      stop
    end if
    
    do i = 1, TSP%NEL
      do j = 1, TSP%NGAUSS**2
        read(nf1,*) TSP%CLE(i,j), TSP%CLP(i,j,:)
      end do
    end do

    ! check if the closest point is between -1 and 1
    if (maxval(abs(TSP%CLP)) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

  ! if the file does not exist, then build it
  else

    if (ismaster) then
      write(*,*) "n2t: Closest-points list does not exist"
      write(*,*) "     Build it:"
    end if

    !------------------------------------------------------------
    ! first read in (if exist) or build the list of closest 
    ! elements (in ascending order of distance)
    !------------------------------------------------------------
    open(nf2, file='n2t_close_elm.dat', status='old', iostat=ier)
    if (ier == 0) then

      if (ismaster) then
        write(*,*) "n2t: Reading the closest-elements list"
      end if
      read(nf2,*) NEL_CLOSE
      allocate(CLOSE_ELM(TSP%NEL,TSP%NGAUSS**2,NEL_CLOSE))
      CLOSE_ELM = 0
      do i = 1, TSP%NEL
        do j = 1, TSP%NGAUSS**2
          read(nf2,*) CLOSE_ELM(i,j,1:NEL_CLOSE)
        end do
      end do

    else

      if (ismaster) then
        write(*,*) "n2t: Closest-elements list does not exist"
        write(*,*) "     Build it:"
      end if

      ! Save the first 5% closest elements into the file
      NEL_CLOSE = int(NRB%NEL*0.05d0)
      if (ismaster) then
        write(*,*) "n2t: Number of closest elements:", NEL_CLOSE
      end if
      allocate(CLOSE_ELM(TSP%NEL,TSP%NGAUSS**2,NEL_CLOSE))
      CLOSE_ELM = 0

      ! we first check and save the 50% closest elements
      ! check the distance between each element's center (gp=0)
      call cpu_time(clock1)
      call n2t_find_elm(TSP, BEZ, NRB, nsd, NEL_CLOSE, CLOSE_ELM)
      call cpu_time(clock2)
      if (ismaster) then
        write(*,*) "time (sec):", clock2-clock1
      end if
      open(nf3, file='n2t_close_elm.dat', status='new')
      write(nf3,*) NEL_CLOSE, TSP%NEL

      write(ch1,'(I10)') NRB%NEL
      write(ch2,'(I10)') len(trim(adjustl(ch1)))+1
      write(ch1,'(I10)') TSP%NEL
      fmat = '('//trim(adjustl(ch1))//'I'//trim(adjustl(ch2))//')'

      do i = 1, TSP%NEL
        do j = 1, TSP%NGAUSS**2
          write(nf3,fmat) CLOSE_ELM(i,j,:)
        end do
      end do

      close(nf3)
    end if
    close(nf2)

    !------------------------------------------------------------
    ! find the closest point based on n% of the closest elements
    ! 2% gave me exactly the same result as without sorting
    ! 1% is not enough...
    !------------------------------------------------------------
    NEL_CLOSE = int(NRB%NEL*0.005d0)
    if (ismaster) then
      write(*,*) "n2t: Number of closest elements:", NEL_CLOSE
    end if
    allocate(n2t_ELM(TSP%NEL,TSP%NGAUSS**2), &
             n2t_CLP(TSP%NEL,TSP%NGAUSS**2,3))
    n2t_ELM = 0; n2t_CLP = 0.0d0

    call cpu_time(clock1)
    call n2t_find_point(TSP, BEZ, NRB, nsd, NEL_CLOSE, &
                           CLOSE_ELM(:,:,1:NEL_CLOSE), &
                           n2t_ELM, n2t_CLP)
    call cpu_time(clock2)
    if (ismaster) then
      write(*,*) "time (sec):", clock2-clock1
    end if

    ! check if the closest point is between -1 and 1
    if (maxval(abs(n2t_CLP(:,:,1:2))) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

    ! output
    open(nf2, file='n2t_close_point.dat', status='new')
    write(nf2,*) TSP%NEL, TSP%NGAUSS**2
    do i = 1, TSP%NEL
      do j = 1, TSP%NGAUSS**2
        write(nf2,'(I8,3F20.15)') n2t_ELM(i,j), n2t_CLP(i,j,:)
      end do
    end do
    close(nf2)

    TSP%CLE = n2t_ELM
    TSP%CLP = n2t_CLP(:,:,1:2)

    deallocate(CLOSE_ELM, n2t_ELM, n2t_CLP)
  end if
  close(nf1)

  !************************************************************
  ! Begin finding the closest points on T-Spline Surface
  !************************************************************
  ! first check if the file exist. If exist, then read the 
  ! closest point from the files instead of building them
  open(nf1, file='t2n_close_point.dat', status='old', iostat=ier)
  if (ier == 0) then

    if (ismaster) then
      write(*,*) "t2n: Reading the closest-points list"
    end if

    read(nf1,*) itmp1, itmp2
    if (itmp1/=NRB%NEL .or. itmp2/=NRB%NGAUSS**2) then
      write(*,*) "ERROR: NRB%NEL or NRB%NGAUSS does not match the"
      write(*,*) "       numbers in n2t_close_point.dat"
      stop
    end if
    
    do i = 1, NRB%NEL
      do j = 1, NRB%NGAUSS**2
        read(nf1,*) NRB%CLE(i,j), NRB%CLP(i,j,:)
      end do
    end do

    ! check if the closest point is between -1 and 1
    if (maxval(abs(NRB%CLP)) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

  ! if the file does not exist, then build it
  else

    if (ismaster) then
      write(*,*) "t2n: Closest-points list does not exist"
      write(*,*) "     Build it:"
    end if

    !------------------------------------------------------------
    ! first read in (if exist) or build the list of closest 
    ! elements (in ascending order of distance)
    !------------------------------------------------------------
    open(nf2, file='t2n_close_elm.dat', status='old', iostat=ier)
    if (ier == 0) then

      if (ismaster) then
        write(*,*) "t2n: Reading the closest-elements list"
      end if
      read(nf2,*) NEL_CLOSE
      allocate(CLOSE_ELM(NRB%NEL,NRB%NGAUSS**2,NEL_CLOSE))
      CLOSE_ELM = 0
      do i = 1, NRB%NEL
        do j = 1, NRB%NGAUSS**2
          read(nf2,*) CLOSE_ELM(i,j,1:NEL_CLOSE)
        end do
      end do

    else

      if (ismaster) then
        write(*,*) "t2n: Closest-elements list does not exist"
        write(*,*) "     Build it:"
      end if

      ! Save the first 5% closest elements into the file
      NEL_CLOSE = int(TSP%NEL*0.05d0)
      if (ismaster) then
        write(*,*) "t2n: Number of closest elements:", NEL_CLOSE
      end if
      allocate(CLOSE_ELM(NRB%NEL,NRB%NGAUSS**2,NEL_CLOSE))
      CLOSE_ELM = 0

      ! we first check and save the 50% closest elements
      ! check the distance between each element's center (gp=0)
      call cpu_time(clock1)
      call t2n_find_elm(NRB, TSP, BEZ, nsd, NEL_CLOSE, CLOSE_ELM) 
      call cpu_time(clock2)
      if (ismaster) then
        write(*,*) "time (sec):", clock2-clock1, (clock2-clock1)/60.0d0
      end if

      open(nf3, file='t2n_close_elm.dat', status='new')
      write(nf3,*) NEL_CLOSE, NRB%NEL

      write(ch1,'(I10)') TSP%NEL
      write(ch2,'(I10)') len(trim(adjustl(ch1)))+1
      write(ch1,'(I10)') NRB%NEL
      fmat = '('//trim(adjustl(ch1))//'I'//trim(adjustl(ch2))//')'

      do i = 1, NRB%NEL
        do j = 1, NRB%NGAUSS**2
          write(nf3,fmat) CLOSE_ELM(i,j,:)
        end do
      end do
      close(nf3)
    end if
    close(nf2)

    !------------------------------------------------------------
    ! find the closest point based on n% of the closest elements
    ! 2% gave me exactly the same result as without sorting
    ! 1% is not enough...
    !------------------------------------------------------------
    NEL_CLOSE = int(TSP%NEL*0.008d0)   ! 0.035
    if (ismaster) then
      write(*,*) "t2n: Number of closest elements:", NEL_CLOSE
    end if

    allocate(t2n_ELM(NRB%NEL,NRB%NGAUSS**2), &
             t2n_CLP(NRB%NEL,NRB%NGAUSS**2,3))
    t2n_ELM = 0; t2n_CLP = 0.0d0

    call cpu_time(clock1)
    call t2n_find_point(NRB, TSP, BEZ, nsd, NEL_CLOSE, &
                        CLOSE_ELM(:,:,1:NEL_CLOSE), &
                        t2n_ELM, t2n_CLP)
    call cpu_time(clock2)
    if (ismaster) then
      write(*,*) "time (sec):", clock2-clock1, (clock2-clock1)/60.0d0
    end if

    ! check if the closest point is between -1 and 1
    if (maxval(abs(t2n_CLP(:,:,1:2))) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

    ! output
    open(nf2, file='t2n_close_point.dat', status='new')
    write(nf2,*) NRB%NEL, NRB%NGAUSS**2
    do i = 1, NRB%NEL
      do j = 1, NRB%NGAUSS**2
    write(nf2,'(I8,3F20.15)') t2n_ELM(i,j), t2n_CLP(i,j,:)
      end do
    end do
    close(nf2)

    NRB%CLE = t2n_ELM
    NRB%CLP = t2n_CLP(:,:,1:2)

    deallocate(CLOSE_ELM, t2n_ELM, t2n_CLP)
  end if
  close(nf1)

end subroutine find_close_point_NRB_TSP
