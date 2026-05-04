!========================================================================
! Main routine to call all the subroutines   
! Find closest points between FEM and T-Spline
!========================================================================
subroutine find_close_point_FEM_TSP(FEM, TSP, BEZ, NSD)
  use mpi
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: FEM, TSP, BEZ
  integer,    intent(in) :: NSD

  integer :: i, j, k, ier, nf1, nf2, nf3, itmp1, itmp2

  integer, allocatable :: CLOSE_ELM(:,:,:)
  integer :: NEL_CLOSE
    
  integer, allocatable :: f2t_ELM(:,:),   t2f_ELM(:,:)
  real(8), allocatable :: f2t_CLP(:,:,:), t2f_CLP(:,:,:)

  real(8) :: clock1, clock2

  character(len=30) :: fname, cname, ch1, ch2, fmat

  nf1 = 11; nf2 = 21; nf3 = 31

  !************************************************************
  ! Begin finding the closest points on FEM surface
  !************************************************************
  ! first check if the file exist. If exist, then read the 
  ! closest point from the files instead of building them
  write(cname,'(I8)') myid + 21
  fname = 'f2t_close_point.'//trim(adjustl(cname))//'.dat'
  open(nf1, file=fname, status='old', iostat=ier)
  if (ier == 0) then

    if (ismaster) then
      write(*,*) "f2t: Reading the closest-points list"
    end if

    read(nf1,*) itmp1, itmp2
    if (itmp1/=TSP%NEL .or. itmp2/=TSP%NGAUSS**2) then
      write(*,*) "ERROR: TSP%NEL or TSP%NGAUSS does not match the"
      write(*,*) "       numbers in f2t_close_point.dat"
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
      write(*,*) "f2t: Closest-points list does not exist"
      write(*,*) "     Build it:"
    end if

    !------------------------------------------------------------
    ! first read in (if exist) or build the list of closest 
    ! elements (in ascending order of distance)
    !------------------------------------------------------------
    write(cname,'(I8)') myid + 21
    fname = 'f2t_close_elm.'//trim(adjustl(cname))//'.dat'
    open(nf2, file=fname, status='old', iostat=ier)
    if (ier == 0) then

      if (ismaster) then
        write(*,*) "f2t: Reading the closest-elements list"
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
        write(*,*) "f2t: Closest-elements list does not exist"
        write(*,*) "     Build it:"
      end if

      ! Save the first 0.1% closest elements into the file
      NEL_CLOSE = int(FEM%NEL*0.001d0)
      if (ismaster) then
        write(*,*) "f2t: Number of closest elements:", NEL_CLOSE
      end if
      allocate(CLOSE_ELM(TSP%NEL,TSP%NGAUSS**2,NEL_CLOSE))
      CLOSE_ELM = 0

      ! we first check and save the 50% closest elements
      ! check the distance between each element's center (gp=0)
      call cpu_time(clock1)
      call f2t_find_elm(TSP, BEZ, FEM, nsd, NEL_CLOSE, CLOSE_ELM)
      call cpu_time(clock2)
      if (ismaster) then
        write(*,*) "time (sec):", clock2-clock1, (clock2-clock1)/60.0d0
      end if

      write(cname,'(I8)') myid + 21
      fname = 'f2t_close_elm.'//trim(adjustl(cname))//'.dat'
      open(nf3, file=fname, status='new')
      write(nf3,*) NEL_CLOSE, TSP%NEL

      write(ch1,'(I10)') FEM%NEL
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
    NEL_CLOSE = int(FEM%NEL*0.0005d0)
    if (ismaster) then
      write(*,*) "f2t: Number of closest elements:", NEL_CLOSE
    end if

    allocate(f2t_ELM(TSP%NEL,TSP%NGAUSS**2), &
             f2t_CLP(TSP%NEL,TSP%NGAUSS**2,3))
    f2t_ELM = 0; f2t_CLP = 0.0d0

    call cpu_time(clock1)
    call f2t_find_point(TSP, BEZ, FEM, nsd, NEL_CLOSE, &
                        CLOSE_ELM(:,:,1:NEL_CLOSE), &
                        f2t_ELM, f2t_CLP)
    call cpu_time(clock2)
    if (ismaster) then    
      write(*,*) "time (sec):", clock2-clock1, (clock2-clock1)/60.0d0
    end if

    ! check if the closest point is between -1 and 1
    if (maxval(f2t_CLP(:,:,1:2)) > 1.0d0 .or. &
        minval(f2t_CLP(:,:,1:2)) < 0.0d0) then
      write(*,*) "ERROR: Point located outside of 0 and 1"
      stop
    end if

    ! output
    write(cname,'(I8)') myid + 21
    fname = 'f2t_close_point.'//trim(adjustl(cname))//'.dat'
    open(nf2, file=fname, status='new')
    write(nf2,*) TSP%NEL, TSP%NGAUSS**2
    do i = 1, TSP%NEL
      do j = 1, TSP%NGAUSS**2
        write(nf2,'(I8,3F20.15)') f2t_ELM(i,j), f2t_CLP(i,j,:)
      end do
    end do
    close(nf2)

    TSP%CLE = f2t_ELM
    TSP%CLP = f2t_CLP(:,:,1:2)

    deallocate(CLOSE_ELM, f2t_ELM, f2t_CLP)
  end if
  close(nf1)


  !************************************************************
  ! Begin finding the closest points on T-Spline Surface
  !************************************************************
  ! first check if the file exist. If exist, then read the 
  ! closest point from the files instead of building them
  write(cname,'(I8)') myid + 21
  fname = 't2f_close_point.'//trim(adjustl(cname))//'.dat'
  open(nf1, file=fname, status='old', iostat=ier)
  if (ier == 0) then

    if (ismaster) then
      write(*,*) "t2f: Reading the closest-points list"
    end if

    read(nf1,*) itmp1, itmp2
    if (itmp1/=FEM%NEL .or. itmp2/=FEM%NGAUSS) then
      write(*,*) "ERROR: FEM%NEL or FEM%NGAUSS does not match the"
      write(*,*) "       numbers in f2t_close_point.dat"
      stop
    end if
    
    do i = 1, FEM%NEL
      do j = 1, FEM%NGAUSS
        read(nf1,*) FEM%CLE(i,j), FEM%CLP(i,j,:)
      end do
    end do

    ! check if the closest point is between -1 and 1
    if (maxval(abs(FEM%CLP)) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

  ! if the file does not exist, then build it
  else

    if (ismaster) then
      write(*,*) "t2f: Closest-points list does not exist"
      write(*,*) "     Build it:"
    end if

    !------------------------------------------------------------
    ! first read in (if exist) or build the list of closest 
    ! elements (in ascending order of distance)
    !------------------------------------------------------------
    write(cname,'(I8)') myid + 21
    fname = 't2f_close_elm.'//trim(adjustl(cname))//'.dat'
    open(nf2, file=fname, status='old', iostat=ier)
    if (ier == 0) then

      if (ismaster) then
        write(*,*) "t2f: Reading the closest-elements list"
      end if
      read(nf2,*) NEL_CLOSE
      allocate(CLOSE_ELM(FEM%NEL,FEM%NGAUSS,NEL_CLOSE))
      CLOSE_ELM = 0
      do i = 1, FEM%NEL
        do j = 1, FEM%NGAUSS
          read(nf2,*) CLOSE_ELM(i,j,1:NEL_CLOSE)
        end do
      end do

    else

      if (ismaster) then
        write(*,*) "t2f: Closest-elements list does not exist"
        write(*,*) "     Build it:"
      end if

      ! Save the first n% closest elements into the file
      NEL_CLOSE = int(TSP%NEL*0.02d0)
      if (ismaster) then
        write(*,*) "t2f: Number of closest elements:", NEL_CLOSE
      end if
      allocate(CLOSE_ELM(FEM%NEL,FEM%NGAUSS,NEL_CLOSE))
      CLOSE_ELM = 0

      ! we first check and save the 50% closest elements
      ! check the distance between each element's center (gp=0)
      call cpu_time(clock1)
      call t2f_find_elm(FEM, TSP, BEZ, nsd, NEL_CLOSE, CLOSE_ELM) 
      call cpu_time(clock2)
      if (ismaster) then
        write(*,*) "time (sec):", clock2-clock1, (clock2-clock1)/60.0d0
      end if

      write(cname,'(I8)') myid + 21
      fname = 't2f_close_elm.'//trim(adjustl(cname))//'.dat'
      open(nf3, file=fname, status='new')
      write(nf3,*) NEL_CLOSE, FEM%NEL

      write(ch1,'(I10)') TSP%NEL
      write(ch2,'(I10)') len(trim(adjustl(ch1)))+1
      write(ch1,'(I10)') FEM%NEL
      fmat = '('//trim(adjustl(ch1))//'I'//trim(adjustl(ch2))//')'

      do i = 1, FEM%NEL
        do j = 1, FEM%NGAUSS
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
      write(*,*) "t2f: Number of closest elements:", NEL_CLOSE
    end if

    allocate(t2f_ELM(FEM%NEL,FEM%NGAUSS), &
             t2f_CLP(FEM%NEL,FEM%NGAUSS,3))
    t2f_ELM = 0; t2f_CLP = 0.0d0

    call cpu_time(clock1)
    call t2f_find_point(FEM, TSP, BEZ, nsd, NEL_CLOSE, &
                        CLOSE_ELM(:,:,1:NEL_CLOSE), &
                        t2f_ELM, t2f_CLP)
    call cpu_time(clock2)
    if (ismaster) then
      write(*,*) "time (sec):", clock2-clock1, (clock2-clock1)/60.0d0
    end if

    ! check if the closest point is between -1 and 1
    if (maxval(abs(t2f_CLP(:,:,1:2))) > 1.0d0) then
      write(*,*) "ERROR: Point located outside of -1 and 1"
      stop
    end if

    ! output
    write(cname,'(I8)') myid + 21
    fname = 't2f_close_point.'//trim(adjustl(cname))//'.dat'
    open(nf2, file=fname, status='new')
    write(nf2,*) FEM%NEL, FEM%NGAUSS
    do i = 1, FEM%NEL
      do j = 1, FEM%NGAUSS
        write(nf2,'(I8,3F20.15)') t2f_ELM(i,j), t2f_CLP(i,j,:)
      end do
    end do
    close(nf2)

    FEM%CLE = t2f_ELM
    FEM%CLP = t2f_CLP(:,:,1:2)

    deallocate(CLOSE_ELM, t2f_ELM, t2f_CLP)
  end if
  close(nf1)

end subroutine find_close_point_FEM_TSP
