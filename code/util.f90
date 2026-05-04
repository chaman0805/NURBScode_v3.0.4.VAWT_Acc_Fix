!======================================================================
! 
!======================================================================
subroutine writeStep(Rstep)
  use mpi
  implicit none

  integer, intent(in) :: Rstep  
  
#if USEMPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)  
#endif
  
  if (ismaster) then
    open (223, file='step.dat', status='replace')
    write(223,*) Rstep
    close(223)
  endif
  
end subroutine writeStep
    
!======================================================================
! 
!======================================================================
subroutine getStep(Rstep)
  
  use mpi
  implicit none
  
  integer, intent(out) :: Rstep
  integer :: ierr
  logical :: stat
  character(len=10) :: cname 
    
  ! Read step on master
  if (ismaster) then 
    open(222, file='step.dat', status='old', iostat=ierr)
    if (ierr == 0) then
      read (222,*) Rstep
      close(222) 
    else  
      Rstep = 0     
    endif  
  endif
      
  ! Communicate    
#if USEMPI
  call MPI_BCAST(Rstep, 1,   MPI_INTEGER, &
                 mpi_master, MPI_COMM_WORLD, mpi_err)
#endif

  ! Check restarts
  inquire(file="restart."//trim(cname(Rstep))//'.'//trim(cname(myid+1)), exist=stat)
  if (.not.stat) Rstep = 0

  ! Synchronize
#if USEMPI
  ierr= Rstep
  call MPI_ALLREDUCE(ierr, Rstep, 1,   MPI_INTEGER, &
                     MPI_MIN, MPI_COMM_WORLD, mpi_err)
#endif 
 
end subroutine getStep

!======================================================================
! Open a file for a global time series
!======================================================================
subroutine openStreamFile (fnum,fname)

 use params  
 use mpi
 implicit none 

 character*(*) :: fname
 integer :: fnum, eof
 real(8) :: ftime
 logical :: stat

 inquire(file=fname, exist=stat)

 if (stat) then
   write(*,*) "Opening existing file: ", fname
   open(fnum, file=fname, status='old')
   ftime = 0d0       
   do while (ftime+1.001d0*Delt < time) 
      read(fnum,*, IOSTAT=eof) ftime
      if (eof < 0) exit   
   enddo 
  else
   write(*,*) "Create new file: ", fname
   open(fnum, file=fname, status='replace')
 endif
      
end subroutine openStreamFile

!======================================================================
! Fucntion converts integer in corresponding character
!======================================================================
function cname(i)

  implicit none
  
  character(len=10) :: cname
  integer, intent(in) :: i

  integer :: il(0:8), ic0, ii, k
  logical :: beg
  character(len=10) :: cc
  
  ic0 = ICHAR("0")
  cc = " "
  ii = i
  
  il(0) = mod(ii,10)
  do k = 1,8
    ii = (ii - il(k-1)) / 10
    il(k) = mod (ii,10)
  enddo
  
  beg = .false.
  
  do k = 8,1,-1
    if (il(k) .ne. 0 .or. beg) then
      beg = .true.
      cc  = TRIM(cc) // CHAR(ic0 + il(k))
    endif
  enddo
  
  cc = TRIM(cc)//CHAR(ic0 + il(0))
  cname = cc
  
end function cname

!----------------------------------------------------------------------
!  Function for levi-civita symbol ==> for outerproduct
!---------------------------------------------------------------------- 
function Eijk(i,j,k)
  implicit none
  real(8)  Eijk
  integer :: i,j,k
  Eijk = real((j-i)*(k-i)*(k-j)/2,8)
end function Eijk
 
!----------------------------------------------------------------------
!  Function checks whether id is part of array
!---------------------------------------------------------------------- 
logical function  contains(size, array, id)

  implicit none
  integer :: size, id
  integer :: array(size)

  integer :: gdof
  integer :: rowvl, rowvh, rowv

  rowvl = 1
  rowvh = size + 1

  do 
    if (rowvh-rowvl > 1) then
      rowv = (rowvh + rowvl) / 2
      if (array(rowv) > id) then 
        rowvh = rowv
      else  
        rowvl = rowv
      endif
    else
      contains = .false.
      if (array(rowvl).eq.id) contains = .true.
      !if (array(rowvh).eq.id) contains = .true.      
      exit
    endif
  enddo

end function contains

!---------------------------------------------------------------------- 
logical function  findIndex(size, array, id,idx)

  implicit none
  integer :: size, id, idx
  integer :: array(size)
  
  findIndex = .true.

  do idx = 1, size
    if (array(idx) == id)  return
  enddo
  
  idx = -1
  findIndex = .false.

end function findIndex

!---------------------------------------------------------------------- 
! subroutine to invert a 3x3 matrix
!---------------------------------------------------------------------- 
subroutine get_inverse_3x3(Amat, Ainv, DetJ)
  implicit none
  real(8), intent(in)  :: Amat(3,3)
  real(8), intent(out) :: Ainv(3,3), DetJ
  real(8) :: tmp

  Ainv = 0.0d0
  Ainv(1,1) = Amat(2,2)*Amat(3,3) - Amat(3,2)*Amat(2,3)
  Ainv(1,2) = Amat(3,2)*Amat(1,3) - Amat(1,2)*Amat(3,3)
  Ainv(1,3) = Amat(1,2)*Amat(2,3) - Amat(1,3)*Amat(2,2)
  
  tmp = 1.0d0/(Ainv(1,1)*Amat(1,1) + Ainv(1,2)*Amat(2,1) + &
               Ainv(1,3)*Amat(3,1))

  Ainv(1,1) = Ainv(1,1) * tmp
  Ainv(1,2) = Ainv(1,2) * tmp
  Ainv(1,3) = Ainv(1,3) * tmp

  Ainv(2,1) = (Amat(2,3)*Amat(3,1) - Amat(2,1)*Amat(3,3)) * tmp
  Ainv(2,2) = (Amat(1,1)*Amat(3,3) - Amat(3,1)*Amat(1,3)) * tmp
  Ainv(2,3) = (Amat(2,1)*Amat(1,3) - Amat(1,1)*Amat(2,3)) * tmp
  Ainv(3,1) = (Amat(2,1)*Amat(3,2) - Amat(2,2)*Amat(3,1)) * tmp
  Ainv(3,2) = (Amat(3,1)*Amat(1,2) - Amat(1,1)*Amat(3,2)) * tmp
  Ainv(3,3) = (Amat(1,1)*Amat(2,2) - Amat(1,2)*Amat(2,1)) * tmp
  
  DetJ = 1.0d0/tmp  
  
end subroutine get_inverse_3x3
