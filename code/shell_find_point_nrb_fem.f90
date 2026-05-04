subroutine find_close_point_nrb_fem(NRB, FEM)

  use mpi
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: NRB, FEM
  integer :: mf, iel, ct, igauss, jgauss
  character(len=30) :: cname, fname

  mf = 12
  write(cname,'(I8)') myid+21
  fname = 'f2n_close_point.'//trim(adjustl(cname))//'.dat'
  open(mf, file=fname, status='old')
    read(mf,*) ct
    do iel = 1, NRB%NEL
      ct = 0
      do igauss = 1, NRB%NGAUSS
        do jgauss = 1, NRB%NGAUSS
          ct = ct + 1
          read(mf,*) NRB%CLE(iel,ct), NRB%CLP(iel,ct,1:2)
        end do
      end do
     end do
  close(mf)

  write(cname,'(I8)') myid+21
  fname = 'n2f_close_point.'//trim(adjustl(cname))//'.dat'
  open(mf, file=fname, status='old')
    read(mf,*) ct
    do iel = 1, FEM%NEL
      do igauss = 1, FEM%NGAUSS
        read(mf,*) FEM%CLE(iel,igauss), FEM%CLP(iel,igauss,1:2)
      end do
     end do
  close(mf)
end subroutine find_close_point_nrb_fem
