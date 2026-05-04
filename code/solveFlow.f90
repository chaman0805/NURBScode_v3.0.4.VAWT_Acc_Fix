!======================================================================
!
!======================================================================
subroutine solflow_stag(istep, SH, NM, NRB_BEA, BEA, mNRB_BEA , Rstep, STRU)
  
  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell
  use types_structure
  use types_beam

  implicit none 

  type(shell_bld), intent(inout) :: SH
  type(shell_nmb), intent(inout) :: NM
  type(mesh_beam), intent(in)    :: NRB_BEA
  type(beam), intent(in)	 :: BEA
  type(mesh_mp_beam), intent(in) :: mNRB_BEA
  type(structure)   :: STRU
  integer :: istep, Rstep, i
!A  real(8), allocatable :: dshalpha(:,:)

  integer :: inewt
  real(8) :: momres0, conres0, convres0, meshres0, value
  logical :: NSConverged, LSConverged
  
  momres0  = -1.0d0
  conres0  = -1.0d0
  convres0 = -1.0d0

  !------------------------------------------------
  ! Solve Navier-Stokes   
  !------------------------------------------------ 
  do inewt = 1, NS_NL_itermax

!A--------------------------E 05/11/13
    if (nonmatch) then
      call solveNavStoDG(inewt, momres0, conres0, NSConverged, NM)
    else
      call solveNavSto(inewt, momres0, conres0, NSConverged)
    end if

!A    do i = 1, NELEM
!A      value = sqrt(abs(ug(i,1)**2 + ug(i,2)**2 +ug(i,3)**2))
!A      if ((value) > 200.0d0) write(*,*) 'Value:', value,'Coordinates:', xg(i,1:3)
!A    end do

    if (NSConverged) exit
    if (shel) then
      call solveKLShell(SH, inewt, istep, STRU, NRB_BEA, BEA, mNRB_BEA)

!A-----------------------------------S 05/06/13
      if (solshell) call get_omega(SH%NRB,SH%Density,Omega, Omegatang)

    if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    if (numnodes > 1) call MPI_BCAST(Omega, 3, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)
    if (numnodes > 1) call MPI_BCAST(Omegatang, 3, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)
    if (numnodes > 1) call MPI_BCAST(xrot, 3, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)

    Rfix = Identity

    if (nonmatch .and. (inewt == 1)) then
      if(ismaster) then
        norm = 0.0d0
        call fem_find_elm_normal(NM%FEM(1), NSD, 1, istep, norm(1,:))
        call fem_find_elm_normal(NM%FEM(2), NSD, 2, istep, norm(2,:))
      end if
      if (numnodes > 1) call MPI_BCAST(norm, 6, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)
      call fem_find_rot_mat(norm, Rfix)
    end if

    if(ismaster) write(*,*) 'Norm:', norm(1,:), norm(2,:)

    call integrateBodyRotation(Rmat,RmatOld, Identity, 0.5d0*(Omega+OmegaOld),1.0d0/Delt)
    call integrateBodyRotation(Rtang,RtangOld, Rfix, 0.5d0*(Omegatang+OmegatangOld),1.0d0/Delt)

    if(ismaster) write(*,*) 'Omega:', Omega, 'Centers:', xrot, xrotOld, 'Rmat:', Rmat

!A-----------------------------------E 05/06/13
    if (nonmatch) then
      call setMeshBCs_tower
      call setNMBCs(NM)
    end if
!A NM%FEM We only do it once per time step at the driver
!A    if (nonmatch) then
!A      do i = 1, 2
!A        allocate(dshAlpha(NM%FEM(i)%NNODE,NSD))
!A        dshAlpha = 0.0d0
!A        dshAlpha = NM%FEM(i)%dshOld + alfi*(NM%FEM(i)%dsh-NM%FEM(i)%dshOld)

!A        NM%FEM(i)%B_NET_D(:,1:3) = NM%FEM(i)%B_NET(:,1:3) + dshAlpha(:,:)
!A        deallocate(dshAlpha)
!A      end do

!A      do i = 1, 2
!A         NM%FEM(i)%CLE = 0; NM%FEM(i)%CLP = 0.0d0
!A      end do

!A      call find_close_point_FEM_FEM(NM%FEM(1), NM%FEM(2), NSD, istep, Rstep, 1000)
!A      call find_close_point_FEM_FEM(NM%FEM(2), NM%FEM(1), NSD, istep, Rstep, 1000)
!A    end if
      call SolveMesh(inewt, meshres0)
    end if

  end do

  !------------------------------------------------
  ! Solve Levelset convection 
  !------------------------------------------------  
  if (conv) then
    if (ismaster) then
      write(*,*) "################################################"
      write(*,*) "Convecting levelset function "
      write(*,*) "################################################"
    end if
    do inewt = 1, LSC_NL_itermax    
      call solveLevelset(inewt, convres0, LSConverged)
      if (LSConverged) exit
    end do
  end if
  
end subroutine solflow_stag


!======================================================================
!
!======================================================================
subroutine solflow_mono(istep)
  
  use aAdjKeep
  use mpi
  use commonvars

  implicit none 
  
  integer :: istep
  integer :: inewt
  real(8) :: momres0, conres0, convres0, meshres0
  real(8) :: normFb0, normMb0
  logical :: NSConverged, LSConverged, RBMConverged
  
  write(*,*) "Are you sure you want to call solflow_mono???"
  stop

  momres0  = -1.0d0
  conres0  = -1.0d0
  convres0 = -1.0d0
  
  !------------------------------------------------     
  ! Get prediction move   
  !------------------------------------------------ 
  if (move) then
    call moveRBMesh(0, normFb0, normMb0, meshres0, RBMConverged)
  else
    RBMConverged = .true.
  end if  
  
  !------------------------------------------------ 
  ! Get prediction convection       
  !------------------------------------------------ 
  do inewt = 1, LSC_pred_step
    call solveLevelset(inewt, convres0, LSConverged)
    if (LSConverged) exit
  end do
  
  !------------------------------------------------ 
  ! Newton loop    
  !------------------------------------------------ 
  do inewt = 1, max(NS_NL_itermax,LSC_NL_itermax)
    if (mono_iter == 0) then
      call solveNavSto  (inewt, momres0, conres0, NSConverged)
      call solveLevelset(inewt, convres0,         LSConverged)
    else if (mono_iter == 1) then
      call solveNavStoLS(inewt, momres0, conres0, NSConverged, &
                         convres0, LSConverged)
    else
      write(*,*) "Undefined iteration procedure"
      write(*,*) "   mono_iter =", mono_iter
    end if
          
    if (move) call moveRBMesh(inewt, normFb0, normMb0, &
                              meshres0, RBMConverged)

    if (NSConverged.and.LSConverged.and.RBMConverged) exit

  end do
    
  !------------------------------------------------ 
  ! Move rigid body and mesh - if not converged     
  !------------------------------------------------   
  if (move.and.(.not.(NSConverged.and.RBMConverged))) then
    call moveRBMesh(inewt, normFb0, normMb0, meshres0, RBMConverged)
  end if
         
end subroutine solflow_mono
