!======================================================================
! MESH SOLVE
!======================================================================
subroutine SolveMesh(inewt, meshres0)
  
  use aAdjKeep
  use mpi
  use commonvars

  implicit none 

  integer, intent(in)    :: inewt
  real(8), intent(inout) :: meshres0

  real(8) :: dgAlpha(NNODE,NSD), meshres, meshresl

!  if (ismaster) then
!    write(*,*)
!    write(*,'(I3,a)') inewt,") Solving Mesh Problem ---"
!    write(*,*)
!  end if

!A-------------------S 05/11/13 move this before calling fluid solver
!A  if (nonmatch) then
!A    call setMeshBCs_tower
!A  else
!A    call setMeshBCs_hawt
!A  end if

  dgAlpha = dgold + alfi*(dg-dgold)

  ! Assemble Interior
  RHSGm  = 0.0d0
  lhsK22 = 0.0d0

  call IntElmMesh_3D(dgAlpha)

  ! Compute Residual   
  if (numnodes > 1) call commu(RHSGm, NSD, 'in ')

  meshresl = 0.0d0
  meshresl = sum(RHSGm(:,1)**2) + sum(RHSGm(:,2)**2) &
           + sum(RHSGm(:,3)**2)
    
  meshres = meshresl
  if (numnodes.gt.1) then
     call MPI_ALLREDUCE(meshresl, meshres, 1, &
                        MPI_DOUBLE_PRECISION,MPI_SUM,  &
                        MPI_COMM_WORLD,mpi_err)	  
  endif

  meshres = sqrt(meshres)

  ! Print Residual  
  if (inewt.eq.1) meshres0 = meshres

  if (ismaster) then
    write(*,*)
    write(*,'(I3,a,x,ES13.6,x,F12.6)') &
        inewt,") Mesh   Res. Norm = ",  &
        meshres,  1.0d2*meshres/meshres0
  end if
  
  ! Solve linear system
  if (meshres < Mesh_NL_tol*meshres0) then
    if (ismaster) then
      write(*,*)
      write(*,*) "      Converged: skip linear solve"
      write(*,*) 
    end if 

  else if (meshres < 1.0d-15) then
    if (ismaster) then
      write(*,*)
      write(*,*) "      Zero Mesh Residual: skip mesh solve"
      write(*,*) 
    end if 

  else
    dgAlpha = 0.0d0
    
    call SparseGMRES_m(col, row, &
                       IBC, IPER, D_FLAG, P_FLAG, rhsGm, &
                       dgAlpha, lhsK22, icnt, &
                       Mesh_GMRES_tol, Mesh_GMRES_itermax,Mesh_GMRES_itermin, &
                       NNODE, maxNSHL, NSD)

    acgm = acgm + dgAlpha
    ugm  = ugm  + gami*Delt*dgAlpha
    dg   = dg   + beti*Delt*Delt*dgAlpha
  end if

  if (numnodes > 1) call MPI_BARRIER (MPI_COMM_WORLD, mpi_err)

end subroutine SolveMesh
