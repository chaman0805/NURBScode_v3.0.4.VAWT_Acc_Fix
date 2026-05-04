!======================================================================
!
!======================================================================
subroutine solveNavStoDG(inewt, MomRes0, ConRes0, NSConverged, NM)  
  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell
  implicit none 

  type(shell_nmb), intent(inout) :: NM

  integer, intent(in)    :: inewt
  real(8), intent(inout) :: MomRes0, ConRes0
  logical, intent(out)   :: NSConverged
 
  integer :: i
 
  real(8) :: dgAlpha(NNODE,NSD),                        &
             ugAlpha(NNODE,NSD),  ugmAlpha(NNODE,NSD),  &
             acgAlpha(NNODE,NSD), acgmAlpha(NNODE,NSD), &
             pgAlpha(NNODE),      phigAlpha(NNODE)
  
  real(8) :: MomRes, ConRes, MomResl, ConResl, meshresl, &
             solu(NNODE,NSD), solp(NNODE),               &
             RHSGu2(NNODE,NSD), RHSGp2(NNODE)             

  real(8) :: vtmp(NSD), TRACT(NNODE,NSD), rtmp 

  character(len=80) :: fomt
  
  ! Set boundary conditions
  call setBCs_CFD()
!!$  call setBCs()

  ! Get quantities at alpha levels (generalized-alpha):
  acgAlpha  = acgold  + almi*(acg  - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha   = ugold   + alfi*(ug   - ugold)
  ugmAlpha  = ugmold  + alfi*(ugm  - ugmold)
  dgAlpha   = dgold   + alfi*(dg   - dgold)
  pgAlpha   = pg

  phigAlpha = phigold + alfi*(phig - phigold)

  ! Assemble Navier-Stokes LHS/RHS
  RHSGu = 0.0d0;
  RHSGp = 0.0d0 
  
  LHSK11 = 0.0d0  
  LHSG   = 0.0d0
  LHSD1  = 0.0d0
  LHSM   = 0.0d0

  RHSGu2 = 0.0d0
  RHSGp2 = 0.0d0 

  call IntElmAss_3D(inewt, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                    acgmAlpha, pgAlpha, phigAlpha, 0)

  call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, 0)

  ! For non-matching
  call FaceAsse_DG_w1t1(dgAlpha, ugAlpha, ugmAlpha, pgAlpha)

  ! For non-local (projected) contributions
  call FaceAsse_DG_w1t2(dgAlpha, ugAlpha, ugmAlpha, pgAlpha, &
                        RHSGu2, RHSGp2, NM)

  ! Compute residuals       
  if (numnodes > 1) then
    call commu(RHSGp, 1  , 'in ')
    call commu(RHSGu, NSD, 'in ')

    call commu(RHSGp2, 1  , 'in ')
    call commu(RHSGu2, NSD, 'in ')
  end if

  RHSGp = RHSGp + RHSGp2
  RHSGu = RHSGu + RHSGu2

  MomResl = sum(RHSGu(:,1)**2) + sum(RHSGu(:,2)**2) + &
            sum(RHSGu(:,3)**2)

  ConResl = sum(RHSGp**2)
  
  MomRes = MomResl
  ConRes = ConResl

  if (numnodes > 1) then
    call MPI_ALLREDUCE(MomResl, MomRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
    call MPI_ALLREDUCE(ConResl, ConRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if

  MomRes = sqrt(MomRes)
  ConRes = sqrt(ConRes)

  ! Print residual  
  if (inewt == 1) then
    MomRes0 = MomRes     
    ConRes0 = ConRes 
  end if
   
  if (ismaster) then
    fomt = "(I3,a,x,ES13.6,x,F12.6)"
    write(*,fomt) inewt, ") Total Mom. Res. Norm = ", &
                  MomRes, 1.0d2*MomRes/MomRes0

    write(*,fomt) inewt, ") Continuity Res. Norm = ", &
                  ConRes, 1.0d2*ConRes/ConRes0
   
    write(*,*) " "
  end if

  if (isnan(MomRes+ConRes)) then
    write(*,*) '!=== Norm is NaN ===!'
    stop
  end if

  ! Check convergence and solve linear system
  if (((MomRes/MomRes0) < NS_NL_Utol).and. &
      ((ConRes/ConRes0) < NS_NL_Ptol) ) then
    NSConverged = .true.
    if (ismaster) then
      write(*,*)
      write(*,*) "    Converged: skip linear solve"
      write(*,*)
    end if

  else
    NSConverged = .false.
    solu = 0.0d0; solp = 0.0d0

    call SparseGMRES_DG(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        RHSGu, RHSGp, solu, solp, &
                        LHSK11, LHSG, LHSD1, LHSM, &
                        icnt, NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
                        NNODE, maxNSHL, NSD, &
                        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                        acgmAlpha, pgAlpha, RHSGu2, RHSGp2, &
                        alfi, gami, Delt, NM)

    acg = acg + solu
    ug  = ug  + gami*Delt*solu   
    pg  = pg  + alfi*gami*Delt*solp
  end if
  
  !------------------------------------------------------------------
  ! compute RHS and torque (weak BC/Surface Int.)
  !------------------------------------------------------------------
  acgAlpha  = acgold  + almi*(acg  - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha   = ugold   + alfi*(ug   - ugold)
  ugmAlpha  = ugmold  + alfi*(ugm  - ugmold)
  dgAlpha   = dgold   + alfi*(dg   - dgold)
  pgAlpha   = pg
  phigAlpha = phigold + alfi*(phig - phigold)

  ! compute RHS and the torque for individual blade
!  call compute_RHS_PhVI(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!                        acgmAlpha, pgAlpha, phigAlpha)

  ! compute individual torque
  call compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                   acgmAlpha, pgAlpha, phigAlpha, 21,    &
                   rtmp, torque1)

!A  call compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!A                   acgmAlpha, pgAlpha, phigAlpha, 22,    &
!A                   rtmp, torque2)

!A  call compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!A                   acgmAlpha, pgAlpha, phigAlpha, 23,    &
!A                   rtmp, torque3)

  !---------------------------------------------------------------
  ! compute RHS and the torque (RHSGu is used afterwards for FSI)
!A  call compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!A                   acgmAlpha, pgAlpha, phigAlpha, 4,     &  ! 1-nrel; 4-5mw
!A                   rtmp, torque4)

  if (numnodes > 1) call MPI_BARRIER (MPI_COMM_WORLD, mpi_err)

end subroutine solveNavStoDG
