!======================================================================
!
!======================================================================  
subroutine allocMatVec(SH, NRB_BEA, NM)
  
  use aAdjKeep
  use commonvars
  use defs_shell
  use types_beam

  implicit none
    
  type(shell_bld), intent(inout) :: SH
  type(mesh_beam), intent(inout) :: NRB_BEA
  type(shell_nmb), intent(inout) :: NM
  integer :: NGAUSSMAX, i

  ! Allocate solution vectors       
  allocate(dg(NNODE,NSD), dgold(NNODE,NSD))

  allocate( ug(NNODE,NSD),  ugold(NNODE,NSD), &
           ugm(NNODE,NSD), ugmold(NNODE,NSD))
   
  allocate( acg(NNODE,NSD),  acgold(NNODE,NSD), &
           acgm(NNODE,NSD), acgmold(NNODE,NSD))

  allocate(phig(NNODE),   phigold(NNODE), &
           rphig(NNODE),  rphigold(NNODE))
  
  allocate(pg(NNODE), pgold(NNODE))
  
  allocate(FPKS(NELEM,NSD,NSD))

  ! allocate averaged solution vectors
  allocate(uavg(NNODE,NSD), pavg(NNODE))

  ! Allocate Residuals     
  allocate(RHSGu(NNODE,NSD), RHSGm(NNODE,NSD), &
           RHSGp(NNODE), RHSGls(NNODE)) 
     
     
  ! Allocate Matrices
  allocate(lhsK11(NSD*NSD,icnt))
  allocate(lhsK12(NSD*NSD,icnt),lhsK22(NSD*NSD,icnt))

  allocate(lhsG(NSD,icnt), lhsD1(NSD,icnt),  &
           lhsD2(NSD,icnt))

  allocate(LHSM(icnt), LHSLS(icnt), LHSmass(icnt))

  allocate(LHSlsu(NSD,icnt),LHSPls(icnt),LHSUls(NSD,icnt))

  NGAUSSMAX = maxval(ELMNGAUSS)

  !--------------------------------------------
  ! allocation for shell
  !--------------------------------------------
  if (solshell) then
    ! allocate the shell RHS and LHS
!  allocate(SH%RHSG(SH%TSP%NNODE,NSD), &
!           SH%RHSG_EXT(SH%TSP%NNODE,NSD), &
!           SH%RHSG_GRA(SH%TSP%NNODE,NSD))


  allocate(SH%RHSG(SH%NRB%NNODE,NSD), &
           SH%RHSG_EXT(SH%NRB%NNODE,NSD), &
           SH%RHSG_GRA(SH%NRB%NNODE,NSD))

    allocate(SH%LHSK(NSD*NSD,SH%icnt))
    SH%RHSG = 0.0d0; SH%RHSG_EXT = 0.0d0; SH%RHSG_GRA = 0.0d0
    SH%LHSK = 0.0d0

    ! Allocate solution variables for SHELL
    allocate(SH%FEM%dsh(SH%FEM%NNODE,NSD), SH%FEM%dshOld(SH%FEM%NNODE,NSD), &
             SH%FEM%ush(SH%FEM%NNODE,NSD), SH%FEM%ushOld(SH%FEM%NNODE,NSD), &
             SH%FEM%ash(SH%FEM%NNODE,NSD), SH%FEM%ashOld(SH%FEM%NNODE,NSD))
    SH%FEM%dsh = 0.0d0; SH%FEM%dshOld = 0.0d0
    SH%FEM%ush = 0.0d0; SH%FEM%ushOld = 0.0d0
    SH%FEM%ash = 0.0d0; SH%FEM%ashOld = 0.0d0

!  allocate(SH%TSP%dsh(SH%TSP%NNODE,NSD), SH%TSP%dshOld(SH%TSP%NNODE,NSD), &
!           SH%TSP%ush(SH%TSP%NNODE,NSD), SH%TSP%ushOld(SH%TSP%NNODE,NSD), &
!           SH%TSP%ash(SH%TSP%NNODE,NSD), SH%TSP%ashOld(SH%TSP%NNODE,NSD))
!  SH%TSP%dsh = 0.0d0; SH%TSP%dshOld = 0.0d0
!  SH%TSP%ush = 0.0d0; SH%TSP%ushOld = 0.0d0
!  SH%TSP%ash = 0.0d0; SH%TSP%ashOld = 0.0d0

    allocate(SH%NRB%dsh(SH%NRB%NNODE,NSD), SH%NRB%dshOld(SH%NRB%NNODE,NSD), &
           SH%NRB%ush(SH%NRB%NNODE,NSD), SH%NRB%ushOld(SH%NRB%NNODE,NSD), &
           SH%NRB%ash(SH%NRB%NNODE,NSD), SH%NRB%ashOld(SH%NRB%NNODE,NSD))
    SH%NRB%dsh = 0.0d0; SH%NRB%dshOld = 0.0d0
    SH%NRB%ush = 0.0d0; SH%NRB%ushOld = 0.0d0
    SH%NRB%ash = 0.0d0; SH%NRB%ashOld = 0.0d0
  end if

  !--------------------------------------------
  ! Allocate solution variables for Non-matching boundaries
  !--------------------------------------------
  if (nonmatch) then
    do i = 1, 2
      allocate(NM%FEM(i)%dsh(NM%FEM(i)%NNODE,NSD), NM%FEM(i)%dshOld(NM%FEM(i)%NNODE,NSD), &
               NM%FEM(i)%ush(NM%FEM(i)%NNODE,NSD), NM%FEM(i)%ushOld(NM%FEM(i)%NNODE,NSD), &
               NM%FEM(i)%ash(NM%FEM(i)%NNODE,NSD), NM%FEM(i)%ashOld(NM%FEM(i)%NNODE,NSD))
      NM%FEM(i)%dsh = 0.0d0; NM%FEM(i)%dshOld = 0.0d0
      NM%FEM(i)%ush = 0.0d0; NM%FEM(i)%ushOld = 0.0d0
      NM%FEM(i)%ash = 0.0d0; NM%FEM(i)%ashOld = 0.0d0
    end do
  end if
end subroutine allocMatVec 

!======================================================================
!
!======================================================================  

subroutine allocate_LRhs(NRB_BEA, NRB_SHL, NSD, STRU)

  use types_beam
  use defs_shell
  use types_structure
  implicit none

  type(mesh_beam),  intent(in)    :: NRB_BEA
  type(mesh),  intent(in)    :: NRB_SHL
  type(structure)   :: STRU

  integer,     intent(in)    :: NSD
  integer :: ier

    allocate(STRU%col(NRB_SHL%NNODE+1), STRU%row(NRB_SHL%NNODE*50*max(NRB_SHL%maxNSHL,NRB_BEA%maxNSHL)), stat=ier)
    if (ier /= 0) stop 'Allocation Error: col'    
    STRU%col = 0; STRU%row = 0  
  !write (*,*) 'icnt=', STRU%icnt  
    call genSparStruc_STRU(NRB_BEA%NEL, NRB_BEA%NNODE, NRB_BEA%maxNSHL, NRB_BEA%IEN, NRB_BEA%NSHL, &
                           NRB_SHL%NEL, NRB_SHL%NNODE, NRB_SHL%maxNSHL, NRB_SHL%IEN, NRB_SHL%NSHL, &
                           STRU%col, STRU%row, STRU%icnt)
  
  !write (*,*) 'icnt=', STRU%icnt
  ! allocate the global RHS and LHS
  allocate(STRU%RHSG(NRB_SHL%NNODE,NSD), &
           STRU%RHSG_EXT(NRB_SHL%NNODE,NSD), &
           STRU%RHSG_GRA(NRB_SHL%NNODE,NSD))
  STRU%RHSG     = 0.0d0
  STRU%RHSG_EXT = 0.0d0
  STRU%RHSG_GRA = 0.0d0

  allocate(STRU%LHSK(NSD*NSD,STRU%icnt))
  STRU%LHSK = 0.0d0

  allocate(STRU%dg(NRB_SHL%NNODE,NSD),  &
           STRU%tg(NRB_SHL%NNODE),     STRU%mg(NRB_SHL%NNODE),  &
           STRU%dl(NRB_SHL%NNODE,NSD), STRU%ddu(NRB_SHL%NNODE,NSD) )
  STRU%dg = 0.0d0; STRU%tg = 0.0d0; STRU%mg = 0.0d0
  STRU%dl = 0.0d0; STRU%ddu = 0.0d0

end subroutine allocate_LRhs




!======================================================================
!
!======================================================================  
subroutine deallocMatVec 
  
  use aAdjKeep
  use mpi
  use commonvars

  implicit none
  
  integer :: i

  ! deallocate mesh and flags     
  deallocate(xg, wg, IEN)
  deallocate(IPER,IBC)

  deallocate(EL_TYP,P_Flag,D_Flag)
  
  if (NPATCH.gt.0) deallocate(EPID, EIJK)
  do i = 1, NPATCH  
    deallocate(patch(i)%U_KNOT)
    deallocate(patch(i)%V_KNOT)
    deallocate(patch(i)%W_KNOT)
  end do  
  deallocate(patch) 
  
  do i = 1, NBOUND  
    deallocate(bound(i)%FACE_IEN)  
    deallocate(bound(i)%F2E)
    deallocate(bound(i)%FACE_OR)     
    deallocate(bound(i)%BNODES)   
  end do    
  deallocate(bound)  
  
  ! deallocate sparseStruct
  deallocate(col,row)
       
  ! dellocate solution vectors       
  deallocate(dg,dgold)
  deallocate(ug,ugold,ugm,ugmold)   
  deallocate(acg,acgold,acgm,acgmold)
  deallocate(phig,phigold,rphig,rphigold)  
  deallocate(pg, pgold)
  deallocate(uavg, pavg)
  deallocate(FPKS)
  
  ! deallocate Residuals     
  deallocate(RHSGu, RHSGm,RHSGp, RHSGls) 
      
  ! deallocate Matrices
  deallocate(lhsK11,lhsK12,lhsK22)
  deallocate(lhsG, lhsD1,lhsD2)
  deallocate(LHSM, LHSLS, LHSmass)
  deallocate(LHSlsu,LHSPls,LHSUls)
  
  ! Conditional deallocations  
  if (numnodes.gt.1) deallocate(ilworku)

end subroutine deallocMatVec 
