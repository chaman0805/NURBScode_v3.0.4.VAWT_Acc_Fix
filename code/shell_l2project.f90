!=======================================================================
! Main subroutine to drive L2 Projection (NURBS -> T-Spline)
!=======================================================================
subroutine n2t_l2project(TSP, BEZ, NRB, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, NRB
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(NRB%NNODE,NSD)
  real(8),    intent(out):: yg(TSP%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(TSP%NNODE+1), row(TSP%NNODE*50*TSP%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(TSP%NEL, TSP%NNODE, TSP%maxNSHL, TSP%IEN, &
                          TSP%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(TSP%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(TSP%NNODE))
  mg = 0.0d0

  call n2t_IntElmAss(TSP, BEZ, NRB, icnt, col, row, nsd, fg, &
                     RHSG, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(TSP%NNODE, TSP%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-15, 0)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine n2t_l2project


!=======================================================================
! Main subroutine to drive L2 Projection (T-Spline -> NURBS)
!=======================================================================
subroutine t2n_l2project(NRB, TSP, BEZ, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: NRB, TSP, BEZ
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(TSP%NNODE,NSD)
  real(8),    intent(out):: yg(NRB%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(NRB%NNODE+1), row(NRB%NNODE*50*NRB%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(NRB%NEL, NRB%NNODE, NRB%maxNSHL, NRB%IEN, &
                          NRB%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(NRB%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(NRB%NNODE))
  mg = 0.0d0

  call t2n_IntElmAss(NRB, TSP, BEZ, icnt, col, row, nsd, fg, &
                     RHSG, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(NRB%NNODE, NRB%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-15, 0)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine t2n_l2project



!=======================================================================
! Main subroutine to drive L2 Projection.
! This is a special type. Project pre-integrated right-hand-side (e.g.
! fluid forces) back to nodal points. Basically just devided by lumped
! mass or solve the linear system with consistent mass.
! fg: pre-integrated values stored at nodal points
! yg: nodal values
!=======================================================================
subroutine n2s_l2project(NRB, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: NRB
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(NRB%NNODE,NSD)
  real(8),    intent(out):: yg(NRB%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(NRB%NNODE+1), row(NRB%NNODE*50*NRB%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(NRB%NEL, NRB%NNODE, NRB%maxNSHL, NRB%IEN, &
                          NRB%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(NRB%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(NRB%NNODE))
  mg = 0.0d0

  ! RHS is the pre-integrated force
  RHSG = fg
  ! LHS is the consistent mass
  call n2s_IntElmAss(NRB, icnt, col, row, nsd, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(NRB%NNODE, NRB%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-15, 0)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine n2s_l2project


!=======================================================================
! Main subroutine to drive L2 Projection (FEM -> T-Spline)
!=======================================================================
subroutine f2t_l2project(TSP, BEZ, FEM, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: TSP, BEZ, FEM
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(FEM%NNODE,NSD)
  real(8),    intent(out):: yg(TSP%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(TSP%NNODE+1), row(TSP%NNODE*50*TSP%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(TSP%NEL, TSP%NNODE, TSP%maxNSHL, TSP%IEN, &
                          TSP%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(TSP%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(TSP%NNODE))
  mg = 0.0d0

  call f2t_IntElmAss(TSP, BEZ, FEM, icnt, col, row, nsd, fg, &
                     RHSG, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(TSP%NNODE, TSP%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-15, 0)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine f2t_l2project


!=======================================================================
! Main subroutine to drive L2 Projection (T-Spline -> FEM)
!=======================================================================
subroutine t2f_l2project(FEM, TSP, BEZ, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: FEM, TSP, BEZ
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(TSP%NNODE,NSD)
  real(8),    intent(out):: yg(FEM%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(FEM%NNODE+1), row(FEM%NNODE*50*FEM%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(FEM%NEL, FEM%NNODE, FEM%maxNSHL, FEM%IEN, &
                          FEM%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(FEM%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(FEM%NNODE))
  mg = 0.0d0

  call t2f_IntElmAss(FEM, TSP, BEZ, icnt, col, row, nsd, fg, &
                     RHSG, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(FEM%NNODE, FEM%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-15, 0)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine t2f_l2project



!=======================================================================
! Main subroutine to drive L2 Projection.
! This is a special type. Project pre-integrated right-hand-side (e.g.
! fluid forces) back to nodal points. Basically just devided by lumped
! mass or solve the linear system with consistent mass.
! fg: pre-integrated values stored at nodal points
! yg: nodal values
!=======================================================================
subroutine f2s_l2project(FEM, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: FEM
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(FEM%NNODE,NSD)
  real(8),    intent(out):: yg(FEM%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(FEM%NNODE+1), row(FEM%NNODE*50*FEM%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  write(*,*) 'Generate SparStruct'
  call genSparStruc_shell(FEM%NEL, FEM%NNODE, FEM%maxNSHL, FEM%IEN, &
                          FEM%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(FEM%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(FEM%NNODE))
  mg = 0.0d0

  ! RHS is the pre-integrated force
  RHSG = fg
  ! LHS is the consistent mass
  call f2s_IntElmAss(FEM, icnt, col, row, nsd, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0 
  call SparseCG_BDIAG_shell(FEM%NNODE, FEM%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-15, 0)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine f2s_l2project




!=======================================================================
! Main subroutine to drive L2 Projection (FEM2 -> FEM1)
!=======================================================================
subroutine f2f_l2project(FEM1, FEM2, NSD, fg, yg)
  use defs_shell
  implicit none

  type(mesh), intent(in) :: FEM1, FEM2
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(FEM2%NNODE,NSD)
  real(8),    intent(out):: yg(FEM1%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(FEM1%NNODE+1), row(FEM1%NNODE*50*FEM1%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(FEM1%NEL, FEM1%NNODE, FEM1%maxNSHL, FEM1%IEN, &
                          FEM1%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(FEM1%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(FEM1%NNODE))
  mg = 0.0d0

  call f2f_IntElmAss(FEM1, FEM2, icnt, col, row, nsd, fg, &
                     RHSG, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(FEM1%NNODE, FEM1%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, 1.0d-15, yg)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine f2f_l2project

!=======================================================================
! Main subroutine to drive L2 Projection (NURBS -> FEM)
!=======================================================================
subroutine n2f_l2project(FEM, NRB, NSD, fg, yg)

  use defs_shell
  implicit none

  type(mesh), intent(in) :: FEM, NRB
  integer,    intent(in) :: NSD
  real(8),    intent(in) :: fg(NRB%NNODE,NSD)
  real(8),    intent(out):: yg(FEM%NNODE,NSD)

  !  Local variables
  integer :: icnt, ier
  integer, allocatable :: col(:), row(:)
  real(8), allocatable :: RHSG(:,:), LHSK(:,:), mg(:)

  ! Build the sparse structure
  allocate(col(FEM%NNODE+1), row(FEM%NNODE*50*FEM%maxNSHL), stat=ier)
  if (ier /= 0) stop 'Allocation Error: col'        
  col = 0; row = 0; icnt = 0
        
  call genSparStruc_shell(FEM%NEL, FEM%NNODE, FEM%maxNSHL, FEM%IEN, &
                          FEM%NSHL, col, row, icnt)

  ! allocate the global RHS and LHS
  allocate(RHSG(FEM%NNODE,NSD), LHSK(NSD*NSD,icnt))
  RHSG = 0.0d0; LHSK = 0.0d0

  ! mg: Lumped mass
  allocate(mg(FEM%NNODE))
  mg = 0.0d0

  call n2f_IntElmAss(FEM, NRB, icnt, col, row, nsd, fg, &
                     RHSG, LHSK, mg)

  ! Do the L2-Projection with consistent mass
  ! Solve Mu = f
  yg = 0.0d0
  call SparseCG_BDIAG_shell(FEM%NNODE, FEM%maxNSHL, NSD, icnt, &
                            col, row, LHSK, RHSG, yg, 1.0d-10, 1)

  deallocate(col, row, RHSG, LHSK, mg)

end subroutine n2f_l2project
