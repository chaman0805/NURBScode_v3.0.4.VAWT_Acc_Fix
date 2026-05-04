subroutine allocate_LRhs(NRB_BEA, NRB_SHL, NSD, STRU)

  use types_beam
  use types_shell
  use types_structure
  implicit none

  type(mesh_beam),  intent(in)    :: NRB_BEA
  type(mesh_shell),  intent(in)    :: NRB_SHL
  type(structure)   :: STRU

  integer,     intent(in)    :: NSD
  integer :: ier

    allocate(STRU%col(NRB_SHL%NNODE+1), STRU%row(NRB_SHL%NNODE*50*max(NRB_SHL%maxNSHL,NRB_BEA%maxNSHL)), stat=ier)
    if (ier /= 0) stop 'Allocation Error: col'    
    STRU%col = 0; STRU%row = 0  
  !write (*,*) 'icnt=', STRU%icnt  
    call genSparStruc(NRB_BEA%NEL, NRB_BEA%NNODE, NRB_BEA%maxNSHL, NRB_BEA%IEN, NRB_BEA%NSHL, &
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
