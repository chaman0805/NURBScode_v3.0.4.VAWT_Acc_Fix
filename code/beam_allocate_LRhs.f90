subroutine beam_allocate_LRhs(FUN, NSD, BEA)

  use types_beam
  implicit none

  type(beam), intent(inout) :: BEA
  type(mesh),  intent(in)    :: FUN
  integer,     intent(in)    :: NSD
  integer :: ier

    allocate(BEA%col(FUN%NNODE+1), BEA%row(FUN%NNODE*50*FUN%maxNSHL), stat=ier)
    if (ier /= 0) stop 'Allocation Error: col'    
    BEA%col = 0; BEA%row = 0    
    call genSparStruc_beam(FUN%NEL, FUN%NNODE, FUN%maxNSHL, FUN%IEN, &
                            FUN%NSHL, BEA%col, BEA%row, BEA%icnt)

  ! allocate the global RHS and LHS
  allocate(BEA%RHSG(FUN%NNODE,NSD), &
           BEA%RHSG_EXT(FUN%NNODE,NSD), &
           BEA%RHSG_GRA(FUN%NNODE,NSD))
  BEA%RHSG     = 0.0d0
  BEA%RHSG_EXT = 0.0d0
  BEA%RHSG_GRA = 0.0d0

    allocate(BEA%LHSK(NSD*NSD,BEA%icnt))
    BEA%LHSK = 0.0d0

  allocate(BEA%dg(FUN%NNODE,NSD),  &
           BEA%tg(FUN%NNODE),     BEA%mg(FUN%NNODE),  &
           BEA%dl(FUN%NNODE,NSD), BEA%ddu(FUN%NNODE,NSD) )
  BEA%dg = 0.0d0; BEA%tg = 0.0d0; BEA%mg = 0.0d0
  BEA%dl = 0.0d0; BEA%ddu = 0.0d0

end subroutine beam_allocate_LRhs
