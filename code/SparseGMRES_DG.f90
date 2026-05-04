!======================================================================
! SparseGMRES with non-local/matrix-free contributions
! Rhs1 is for local     contribution
! Rhs2 is for non-local contribution
! LHS for non-local contribution is NOT build.
! Non-local part is replaced with matrix-free method
!======================================================================
subroutine SparseGMRES_DG(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                          RHSGu, RHSGp, solu, solp, &
                          LHSK11, LHSG, LHSD1, LHSM, &
                          icntu, Utol, Kspaceu, Kspaceu_mn, &
                          NNODE, NSHL, NSD, &
                          dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                          acgmAlpha, pgAlpha, RHSGu2, RHSGp2, &
                          alfi, gami, Delt, NM) 
  use mpi
  use defs_shell
  implicit none

  type(shell_nmb), intent(inout) :: NM

  integer, intent(in) :: icntu, NNODE, NSHL, NSD, &
                         col(NNODE+1), row(icntu), &
                         IBC(NNODE,2*NSD+1), IPER(NNODE), &
                         D_FLAG(NNODE), P_FLAG(NNODE), &
                         Kspaceu, Kspaceu_mn

  real(8), intent(in) :: RHSGu(NNODE,NSD), RHSGp(NNODE), &
                         LHSK11(NSD*NSD,icntu), LHSG(NSD,icntu), &
                         LHSD1(NSD,icntu), LHSM(icntu), &
                         dgAlpha(NNODE,NSD), &
                         ugAlpha(NNODE,NSD),  ugmAlpha(NNODE,NSD), &
                         acgAlpha(NNODE,NSD), acgmAlpha(NNODE,NSD), &
                         pgAlpha(NNODE),  &
                         RHSGu2(NNODE,NSD), RHSGp2(NNODE), &
                         alfi, gami, Delt

  real(8), intent(inout) :: solu(NNODE,NSD), solp(NNODE)

  integer :: n, i, j, k, iK, iKs, jK, lK, itask, is, lenseg
  real(8) :: tmp2, abctmp(9)
    
  real(8) :: LHSKBdiagu(NNODE,NSD*NSD), LHSKBdiagm(NNODE,NSD*NSD), &
             LHSMdiag(NNODE)
  
  real(8) :: RHStmpu(NNODE,NSD), RHStmpm(NNODE,NSD), &
             RHStmpp(NNODE), temp1u(NNODE,NSD), &
             temp1m(NNODE,NSD), temp1p(NNODE), &
             RHSGu3(NNODE,NSD), RHSGp3(NNODE)            
  
  real(8) :: ugtemp(NNODE,NSD), pgtemp(NNODE)

  real(8) :: uBrgu(NNODE,NSD,Kspaceu+1), uBrgm(NNODE,NSD,Kspaceu+1), &
             uBrgp(NNODE,Kspaceu+1)
  
  real(8) :: HBrg(Kspaceu+1,Kspaceu), eBrg(Kspaceu+1), &
             yBrg(Kspaceu), &
             Rcos(Kspaceu), Rsin(Kspaceu), rr, unorm, &
             epsnrm, beta, ercheck, tmp, rr0, Utol, &
             Binv(NSD,NSD), flrate, unorm_ref, rrglob

  real(8), parameter :: RHSeps = 1.0d0

  ! Algorithm from Shakib/Johan's theses
  RHStmpu(:,:) = RHSGu(:,:)
  RHStmpp(:)   = RHSGp(:)
  
  ! Precondition Residual
  LHSKBdiagu = 0.0d0
  LHSMdiag   = 0.0d0
  
  do i = 1, NNODE     
    do j = col(i), col(i+1)-1
      n = row(j)
      if (n == i) then   
        LHSKBdiagu(i,:) = LHSK11(:,j)
        LHSMdiag(i)     = LHSM(j)      
      end if
    end do     
  end do  

  ! communicate the block-diagonal ! LGGL    
  if (numnodes .gt. 1) then
    call commu(LHSKBdiagu, NSD*NSD, 'in ')
    call commu(LHSMdiag,   1,       'in ')
     	 
!!!    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
	 
    call commu(LHSKBdiagu, NSD*NSD, 'out')
    call commu(LHSMdiag,   1,       'out')
  end if
  
  ! invert block-diagonal: momentum  
  Binv = 0.0d0
  do i = 1, NNODE
	  	   
    Binv(1,1) = LHSKBdiagu(i,5)*LHSKBdiagu(i,9) &
              - LHSKBdiagu(i,8)*LHSKBdiagu(i,6)
    Binv(1,2) = LHSKBdiagu(i,8)*LHSKBdiagu(i,3) &
              - LHSKBdiagu(i,2)*LHSKBdiagu(i,9)
    Binv(1,3) = LHSKBdiagu(i,2)*LHSKBdiagu(i,6) &
              - LHSKBdiagu(i,3)*LHSKBdiagu(i,5)

    tmp2 = Binv(1,1)*LHSKBdiagu(i,1) &
         + Binv(1,2)*LHSKBdiagu(i,4) &
         + Binv(1,3)*LHSKBdiagu(i,7)

    tmp = 1.0d0/tmp2
     
    Binv(1,1) = Binv(1,1)*tmp
    Binv(1,2) = Binv(1,2)*tmp
    Binv(1,3) = Binv(1,3)*tmp

    Binv(2,1) = (LHSKBdiagu(i,6)*LHSKBdiagu(i,7) &
               - LHSKBdiagu(i,4)*LHSKBdiagu(i,9)) * tmp
    Binv(2,2) = (LHSKBdiagu(i,1)*LHSKBdiagu(i,9) &
               - LHSKBdiagu(i,7)*LHSKBdiagu(i,3)) * tmp
    Binv(2,3) = (LHSKBdiagu(i,4)*LHSKBdiagu(i,3) &
               - LHSKBdiagu(i,1)*LHSKBdiagu(i,6)) * tmp
    Binv(3,1) = (LHSKBdiagu(i,4)*LHSKBdiagu(i,8) &
               - LHSKBdiagu(i,5)*LHSKBdiagu(i,7)) * tmp
    Binv(3,2) = (LHSKBdiagu(i,7)*LHSKBdiagu(i,2) &
               - LHSKBdiagu(i,1)*LHSKBdiagu(i,8)) * tmp
    Binv(3,3) = (LHSKBdiagu(i,1)*LHSKBdiagu(i,5) &
               - LHSKBdiagu(i,2)*LHSKBdiagu(i,4)) * tmp
     	  	 
    LHSKBdiagu(i,1) = Binv(1,1)
    LHSKBdiagu(i,2) = Binv(1,2)
    LHSKBdiagu(i,3) = Binv(1,3)
    LHSKBdiagu(i,4) = Binv(2,1)
    LHSKBdiagu(i,5) = Binv(2,2)
    LHSKBdiagu(i,6) = Binv(2,3)
    LHSKBdiagu(i,7) = Binv(3,1)
    LHSKBdiagu(i,8) = Binv(3,2)
    LHSKBdiagu(i,9) = Binv(3,3)   	  
  end do

  do i = 1, NNODE
    if (LHSMdiag(i) == 0.0d0) then
      LHSMdiag(i) = 1.0d0
    else   
      LHSMdiag(i) = 1.0d0/LHSMdiag(i)
    end if
  end do

  ! Precondition residual
  do i = 1, NNODE
    uBrgu(i,1,1) = LHSKBdiagu(i,1)*RHStmpu(i,1) + &
                   LHSKBdiagu(i,2)*RHStmpu(i,2) + &
                   LHSKBdiagu(i,3)*RHStmpu(i,3)
    uBrgu(i,2,1) = LHSKBdiagu(i,4)*RHStmpu(i,1) + &
                   LHSKBdiagu(i,5)*RHStmpu(i,2) + &
                   LHSKBdiagu(i,6)*RHStmpu(i,3)
    uBrgu(i,3,1) = LHSKBdiagu(i,7)*RHStmpu(i,1) + &
                   LHSKBdiagu(i,8)*RHStmpu(i,2) + &
                   LHSKBdiagu(i,9)*RHStmpu(i,3)

    uBrgp(i,1) = RHStmpp(i)*LHSMdiag(i)
  end do
    
  RHStmpu(:,:) = uBrgu(:,:,1)
  RHStmpp(:)   = uBrgp(:,1)
  
  !----------------------------------------------  
  ! calculate norm of RHS
  rr = sum(RHStmpp(:)*RHStmpp(:))
  do i = 1, NSD
    rr = rr + sum(RHStmpu(:,i)*RHStmpu(:,i))
  end do

  rrglob = rr

  if (numnodes > 1) then
    call MPI_ALLREDUCE(rrglob, rr, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if

  unorm = sqrt(rr)

!!!  if (ismaster) 
!!!     &     write(*,*) "After-Prec. Residual L_2 Norm is...", unorm

  iKs = 0
  
  ! set up tolerance
  epsnrm = Utol*unorm
  unorm_ref = 1.0d2/unorm
  
  ! set up RHS of the Hessenberg's problem
  eBrg(:) = 0.0d0
  eBrg(1) = unorm
  
  ! normalize the first Krylov vector
  uBrgu(:,:,1) = uBrgu(:,:,1) / unorm
  uBrgp(:,1)   = uBrgp(:,1)   / unorm
  
  ! loop through GMRES iterations
!!!  do 1000 iK = 1, Kspaceu
  do iK = 1, Kspaceu
  
    iKs = iK
     
    ! matrix-vector product
    temp1u(:,:) = uBrgu(:,:,iKs)
    temp1p(:)   = uBrgp(:,iKs)
	
    !----------------------------------------------------------------
    ! Periodicity (Slave = Master) - GL
    if (numnodes > 1) then
      call commu(temp1u, NSD, 'out')
      call commu(temp1p,   1, 'out')
    end if
         
!!$    do i = 1, NNODE
!!$      if ((IBC(i,1)==3).or.(IBC(i,2)==3).or.(IBC(i,3)==3)) then
!!$        temp1u(i,:) = temp1u(IPER(i),:) ! Slave = Master
!!$      end if 
!!$  
!!$      if (IBC(i,7)==3) then
!!$        temp1p(i)   = temp1p(IPER(i)) ! Slave = Master     
!!$      end if
!!$    end do

    ! Product
    call SparseProdUP_3D(col, row, LHSK11, LHSG, LHSD1, LHSM, &
                         temp1u, temp1p, uBrgu(:,:,iKs+1), &
                         uBrgp(:,iKs+1), D_FLAG, P_FLAG, &
                         NNODE, NSHL, icntu, NSD)
     
!!$    ! Communicate to Masters, Zero out Slaves - LG 
!!$    do i = 1, NNODE    
!!$      if ((IBC(i,1)==3).or.(IBC(i,2)==3).or.(IBC(i,3)==3)) then
!!$    uBrgu(IPER(i),:,iKs+1) = uBrgu(IPER(i),:,iKs+1) +
!!$     &           uBrgu(i,:,iKs+1) ! Master = Master+Slave
!!$    uBrgu(i,:,iKs+1) = 0.0d0 ! Slave = zero
!!$      end if
!!$
!!$      if (IBC(i,7)==3) then
!!$    uBrgp(IPER(i),iKs+1) = uBrgp(IPER(i),iKs+1) +
!!$     &             uBrgp(i,iKs+1) ! Master = Master+Slave
!!$    uBrgp(i,iKs+1) = 0d+0 ! Slave = zero
!!$      end if
!!$    end do
   
    !--------------------------------------
    ! matrix free for non-local part
    !---------------------------------------
!!!    if (iKs <= 50) then
      RHSGu3 = 0.0d0; RHSGp3 = 0.0d0
      ugtemp = ugAlpha + RHSeps*alfi*gami*Delt*temp1u
!!!      ugtemp = ugAlpha + RHSeps*gami*Delt*temp1u
      pgtemp = pgAlpha + RHSeps*alfi*gami*Delt*temp1p

      call FaceAsse_DG_w1t2(dgAlpha, ugtemp, ugmAlpha, pgtemp, &
                            RHSGu3, RHSGp3, NM)
!!!    end if     

    ! add non-local/matrix-free part
    RHStmpu = uBrgu(:,:,iKs+1) - (RHSGu3-RHSGu2)/RHSeps
    RHStmpp = uBrgp(:,iKs+1)   - (RHSGp3-RHSGp2)/RHSeps

    if (numnodes > 1) then     
      call commu(RHStmpu, NSD, 'in ')
      call commu(RHStmpp,   1, 'in ')
    end if

    !----------------------------------------------------------------------

!    RHStmpu(:,:) = uBrgu(:,:,iKs+1)
!    RHStmpp(:)   = uBrgp(:,iKs+1)
     
    ! Precondition product
    do i = 1, NNODE    
      uBrgu(i,1,iKs+1) = LHSKBdiagu(i,1)*RHStmpu(i,1) + &
                         LHSKBdiagu(i,2)*RHStmpu(i,2) + &
                         LHSKBdiagu(i,3)*RHStmpu(i,3)
      uBrgu(i,2,iKs+1) = LHSKBdiagu(i,4)*RHStmpu(i,1) + &
                         LHSKBdiagu(i,5)*RHStmpu(i,2) + &
                         LHSKBdiagu(i,6)*RHStmpu(i,3)
      uBrgu(i,3,iKs+1) = LHSKBdiagu(i,7)*RHStmpu(i,1) + &
                         LHSKBdiagu(i,8)*RHStmpu(i,2) + &
                         LHSKBdiagu(i,9)*RHStmpu(i,3)
        
      uBrgp(i,iKs+1) = RHStmpp(i)*LHSMdiag(i)
    end do
       
    ! orthogonalize and get the norm
    do jK = 1, iKs+1  

      if (jK == 1) then
     
        rr = sum(uBrgp(:,iKs+1)*uBrgp(:,1)) !{u_{i+1}*u_1} vector
        do i = 1, NSD
          rr = rr + sum(uBrgu(:,i,iKs+1)*uBrgu(:,i,1))
        end do
        rrglob = rr 
        if (numnodes > 1) then
          call MPI_ALLREDUCE(rrglob, rr, 1, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, &
                             MPI_COMM_WORLD, mpi_err)  
        end if
       
        beta = rr
!!!        if (ismaster) write(*,*) iKs, jK, beta   
      else
  
        ! project off jK-1 vector
        uBrgu(:,:,iKs+1) = uBrgu(:,:,iKs+1) - beta*uBrgu(:,:,jK-1)
        uBrgp(:,iKs+1)   = uBrgp(:,iKs+1)   - beta*uBrgp(:,jK-1)
     
        rr = sum(uBrgp(:,iKs+1)*uBrgp(:,jK)) !{u_{i+1}*u_j} vector
       
        do i = 1, NSD
          rr = rr + sum(uBrgu(:,i,iKs+1)*uBrgu(:,i,jK))
        end do
        rrglob = rr
        if (numnodes > 1) then
          call MPI_ALLREDUCE(rrglob, rr, 1, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, &
                             MPI_COMM_WORLD, mpi_err)
        end if
       
        beta = rr
!!!    if (ismaster) write(*,*) iKs, jK, beta   
      end if
    
      HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix

    end do    
     
    unorm = sqrt(beta)
    HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band
    
    ! normalize the Krylov vector
    uBrgu(:,:,iKs+1) = uBrgu(:,:,iKs+1) / unorm
    uBrgp(:,iKs+1)   = uBrgp(:,iKs+1)   / unorm
    
    ! construct and reduce the Hessenberg Matrix
    ! since there is only one subdiagonal we can use a Givens rotation
    ! to rotate off each subdiagonal AS IT IS FORMED. We do this because it
    ! allows us to check progress of solution and quit when satisfied.  Note
    ! that all future K vects will put a subdiagonal in the next column so
    ! there is no penalty to work ahead as  the rotation for the next vector
    ! will be unaffected by this rotation.
     
    ! H Y = E ========>   R_i H Y = R_i E
    do jK = 1, iKs-1
      tmp = Rcos(jK)*HBrg(jK,iKs) + Rsin(jK)*HBrg(jK+1,iKs)

      HBrg(jK+1,iKs) = -Rsin(jK)*HBrg(jK,  iKs) + &
                        Rcos(jK)*HBrg(jK+1,iKs)
      HBrg(jK,  iKs) = tmp
    end do
     
    tmp     = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
    Rcos(iKs)   = HBrg(iKs,  iKs) / tmp
    Rsin(iKs)   = HBrg(iKs+1,iKs) / tmp
    HBrg(iKs,  iKs) = tmp
    HBrg(iKs+1,iKs) = 0.0d0

    ! rotate eBrg    R_i E
    tmp = Rcos(iKs)*eBrg(iKs) + Rsin(iKs)*eBrg(iKs+1)
    eBrg(iKs+1) = -Rsin(iKs)*eBrg(iKs) + Rcos(iKs)*eBrg(iKs+1)
    eBrg(iKs)   =  tmp

    ! check for convergence
    ercheck = abs(eBrg(iKs+1))

    if ((ercheck <= epsnrm).and.(iKs >= Kspaceu_mn)) exit
        
    if (mod(iKs,100) == 0) then
      if (ismaster) then
        write(*,8000) iKs, ercheck, epsnrm, ercheck*unorm_ref
      end if
    end if 

  end do    ! end of GMRES iteration loop
!!! 1000 continue

  ! -- Solution -------------------------------------------------
  ! if converged or end of Krylov space
  ! solve for yBrg
  do jK = iKs, 1, -1
    yBrg(jK) = eBrg(jK)/HBrg(jK,jK)
    do lK = 1, jK-1
      eBrg(lK) = eBrg(lK) - yBrg(jK)*HBrg(lK,jK)
    end do
  end do
  
  ! update du, dp
  do jK = 1, iKs
    solu(:,:) = solu(:,:) + yBrg(jK)*uBrgu(:,:,jK)
    solp(:)   = solp(:)   + yBrg(jK)*uBrgp(:,jK)
  end do

  ! communicate solution - GL
  if (numnodes > 1) then 
    call commu(solu, NSD, 'out')
    call commu(solp,   1, 'out')
  end if

  !---------------------------------------------------------------

  if (ismaster) write(*,8500) iKs, ercheck*unorm_ref
  
  8000 format(3x,i4,') Residual, Goal, Reduction= (%):', &
              1x,e10.4,1x,e10.4,1x,F12.6)
  8500 format(10x,' Iterations:', 2x,i4,'  Reduction:', &
              2x,f10.6,"(%)",/)
    
end subroutine SparseGMRES_DG
