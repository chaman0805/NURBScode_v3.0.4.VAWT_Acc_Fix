subroutine moveRBMesh(inewt, normFb0, normMb0, normFb, normMb,&
                      meshres0, RBMConverged)
  
  use aAdjKeep
  use mpi
  use commonvars

  implicit none 

  integer :: inewt
  real(8) :: normFb0,normMb0,meshres0
  logical :: RBMConverged

  real(8) :: alphaRB,dgAlpha(NNODE,NSD), &
             ugAlpha(NNODE,NSD), ugmAlpha(NNODE,NSD),&
             acgAlpha(NNODE,NSD), acgmAlpha(NNODE,NSD),&
             pgAlpha(NNODE),phigAlpha(NNODE)

  real(8) :: FbJ(NSD,NSD), MbJ(NSD,NSD)

  real(8) :: relFac(2),normFb,normMb,meshres,meshresl
    
  RBMConverged = .true.	    
! Get values at alpha level
  alphaRB = 0.5d0
  
  call setBCs()
  
  acgAlpha  = acgold  + alphaRB*(acg  - acgold)
  acgmAlpha = acgmold + alphaRB*(acgm - acgmold)
  ugAlpha   = ugold   + alphaRB*(ug   - ugold)
  ugmAlpha  = ugmold  + alphaRB*(ugm  - ugmold)
  dgAlpha   = dgold   + alphaRB*(dg   - dgold)
  pgAlpha   = pgold   + alphaRB*(pg   - pgold)

  phigAlpha = phigold + alphaRB*(phig - phigold)
    
! Assemble Force and Jacobian
  Fb = 0.0d0  
  Mb = 0.0d0 
  FbJ = 0.0d0 
  MbJ = 0.0d0 	 
  call ForceAssembly_3D(dgAlpha, ugAlpha, ugmAlpha,&
                        acgAlpha, acgmAlpha, pgAlpha, phigAlpha,&
                        Fb, Mb, FbJ, MbJ,alphaRB)
     
  Fb = Fb  + massb*gravvec    
  FbJ = 0d0   
  MbJ = 0d0
	    
  ! Move rigid body
  if (inewt.eq.1) then	  
    relFac(1)  = 0.5d0
    relFac(2)  = 0.5d0
  else
    relFac(1)  = 1d0/(1d0+(normFb/normFb0))
    relFac(2)  = 1d0/(1d0+(normMb/normMb0))   
  endif  
    
  call SingleRigidBody3D(massb,Ibhat,Rn0,Delt,Fb,Mb,FbJ,MbJ,&
                         vbn0,dbn0,wbn0,VBC,MBC,1,     &
                         vbn1,dbn1,wbn1,Rn1,normFb,normMb,relFac)

  ! Print Residual  		
  if (inewt.eq.1) then
    normFb0 = normFb
    normMb0 = normMb
  endif

  if (ismaster) then
    if (sum(VBC).ne.3) then
      write(*,'(I2,x,a,x,ES12.4,x,F12.6)') &
        inewt,") Rigid Body Lin. Mom. = ", &
        normFb,1d2*normFb/normFb0  
    endif  
    if (sum(MBC).ne.3) then
      write(*,'(I2,x,a,x,ES12.4,x,F12.6)') &
           inewt,") Rigid Body Ang. Mom. = ", &
           normMb,1d2*normMb/normMb0	
    endif 
  endif

  if (sum(VBC).ne.3) then
     if (normFb.ge.RB_NL_Ftol*normFb0) RBMConverged = .false. 
  endif  
  if (sum(MBC).ne.3) then
	 if (normMb.ge.RB_NL_Mtol*normMb0) RBMConverged = .false. 
  endif 
	    
!-------------------------------------------------------------------------
! MESH SOLVE
!-------------------------------------------------------------------------
  call setMeshBCs()
    	    
  acgAlpha  = acgold  + almi*(acg  - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha   = ugold   + alfi*(ug   - ugold)
  ugmAlpha  = ugmold  + alfi*(ugm  - ugmold)
  dgAlpha   = dgold   + alfi*(dg   - dgold)
  pgAlpha   = pg

  ! Assemble Interior
  RHSGm  = 0d0	 
  lhsK22 = 0d0
  call IntElmMesh_3D(dgAlpha)

  ! Compute Residual   
  if (numnodes .gt. 1) call commu(RHSGm, NSD, 'in ')
       	       
  meshresl = sum(RHSGm(:,1)*RHSGm(:,1)) &
           + sum(RHSGm(:,2)*RHSGm(:,2)) &
           + sum(RHSGm(:,3)*RHSGm(:,3))
    
  meshres = meshresl
  if (numnodes.gt.1) then		
     call MPI_ALLREDUCE(meshresl, meshres, 1, &
                        MPI_DOUBLE_PRECISION,MPI_SUM,  &
                        MPI_COMM_WORLD,mpi_err)	  
  endif

  meshres  = sqrt(meshres)

  ! Print Residual  
  if (inewt.eq.1) meshres0 = meshres

  if (ismaster) then
    write(*,'(I2,x,a,x,ES12.4,x,F12.6)') &
        inewt,") Mesh   Res. Norm = ",  &
        meshres,  1d2*meshres/meshres0
  end if
  
  ! Solve linear system
  if (meshres.lt.Mesh_NL_tol*meshres0) then
    if (ismaster) then
      write(*,*)
      write(*,*) "     Converged: skip linear solve"
      write(*,*) 
    end if 
  else
    RBMConverged = .false. 
    acgmAlpha = 0.0d0
    
    call SparseGMRES_m(col, row, &
                       IBC, IPER, D_FLAG, P_FLAG, rhsGm, &
                       acgmAlpha, lhsK22, icnt, &
                       Mesh_GMRES_tol, Mesh_GMRES_itermax,Mesh_GMRES_itermin, &
                       NNODE, maxNSHL, NSD)

    acgm = acgm + acgmAlpha
    ugm  = ugm  + gami*Delt*acgmAlpha
    dg   = dg   + beti*Delt*Delt*acgmAlpha
  end if	

  if (numnodes.gt.1) call MPI_BARRIER (MPI_COMM_WORLD,mpi_err)

end subroutine moveRBmesh
