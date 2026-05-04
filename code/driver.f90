!======================================================================
! Main routine to call all the subroutines                               
!======================================================================
program NURBScode  

  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell
  use types_beam
  use types_structure
  implicit none
  
  type(shell_bld) :: SH       ! shell for the blade  
  type(beam) :: BEA           ! beam
  type(mesh_mp_beam) :: mNRB_BEA
  type(mesh_beam)    :: NRB_BEA
  type(shell_nmb) :: NM       ! shell for non-matching (zones)
  type(structure)   :: STRU

  integer :: i, j, ii, k, ier, istep, Rstep, nn, dd, ibld, avgstepold, avgstep, b
  real(8) :: Rmt(3,3),    Rmdt(3,3),    Rmddt(3,3), xh_E_old(3)!,    &
  real(8) :: Center_Rot(3), &
             RmtOld(3,3), RmdtOld(3,3), RmddtOld(3,3), Rtmp(3,3)
!A  real(8), allocatable :: dshalpha(:,:)
  real(8), allocatable :: NRmat(:,:,:), NRdot(:,:,:), NRddt(:,:,:), &
                          NRmatOld(:,:,:), NRdotOld(:,:,:), NRddtOld(:,:,:), dshAlpha(:,:)

  character(len=30) :: iname,fname,cname

  ! Initialize MPI
  call MPI_INIT(                               mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid    , mpi_err)
  ismaster = myid.eq.mpi_master

  Center_Rot(1) = 0.0d0
  Center_Rot(2) = 0.0d0
  Center_Rot(3) = 0.0d0

!  solshell = (myid.eq.0)!.or.(myid.eq.1).or.(myid.eq.2)
 ! solshell = .false.
  NBlade = 1
  solshell = myid.eq.0

  ! flag for non-matching computation
  nonmatch = .true.

  ! Read mesh and MPI-communication Data
  if (ismaster) write(*,*) "Read mesh and communication data" 
  call input(myid+1)    
  if (numnodes > 1) call ctypes   

  !------------------------------------------------------------
  ! Only the first three processors are used to solve the blade 
  ! problem (there are three blades...)
  if (solshell) then
    ! Read shell mesh
    if (ismaster) write(*,*) "Read shell and beam mesh"
    call input_shell_blade(NSD, SH, BEA, mNRB_BEA, NRB_BEA)
    if (ismaster) write(*,*) "Finish reading shell and beam mesh"
  end if

  ! For the first three processors, get the number of
  ! fem shell nodes of the blades
  blade(:)%NNODE = 0 
  if (solshell) then
    blade(myid+1)%NNODE = SH%FEM%NNODE
  end if


  ! Broadcast them to all processors for the purpose of 
  ! allocating arrays
  do ibld = 1, NBlade
    call MPI_BCAST(blade(ibld)%NNODE, 1, MPI_INTEGER, &
                   ibld-1, MPI_COMM_WORLD, mpi_err)
  end do

  !----------------------------------------------------
  ! find how many surfaces numbered before the object
  !----------------------------------------------------
  SH%bmap = 0
  do i = 1, NBOUND
    if (bound(i)%Face_ID <= 20) SH%bmap = SH%bmap + 1
  end do


  !-------------------------------------------------------------
  ! non-matching
  !-------------------------------------------------------------
  if (nonmatch) then
    call input_shell_nmb(NSD, NM, Center_Rot)
  end if

  ! Get run parameters  
  if (ismaster) write(*,*) "Get run parameters"
  call getparam

  ! Get element gauss points (this needs to be after "getparam")
  allocate(ELMNGAUSS(NELEM))
  do i = 1, NELEM
    if (ELMNSHL(i) == 4) then
      ELMNGAUSS(i) = NGAUSSTET
    else if (ELMNSHL(i) == 6) then
      ELMNGAUSS(i) = NGAUSSPRI
    else
      write(*,*) "ERROR: Undefined ELMNSHL for ELMNGAUSS", ELMNSHL
    end if
  end do

  ! parameter for shell
  allocate(SH%Nnewt(NS_NL_itermax))
  SH%Nnewt = 1

  ! Generate Sparse Structures 
  if (ismaster) write(*,*) "Generating sparse structure"    
  call genSparStruc
!A allocate type structure which is combination of BEA and SH
  call allocate_LRhs(NRB_BEA, SH%NRB, NSD, STRU)
  
  ! Allocate Matrices and Vectors     
  if (ismaster) write(*,*) "Allocating matrices and vectors"
  call allocMatVec(SH, NRB_BEA, NM) !05/02/13 Allocate BEA and SH

!!!!!  call MPI_BARRIER (MPI_COMM_WORLD, mpi_err)
!!!!!  write(*,*) 'here', myid


!!$  ! Compute mass matrix and mesh data
!!$  if (ismaster) write(*,*) "Compute mass matrix and obtain meshsize"
!!$  call IntElmAss_mass
!!$  call FaceAssembly_area

  ! Read in restart files
  call readStep(Rstep)


!A--------------------S 05/06/13

    Identity = 0.0d0
    Identity(1,1) = 1.0d0
    Identity(2,2) = 1.0d0
    Identity(3,3) = 1.0d0

  if(solshell) then
    if (Rstep == 0) then
      RmatOld = Identity
      RtangOld = Identity
      OmegaOld = 0.0d0 !Full angular velocity 
      OmegatangOld = 0.0d0 !Tangential angular velocity 
    end if

      do i = 1, SH%NRB%NNODE
        if((abs(SH%NRB%B_NET(i,3)- 3.0d0) <=1.0d-4).and.(abs(SH%NRB%B_NET(i,1)- &
          0.0d0) <=1.0d-4).and.(abs(SH%NRB%B_NET(i,2)-0.0d0) <=1.0d-4)) then
          xrotOld(:) = SH%NRB%B_NET(i,1:3) !Initial center of rotation
        elseif((abs(SH%NRB%B_NET(i,3)- 9.0d0) <=1.0d-4).and.(abs(SH%NRB%B_NET(i,1)- &
        0.0d0) <=1.0d-4.and.(abs(SH%NRB%B_NET(i,2)-0.0d0) <=1.0d-4))) then
          xh_E_old(:) = SH%NRB%B_NET(i,1:3) !End of the hub
        end if
      end do
      write(*,*) 'XrotOld', xrotOld
      write(*,*) 'xh_E', xh_E_old
  end if
  
  if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  if (numnodes > 1) call MPI_BCAST(xrotOld, 3, MPI_DOUBLE_PRECISION, 0, &
                                   MPI_COMM_WORLD, mpi_err)
!A--------------------E 05/06/13)
  
  ! Get initial condition
  if (Rstep == 0) then   
    call generateIC(SH, NM, Center_Rot)
    call writeSol(Rstep)

    if(solshell) then
      do ii = 1, SH%NPT
        write(cname,'(I8.2)') ii
        write(iname,'(I8.2)') 0
        fname ='DISP/acel_SHL.'//trim(adjustl(iname))//'.'//trim(adjustl(cname))//'.dat'
        open(61, file=fname, status='replace')

        write(61,*) 0
        do i = 1, SH%mNRB%NNODE(ii)
          write(61,*) SH%NRB%FORCE(SH%mNRB%MAP(ii,i),1:3)
        end do
        close(61)
      end do


    do ii = 1, SH%NPT
      write(cname,'(I8.2)') ii
      write(iname,'(I8.2)') 0
      fname ='DISP/dis.'//trim(adjustl(iname))//'.'//trim(adjustl(cname))//'.dat'
      open(71, file=fname, status='replace')
      do i = 1, SH%mNRB%NNODE(ii)
        write(71,*) (SH%NRB%dsh(SH%mNRB%MAP(ii,i),k), k = 1, 4)
      end do
      write(71,*) 0
      close(71)
    end do

    do ii = 1, BEA%NPT
      write(cname,'(I8.2)') ii
      write(iname,'(I8.2)') 0
      fname ='DISP/beam.dis.'//trim(adjustl(iname))//'.'//trim(adjustl(cname))//'.dat'
      open(71, file=fname, status='replace')
      do i = 1, mNRB_BEA%NNODE(ii)
        write(71,*) (NRB_BEA%B_NET_D(mNRB_BEA%MAP(ii,i),k), k = 1, 4)
      end do
      write(71,*) 0
      close(71)
    end do
    end if ! endif solshell



 
!!!    call writeRB(Rstep)
  else
    call readSol(Rstep)
!!!    call readAvgSol(avgstepold)

    if (nonmatch) call readShellSol_NM(Rstep, NM)

    if (solshell) then     
!      SH%TSP%dshOld = 0.0d0; SH%FEM%dshOld = 0.0d0
!      SH%TSP%ushOld = 0.0d0; SH%FEM%ushOld = 0.0d0
!      SH%TSP%ashOld = 0.0d0; SH%FEM%ashOld = 0.0d0
!      SH%TSP%dsh = SH%TSP%dshOld; SH%FEM%dsh = SH%FEM%dshOld
!      SH%TSP%ush = SH%TSP%ushOld; SH%FEM%ush = SH%FEM%ushOld
!      SH%TSP%ash = SH%TSP%ashOld; SH%FEM%ash = SH%FEM%ashOld

      call readShellSol(Rstep, SH)
    end if
  end if ! end Rstep=0


  !------------------------------------------
  ! Loop over time steps
  !------------------------------------------
  avgstep = 0
  do istep = Rstep+1, Nstep

    avgstep = avgstep + 1

    time = time + Delt     

    if (ismaster) then
      write(*,'(60("="))')
      write(*,"(a,x,I8,x,ES14.6)") "Time Step Number:", istep, time
      write(*,'(60("="))')
    end if

    ! Predictor: Same Velocity
    ug  = ugold
    acg = (gami-1.0d0)/gami*acgold
    pg  = pgold
 
    phig  = phigold
    rphig = (gami-1.0d0)/gami*rphigold

    vbn1 = vbn0
    dbn1 = dbn0 + Delt*vbn0
    wbn1 = wbn0
!!!    call integrateRotation(Rn1,Rn0,wbn0,Delt)
!A----------------------------S 05/11/13
    if (solshell) then
!A Predictor for mesh
      ugm  = ugmold
      acgm = (gami-1.0d0)/gami*acgmold
      dg   = dgold + Delt*ugmold + &
             (gami-2.0d0*beti)/(2.0d0*gami)*Delt*Delt*acgmold

!A Predictor for SH%FEM

      SH%FEM%ush  = SH%FEM%ushold
      SH%FEM%ash = (gami-1.0d0)/gami*SH%FEM%ashold
      SH%FEM%dsh   = SH%FEM%dshold + Delt*SH%FEM%ushold + &
           (gami-2.0d0*beti)/(2.0d0*gami)*Delt*Delt*SH%FEM%ashold

!A Predictor for SH%NRB

      SH%NRB%ush  = SH%NRB%ushold
      SH%NRB%ash = (gami-1.0d0)/gami*SH%NRB%ashold
      SH%NRB%dsh   = SH%NRB%dshold + Delt*SH%NRB%ushold + &
           (gami-2.0d0*beti)/(2.0d0*gami)*Delt*Delt*SH%NRB%ashold

!A Acceleration
!A      if(time < Tacc) then
!A        thetd = maxthetd/Tacc*time
!A        thedd = maxthetd/Tacc
!A        theta = maxthetd/Tacc*time*time/2.0d0
!A      else
!A        thetd = -0.0d0
!A        thedd = 0.0d0
!A        theta = thetaold + Delt*thetdold + Delt*Delt/2.0d0* &
!A                         ((1.0d0-2.0d0*beti)*theddold+2.0d0*beti*thedd)
!A      end if

!A      call get_Rmat_Z(theta, thetd, thedd, Rmt, Rmdt, Rmddt)

!A      do i = 1, SH%NRB%NNODE
!A       do j = 1, 3
!A         if((abs(SH%NRB%B_NET_U(i,3)-90.0d0)<=1.0d-12).and.((SH%NRB%B_NET_U(i,1) )<4.0d0)) then

!A           SH%NRB%dsh(i,j) = sum(Rmt(j,:)*SH%NRB%B_NET(i,1:3))
!A           SH%NRB%ush(i,j) = sum(Rmdt(j,:)*SH%NRB%B_NET(i,1:3))
!A           SH%NRB%ash(i,j) = sum(Rmddt(j,:)*SH%NRB%B_NET(i,1:3))
     
!A         end if
!A       end do
!A      end do

      SH%NRB%B_NET_D(:,1:NSD) = SH%NRB%B_NET(:,1:NSD) + SH%NRB%dsh(:,1:3)


!A Predictor for BEA
      NRB_BEA%B_NET_Dt = SH%NRB%ush!NRB_BEA%B_NET_Dt_old
      NRB_BEA%B_NET_DDt = SH%NRB%ash!(gami-1.0d0)/gami*NRB_BEA%B_NET_DDt_old
      NRB_BEA%B_NET_D(:,1:NSD) = SH%NRB%B_NET_D(:,1:3)!NRB_BEA%B_NET_D_old(:,1:NSD)+Delt*NRB_BEA%B_NET_Dt_old(:,:) &
!A                                + Delt**2.0d0/2.0d0* ((1.0d0-2.0d0*beti)*NRB_BEA%B_NET_DDt_old(:,:) &
!A                                + 2.0d0*beti*NRB_BEA%B_NET_DDt(:,:))
!A    NRB_BEA%B_NET_D(:,NSD+1)=NRB_BEA%B_NET_D_old(:,NSD+1)
    end if !end solshel


!A--------------------S 05/06/13

    if(solshell) then
      call get_omega(SH%NRB,SH%Density,Omega, Omegatang)
    end if

    if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    if (numnodes > 1) call MPI_BCAST(Omega, 3, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)
    if (numnodes > 1) call MPI_BCAST(Omegatang, 3, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)
    if (numnodes > 1) call MPI_BCAST(xrot, 3, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)

    if(solshell) then
    Rfix = Identity
    if (nonmatch .and. mod(istep,1) == 3002) then
      if(ismaster) then
      norm = 0.0d0
        call fem_find_elm_normal(NM%FEM(1), NSD, 1, istep, norm(1,:))
        call fem_find_elm_normal(NM%FEM(2), NSD, 2, istep, norm(2,:))
      end if
      if (numnodes > 1) call MPI_BCAST(norm, 6, MPI_DOUBLE_PRECISION, 0, &
                                     MPI_COMM_WORLD, mpi_err)
      call fem_find_rot_mat(norm, Rfix)
    end if


    call integrateBodyRotation(Rmat,RmatOld, Identity,0.5d0*(Omega+OmegaOld),1.0d0/Delt)
    call integrateBodyRotation(Rtang,RtangOld, Identity,0.5d0*(Omegatang+OmegatangOld),1.0d0/Delt)

    if(ismaster) then
      write(*,*) 'Norm:', norm(1,:), norm(2,:)
      write(*,*) 'Omega:', Omega, 'Rmat:', Rmat
      write(*,*) 'Omegatang:', Omegatang, 'Rtang:', Rtang
    end if 
    
    if (nonmatch) then
      call setMeshBCs_tower !also includes predictor for NM, based on Rt matrix
      call setNMBCs(NM)
    end if
    
    else ! extract standard rotation matrix
!A----------------------------E 05/11/13

    !-- Algorithm for rotation ----------------------------------
    ! apply the exact rotation to the predictor, then solve
    ! for the deflection.
!A----------------------------S 05/06/13
    thetd = maxthetd
    thedd = 0.0d0
    theta = thetaold + Delt*thetdold + Delt*Delt/2.0d0* &
                       ((1.0d0-2.0d0*beti)*theddold+2.0d0*beti*thedd)

    call get_Rmat_Z(theta, thetd, thedd, Rmt, Rmdt, Rmddt)
    call get_Rmat_Z(thetaold, thetdold, theddold, RmtOld, RmdtOld, &
                  RmddtOld)


    allocate(NRmat(3,3,NNODE),NRdot(3,3,NNODE),NRddt(3,3,NNODE),&
             NRmatOld(3,3,NNODE),NRdotOld(3,3,NNODE),NRddtOld(3,3,NNODE))
    NRmat = 0.0d0; NRdot = 0.0d0; NRddt = 0.0d0
    NRmatOld = 0.0d0; NRdotOld = 0.0d0; NRddtOld = 0.0d0

    do nn = 1, NNODE
      if (NodeID(nn) >= 32) then
        NRmat(:,:,nn) = Rmt
        NRdot(:,:,nn) = Rmdt
        NRddt(:,:,nn) = Rmddt
        NRmatOld(:,:,nn) = RmtOld
        NRdotOld(:,:,nn) = RmdtOld
        NRddtOld(:,:,nn) = RmddtOld
      end if
    end do
 
    ! Predictor: rotation
    forall (nn = 1:NNODE, dd = 1:NSD)
      ugm(nn,dd) = ugmold(nn,dd) + &
                   sum((NRdot(dd,:,nn)-NRdotOld(dd,:,nn))*(xg(nn,1:3)-Center_Rot(1:3)))

      acgm(nn,dd) = sum(NRddt(dd,:,nn)*(xg(nn,1:3)-Center_Rot(1:3))) + &
                    (gami-1.0d0)/gami*(acgmold(nn,dd) - &
                     sum(NRddtOld(dd,:,nn)*(xg(nn,1:3)-Center_Rot(1:3))))

      dg(nn,dd) = dgold(nn,dd) + &
                  sum((NRmat(dd,:,nn)-NRmatOld(dd,:,nn))*(xg(nn,1:3)-Center_Rot(1:3))) + &
                 Delt*(ugmold(nn,dd)-sum(NRdotOld(dd,:,nn)*(xg(nn,1:3)-Center_Rot(1:3))))+ &
                  Delt*Delt/2.0d0*((1.0d0-2.0d0*beti)* &
                  (acgmold(nn,dd)-sum(NRddtOld(dd,:,nn)*(xg(nn,1:3)-Center_Rot(1:3)))) + &
                  2.0d0*beti*(acgm(nn,dd)-sum(NRddt(dd,:,nn)*(xg(nn,1:3)-Center_Rot(1:3)))))

    end forall

    deallocate(NRmat, NRdot, NRddt, NRmatOld, NRdotOld, NRddtOld)

! Will decide later what to do with FSI when NM is stationary
! Ideally when FSI NM is non-stationary

!A    ! Only the processors that solve the shell problem will need this
!A    ! prediction
!A    if (solshell) then
!A      ! The blade motion has to be the same as the corresponding
!A      ! fluid mesh, i.e. previous paragraph
!A      call predictor_rot_shell(SH%FEM, Rmat, Rdot, Rddt, RmatOld, &
!A                               RdotOld, RddtOld, Delt, beti, gami)
!A
!A      ! The blade motion has to be the same as the corresponding
!A      ! fluid mesh, i.e. previous paragraph
!A      call predictor_rot_shell(SH%NRB, Rmat, Rdot, Rddt, RmatOld, &
!A                               RdotOld, RddtOld, Delt, beti, gami)
!A    end if


    !-------------------------------------------------------
    ! for non-matching boundaries (only the inner part)
    !-------------------------------------------------------
    if (nonmatch) then
      call predictor_rot_shell(NM%FEM(1), Rmt, Rmdt, Rmddt, RmtOld, &
                               RmdtOld, RmddtOld, Delt, beti, gami, Center_Rot)

      ! the outer part is non-rotating
      Rtmp = 0.0d0
      call predictor_rot_shell(NM%FEM(2), Rtmp, Rtmp, Rtmp, Rtmp, &
                               Rtmp, Rtmp, Delt, beti, gami, Center_Rot)
    end if

    if (ismaster) then
      write(*,"(a,3ES15.6)") "Theta, Thetad, Thetadd =", theta, thetd, thedd
    end if

!A--------------------------------------------E 05/06/13

    !-- end rotation predictor ------------------------------------------------

    end if ! end solshell 
!A-----------------------------------E 05/06/13

!A    if(ismaster) then
!A            write(*,*) 'Centers'
!A	    write(*,*) xrot, xrotOld
!A            write(*,*) 'Rmat'
!A	    write(*,*) Rmat
!A	    write(*,*) RmatOld
!A            write(*,*) 'Rtan'
!A	    write(*,*) Rtang
!A	    write(*,*) RtangOld
!A            write(*,*) 'Omega'
!A	    write(*,*) Omega
!A	    write(*,*) OmegaOld
!A	    write(*,*) OmegaTang
!A	    write(*,*) OmegaTangOld
            !stop
!A    end if
! Mesh BC without NM
!A      call setMeshBCs_hawt

!A    do i = 1, NM%FEM(2)%NNODE
!A      do j = 1, 3
!A        if(ismaster) write(*,*) NM%FEM(2)%dsh(i,j)
!A        if((ismaster) .and. abs(NM%FEM(2)%dsh(i,j)) > 1e-8) write(*,*) 'Oops'
!A      end do
!A    end do
!A    do i = 1, NM%FEM(1)%NNODE
!A      do j = 1, 3
!A        if(ismaster) write(*,*) NM%FEM(1)%dsh(i,j)
!A        if((ismaster) .and. abs(NM%FEM(2)%dsh(i,j)) > 1e-8) write(*,*) 'Oops'
!A      end do
!A    end do
!    NM%FEM(1)%dsh(:,:) = -1.0d-6
!    NM%FEM(2)%dsh(:,:) = 1.0d-6


!A NM%FEM
    ! for non-matching 
    ! find the closest points (including finding closest elements)
    if (nonmatch) then
      do i = 1, 2
        allocate(dshAlpha(NM%FEM(i)%NNODE,NSD))
        dshAlpha = 0.0d0
        dshAlpha = NM%FEM(i)%dshOld + alfi*(NM%FEM(i)%dsh-NM%FEM(i)%dshOld)

        NM%FEM(i)%B_NET_D(:,1:3) = NM%FEM(i)%B_NET(:,1:3) + dshAlpha(:,:)
        deallocate(dshAlpha)
      end do

      do i = 1, 2
         NM%FEM(i)%CLE = 0; NM%FEM(i)%CLP = 0.0d0
      end do

      call find_close_point_FEM_FEM(NM%FEM(1), NM%FEM(2), NSD, istep, Rstep, 50)
      call find_close_point_FEM_FEM(NM%FEM(2), NM%FEM(1), NSD, istep, Rstep, 50)
    end if




    !--------------------------------------------
    ! Solve problems
    !--------------------------------------------
!!!    call setMeshBCs()
     
    ! Flags
    move = (time >= move_time)
    mono = (time >= mono_time).or.move
    conv = (time >= conv_time).or.mono

    shel = (time >= shel_time)      

    ! solve flow
    if (mono) then  
      call solflow_mono(istep)
    else
      call solflow_stag(istep, SH, NM, NRB_BEA, BEA, mNRB_BEA, Rstep, STRU)
    end if

    if (conv) then
      call redistance(istep)
      call massfix(istep) 
    end if  

    !--------------------------------------------
    ! Output results 
    !--------------------------------------------
    if (mod(istep,ifq_sh) == 0 .and. solshell) then

      do ii = 1, SH%NPT
        write(cname,'(I8.2)') ii
        write(iname,'(I8.2)') istep
        fname = 'DISP/acel_SHL.'//trim(adjustl(iname))//'.'//trim(adjustl(cname))//'.dat'
        open(61, file=fname, status='replace')
        write(61,*) istep
        do i = 1, SH%mNRB%NNODE(ii)
          write(61,*) SH%NRB%FORCE(SH%mNRB%MAP(ii,i),1:3)
        end do
        close(61)
      end do

    ! output displacement for structure mesh:
    do ii = 1, SH%NPT
      write(cname,'(I8.2)') ii
      write(iname,'(I8.2)') istep
      fname = 'DISP/dis.'//trim(adjustl(iname))//'.'//trim(adjustl(cname))//'.dat'
        open(71, file=fname, status='replace')
      do i = 1, SH%mNRB%NNODE(ii)
        write(71,*) (SH%NRB%dsh(SH%mNRB%MAP(ii,i),k), k = 1, 3)
      end do
      write(71,*) istep
      close(71)
    end do

    do ii = 1, BEA%NPT
      write(cname,'(I8.2)') ii
      write(iname,'(I8.2)') istep
      fname = 'DISP/beam.dis.'//trim(adjustl(iname))//'.'//trim(adjustl(cname))//'.dat'
        open(71, file=fname, status='replace')
      do i = 1, mNRB_BEA%NNODE(ii)
        write(71,*) (NRB_BEA%B_NET_D(mNRB_BEA%MAP(ii,i),k), k = 1, 4)
      end do
      write(71,*) istep
      close(71)
    end do

      call writeShellSol(istep, SH)
    end if

    if (mod(istep,ifq_tq) == 0 .and. ismaster) &
      call writetq(istep)

    if (mod(istep,ifq_tq) == 0 .and. solshell) &
      call writeTip(istep, SH)

    if (mod(istep,ifq) == 0 .and. ismaster .and. nonmatch) &
      call writeShellSol_NM(istep, NM)

    if (mod(istep,ifq) == 0) then
      call writeSol(istep)
!!!      call writeRB (istep)
    end if  

    !--------------------------------------------
    ! Update Old Quantities
    !--------------------------------------------
    call update_sol(SH, NM, NRB_BEA)

!    !--------------------------------------------
!    ! Get the averaged (in time) solutions
!    !--------------------------------------------
!    ! averaged relative velocity and pressure
!    do i = 1, NNODE
!      utmp = 0.0d0; umtmp = 0.0d0
!      call rot_vec_z(NSD, theta,  ug(i,:),  utmp)
!      call rot_vec_z(NSD, theta, ugm(i,:), umtmp)
!      uavg(i,:) = uavg(i,:) + (utmp-umtmp)
!    end do
!    pavg = pavg + pg
!    ! use the same ifq as solution outputs so that when restarted 
!    ! you won't average the same solution again
!    if (mod(istep,ifq)==0) call writeAvgSol(avgstepold+avgstep)

  end do
  
  ! Deallocate Matrices and Vectors     
  if (ismaster) write(*,*) "Deallocating matrices and vectors"
  call deallocMatVec

  ! Finalize MPI
  call MPI_FINALIZE(mpi_err)

end program NURBScode



!======================================================================
! subroutine to rotate the vector back to reference      
!======================================================================
subroutine rot_vec_z(nsd, theta, sol, rot)
  implicit none
  integer, intent(in)  :: nsd
  real(8), intent(in)  :: theta, sol(nsd)
  real(8), intent(out) :: rot(nsd)

  real(8) :: tmpx, tmpy

  tmpx = sol(1)
  tmpy = sol(2)
  rot(1) = cos(-theta)*tmpx - sin(-theta)*tmpy
  rot(2) = sin(-theta)*tmpx + cos(-theta)*tmpy
  rot(3) = sol(3)
end subroutine rot_vec_z



!======================================================================
! subroutine to update the solution
!======================================================================
subroutine update_sol(SH, NM, NRB_BEA)
  use aAdjKeep
  use commonvars
  use mpi
  use defs_shell
  use types_beam
  implicit none
 
  type(shell_bld), intent(inout) :: SH
  type(shell_nmb), intent(inout) :: NM
  type(mesh_beam), intent(inout) :: NRB_BEA
  integer :: i

  ! Update Old Quantities     
  dgold    = dg
  ugold    = ug
  acgold   = acg
  ugmold   = ugm
  acgmold  = acgm
  pgold    = pg
  rphigold = rphig
  phigold  = phig

  vbn0 = vbn1
  dbn0 = dbn1
  wbn0 = wbn1
  Rn0  = Rn1

  thetaOld = theta
  thetdOld = thetd
  theddOld = thedd

  OmegaOld     = Omega
  OmegatangOld = Omegatang
  RmatOld      = Rmat
  RtangOld     = Rtang

  if (solshell) then
!    SH%TSP%dshOld = SH%TSP%dsh
!    SH%TSP%ushOld = SH%TSP%ush
!    SH%TSP%ashOld = SH%TSP%ash

    SH%NRB%dshOld = SH%NRB%dsh
    SH%NRB%ushOld = SH%NRB%ush
    SH%NRB%ashOld = SH%NRB%ash

    SH%FEM%dshOld = SH%FEM%dsh
    SH%FEM%ushOld = SH%FEM%ush
    SH%FEM%ashOld = SH%FEM%ash

    NRB_BEA%B_NET_D_old = NRB_BEA%B_NET_D
    NRB_BEA%B_NET_Dt_old = NRB_BEA%B_NET_Dt
    NRB_BEA%B_NET_DDt_old = NRB_BEA%B_NET_DDt
  end if

  if (nonmatch) then
    do i = 1, 2
      NM%FEM(i)%dshOld = NM%FEM(i)%dsh
      NM%FEM(i)%ushOld = NM%FEM(i)%ush
      NM%FEM(i)%ashOld = NM%FEM(i)%ash
    end do
  end if
end subroutine update_sol



!======================================================================
! apply the exact rotation to the predictor
!======================================================================
subroutine predictor_rot_shell(SHL, Rmat, Rdot, Rddt, RmatOld, RdotOld, &
                               RddtOld, Delt, beti, gami, Center_Rot)
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: SHL
  real(8),    intent(in)    :: Rmat(3,3),    Rdot(3,3),    Rddt(3,3),   &
                               RmatOld(3,3), RdotOld(3,3), RddtOld(3,3),&
                               Delt, beti, gami, Center_Rot(3)
  integer :: nn, dd

  forall (nn = 1:SHL%NNODE, dd = 1:3)
    SHL%ush(nn,dd) = SHL%ushold(nn,dd)+  &
      sum((Rdot(dd,:)-RdotOld(dd,:))*(SHL%B_NET(nn,1:3)-Center_Rot(1:3)))

    SHL%ash(nn,dd) = sum(Rddt(dd,:)*(SHL%B_NET(nn,1:3)-Center_Rot(1:3))) + &
      (gami-1.0d0)/gami*(SHL%ashold(nn,dd) - &
      sum(RddtOld(dd,:)*(SHL%B_NET(nn,1:3)-Center_Rot(1:3))))

    SHL%dsh(nn,dd) = SHL%dshold(nn,dd) + &
          sum((Rmat(dd,:)-RmatOld(dd,:))*(SHL%B_NET(nn,1:3)-Center_Rot(1:3))) + &
          Delt*(SHL%ushold(nn,dd)-sum(RdotOld(dd,:)*(SHL%B_NET(nn,1:3)-Center_Rot(1:3))))+ &
          Delt*Delt/2.0d0*((1.0d0-2.0d0*beti)* &
          (SHL%ashold(nn,dd)-sum(RddtOld(dd,:)*(SHL%B_NET(nn,1:3)-Center_Rot(1:3)))) + &
          2.0d0*beti*(SHL%ash(nn,dd)-sum(Rddt(dd,:)*(SHL%B_NET(nn,1:3)-Center_Rot(1:3)))))
  end forall
end subroutine predictor_rot_shell

