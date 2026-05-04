module elasticBody
  !-------------------------------------------------------------------------------------------------------  
  type bodyState
    real(8), allocatable :: lambda (:), dlambda(:) , Fa(:)
    real(8) :: vb(3),db(3),wb(3),Qb(3,3), Fd(3), Fw(3)    
    
    real(8) :: Iij(3,3)
    real(8), allocatable :: Nai(:,:)
    
    integer :: NMODES = 0   
    integer :: MPItype
    integer :: force = -1
    
  end type bodyState 

  !-------------------------------------------------------------------- 
  type bodyParams
    real(8) :: mass
    real(8) :: xcg(3)
    real(8) :: Jij(3,3)
  
    integer :: VBC(3),MBC(3),EBC 
    
    integer  NBOUND
    integer, allocatable :: Bound_ID (:)
    integer, allocatable :: DofMap   (:)
        
    integer ::  NMODES,NDOFMODE,NFACEMODE

    integer, allocatable :: IEN_MODE (:,:)   
    real(8), allocatable :: xg_mode  (:,:) 
    integer, allocatable :: modemap  (:)   
    real(8), allocatable :: modes    (:,:,:) 
    
    real(8), allocatable :: Mabij(:,:,:,:)  
    real(8), allocatable :: Saij (:,:,:)
    real(8), allocatable :: Ka   (:)  
    real(8), allocatable :: Ma   (:)
  end type bodyParams  
  
  !--------------------------------------------------------------------  
  contains      
  !--------------------------------------------------------------------
  subroutine allocBodyState(st,NMODES)
    use mpi
    implicit none
    integer, intent(in) ::  NMODES
    type(bodyState), intent(inout) :: st  

    if (st%NMODES /= NMODES) then  
      call deallocBodyState(st) 
      allocate(st%lambda(NMODES),st%dlambda(NMODES),st%Fa(NMODES))       
      allocate(st%Nai (NMODES,3))
      st%NMODES = NMODES 
    endif
    
!#if USEMPI
!    call genBodyType(st)  
!#endif

  end subroutine allocBodyState
  
  !--------------------------------------------------------------------
  subroutine deallocBodyState(st)
    implicit none
    type(bodyState), intent(inout) :: st  

    if (st%NMODES /= 0) deallocate(st%lambda,st%dlambda,st%Fa,st%Nai) 
    st%NMODES = 0

  end subroutine deallocBodyState  
   
  !--------------------------------------------------------------------
#if USEMPI 
  subroutine genBodyType(st)
  
    use mpi 
    implicit none
  
    type(bodyState), intent(inout) :: st
    integer blocklength(6),typelist(6), reloc(6), add0,add
    integer ierr
    
    ! Generate   
    blocklength(1:3) = 3
    blocklength(4)   = 9
    blocklength(5:6) = st%NMODES
  
    typelist  = MPI_DOUBLE_PRECISION

    call MPI_Address(st, add0, mpi_err)  
     
    call MPI_Address(st%vb, add, mpi_err)   
    reloc(1) = add - add0
 
    call MPI_Address(st%db, add, mpi_err)   
    reloc(2) = add - add0
  
    call MPI_Address(st%wb, add, mpi_err)   
    reloc(3) = add - add0 
    
    call MPI_Address(st%Qb, add, mpi_err)   
    reloc(4) = add - add0   

    if (st%NMODES /= 0) then 
      call MPI_Address(st%lambda, add, mpi_err)  
      reloc(5) = add - add0

      call MPI_Address(st%dlambda, add, mpi_err)  
      reloc(6) = add - add0    
      call MPI_Type_struct(6, blocklength, reloc, typelist, st%MPItype, mpi_err)  
    else
      call MPI_Type_struct(4, blocklength, reloc, typelist, st%MPItype, mpi_err)        
    endif
     
    call MPI_Type_commit(st%MPItype, mpi_err)
            
  end subroutine genBodyType
#endif  

  !--------------------------------------------------------------------   
#if USEMPI 
  subroutine syncBodyState(st)    
    use mpi 
    implicit none
 
    type(bodyState), intent(inout) :: st
           
    !call MPI_BCAST(st,1,st%MPItype, &
    !               mpi_master,MPI_COMM_WORLD,mpi_err)
                   
    call MPI_BCAST(st%vb,3,MPI_DOUBLE_PRECISION,mpi_master,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(st%db,3,MPI_DOUBLE_PRECISION,mpi_master,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(st%wb,3,MPI_DOUBLE_PRECISION,mpi_master,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(st%Qb,9,MPI_DOUBLE_PRECISION,mpi_master,MPI_COMM_WORLD,mpi_err)
    if (st%NMODES /= 0) then
      call MPI_BCAST(st%lambda,st%NMODES,MPI_DOUBLE_PRECISION,mpi_master,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(st%dlambda,st%NMODES,MPI_DOUBLE_PRECISION,mpi_master,MPI_COMM_WORLD,mpi_err)
    endif
                     
  end subroutine syncBodyState   
#endif
   
  !--------------------------------------------------------------------   
  subroutine initBodyState(st)
    implicit none

    type(bodyState), intent(inout) :: st
 
    st%Fd = 0d0  
    st%vb = 0d0 
    st%db = 0d0

    st%Fw = 0d0
    st%wb = 0d0     
    st%Qb = 0d0    
    st%Qb(1,1) = 1d0
    st%Qb(2,2) = 1d0
    st%Qb(3,3) = 1d0
    
    if (st%NMODES /= 0 ) then     
      st%Fa      = 0d0
      st%lambda  = 0d0   
      st%dlambda = 0d0         
    endif
                
  end subroutine initBodyState  
  !--------------------------------------------------------------------       
  subroutine readBodyState(st, Rstep, b) 
    use mpi

    implicit none
  
    type(bodyState), intent(inout) :: st
    integer,         intent(in)    :: Rstep, b 

    integer :: sfile, ierr, step1, step2
    integer :: NMODES
       
    character(len=30) :: fname
    character(len=10) :: cname1  
    character(len=10) :: cname2
    if (ismaster) then 
      sfile = 15
      
      write(cname1,'(I8)') b
      write(cname2,'(I8)') Rstep

      fname = 'body.' // trim(adjustl(cname1))//'.'// trim(adjustl(cname2))
      open(sfile, file=fname, status='old', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "File does not exist:", fname, " ==> revert to default"
        call initBodyState(st)
             
      else
    
      ! Read header
      read(sfile,*) step1
      read(sfile,*) NMODES
    
      if ((step1  /= Rstep   ).or. &
          (NMODES /= st%NMODES  ) ) then
        write(*,*) "ERROR IN READING", fname
        write(*,*) " -STEP  = ",  step1,Rstep   
        write(*,*) " -MODES = ",  NMODES,st%NMODES 
        stop
      endif 
       
      ! Read mesh function data        
      st%Fd = 0d0
      read(sfile,*) st%vb
      read(sfile,*) st%db
      st%Fw = 0d0
      read(sfile,*) st%wb
      read(sfile,*) st%Qb

      if (NMODES /= 0) then
        st%Fa = 0d0
        read(sfile,*) st%dlambda 
        read(sfile,*) st%lambda
      endif
    
      read(sfile,*) step2
      close(sfile)
  
      if (step2 /= Rstep   ) then
        write(*,*) "ERROR IN READING", fname
        write(*,*) " -STEP  = ",  Rstep , step2  
        write(*,*) " FILE END CHECK FAILED"
        stop
      endif
      endif    
    endif

#if USEMPI 
    call syncBodyState(st)    
#endif    
 
  end subroutine readBodyState
  
  !--------------------------------------------------------------------      
  subroutine writeBodyState(st, Rstep, b)
    use mpi
    
    implicit none
  
    type(bodyState), intent(in) :: st
    integer,         intent(in) :: Rstep, b 

    integer :: sfile, ierr
       
    character(len=30) :: fname
    character(len=10) :: cname1  
    character(len=10) :: cname2
      
    if (ismaster) then  
      sfile = 15
      
      write(cname1,'(I8)') b    
      write(cname2,'(I8)') Rstep

      fname = 'body.' // trim(adjustl(cname1))//'.'// trim(adjustl(cname2))
      open(sfile, file=fname, status='replace', iostat=ierr)
    
      ! Read header
      write(sfile,*) Rstep
      write(sfile,*) st%NMODES
    
      ! Read mesh function data        
      write(sfile,*) st%vb
      write(sfile,*) st%db
      write(sfile,*) st%wb
      write(sfile,*) st%Qb

      if (st%NMODES /= 0) then
        write(sfile,*) st%dlambda 
        write(sfile,*) st%lambda
      endif
    
      write(sfile,*) Rstep
      close(sfile)
    endif
  
  end subroutine writeBodyState   
                   
  !--------------------------------------------------------------------      
  subroutine printBodyState(st,b)
    use mpi
    use params
    implicit none
  
    type(bodyState), intent(inout) :: st
    integer,         intent(in)    :: b
    integer :: i, ii
    
    character(len=*), parameter :: fmt1 ='ES12.4'
    character(len=*), parameter :: fmt2 ='F12.06'
    character(len=30) :: fname
    character(len=10) :: cname1  
    character(len=10) :: cname2
            
    if (.not.ismaster) return

    write(*,'(15(" --"))')
    if (st%force==1) write(*,'(4x,a,1x,3ES12.4)') " Force        = ", st%Fd
    write(*,'(4x,a,1x,3ES12.4)') " Velocity     = ", st%vb
    write(*,'(4x,a,1x,3ES12.4)') " Displacement = ", st%db
    
    if (st%force==1) write(*,'(4x,a,1x,3ES12.4)') " Moment       = ", st%Fw
    write(*,'(4x,a,1x,3ES12.4)') " Omega        = ", st%wb
    write(*,'(4x,a,1x,3ES12.4)') " Rotation     = ", st%Qb(1,:)
    write(*,'(21x,3Es12.4)')                         st%Qb(2,:)
    write(*,'(21x,3Es12.4)')                         st%Qb(3,:)  
                             
    do i = 1, st%NMODES/8
      ii = min(st%NMODES, 8*i)
      if (i == 1) then            
        write(*,'(4x,a,1x,8Es12.4)') " Gen. Force   = ", st%Fa(1:ii) 
      else
        write(*,'(21x,8Es12.4)') st%Fa(8*(i-1)+1:ii) 
      endif

    enddo    
    do i = 1, st%NMODES/8
      ii = min(st%NMODES, 8*i)
      if (i == 1) then
        write(*,'(4x,a,1x,8Es12.4)') " dot(L)       = ", st%dlambda(1:ii) 
      else
        write(*,'(21x,8Es12.4)') st%dlambda(8*(i-1)+1:ii) 
      endif

    enddo
    do i = 1, st%NMODES/8
      ii = min(st%NMODES, 8*i)
      if (i == 1) then     
        write(*,'(4x,a,1x,8Es12.4)') " Lambda       = ",  st%lambda(1:ii)
      else
        write(*,'(21x,8Es12.4)')  st%lambda(8*(i-1)+1:ii)
      endif

    enddo      
    
    write(*,'(15(" --"))')
    
    if (bodyfile(b) == -1) then
      bodyfile(b) = 455 + 5*b
      call openStreamFile (bodyfile(b)+1,"vel."//trim(bodyname(b)))          
      call openStreamFile (bodyfile(b)+2,"disp."//trim(bodyname(b)))       
      call openStreamFile (bodyfile(b)+3,"force."//trim(bodyname(b))) 
    endif
    
    if (st%NMODES ==0) then      
      write(bodyfile(b)+1, '(7ES12.4)') time, st%vb, st%wb
      write(bodyfile(b)+2,'(13ES12.4)') time, st%db, st%Qb
      if (st%force==1) write(bodyfile(b)+3, '(7ES12.4)') time, st%Fd, st%Fw
    else
      write(bodyfile(b)+1,'(128ES12.4)') time, st%vb, st%wb, st%dlambda 
      write(bodyfile(b)+2,'(128ES12.4)') time, st%db, st%Qb, st%lambda  
      if (st%force==1) write(bodyfile(b)+3,'(128ES12.4)') time, st%Fd, st%Fw, st%Fa    
    endif    

     
  end subroutine printBodyState  
  !--------------------------------------------------------------------
  subroutine readBodyParams(bp, fname)

    use mpi
    implicit none

    type(bodyParams) bp
    character(len=*) :: fname  
    integer :: mfid,a,b,i,j,RBM
    real(8) :: rho, Emod
    
    ! Read mode file  
    if (ismaster) then
      ! Open mode file
      mfid = 41 
      write(*,*) "Reading body: ", fname
      open(mfid, file=fname, status='old') 

      ! Read rigid body params  
      read(mfid,*) bp%mass,rho,Emod   
    
      write(*,"(a16,1x,' = ',1x,ES12.4)")  "mass", bp%mass
      write(*,"(a16,1x,' = ',1x,ES12.4)")  "rho",  rho     
      write(*,"(a16,1x,' = ',1x,ES12.4)")  "Emod", Emod  
    
      read(mfid,*) bp%xcg       
      read(mfid,*) bp%Jij     
      bp%Jij = rho*bp%Jij 

            
      ! Read rigid body constraints    
      read(mfid,*) bp%VBC 
      read(mfid,*) bp%MBC  
      read(mfid,*) bp%EBC      

      ! Read size    
      read(mfid,*) bp%NBOUND    
      allocate(bp%Bound_ID(bp%NBOUND )) 
      read(mfid,*) bp%Bound_ID             
              
      ! Read size    
      read(mfid,*) bp%NMODES
            
      ! Remove Rigid body modes     
      RBM = min(6, bp%NMODES)
      bp%NMODES = bp%NMODES-RBM 
    
      if (bp%NMODES.eq.0) bp%EBC = 1 
     
      ! Allocate 
      allocate(bp%Mabij(bp%NMODES,bp%NMODES,3,3) ) 
      allocate(bp%Saij (bp%NMODES,3,3))
      allocate(bp%Ka   (bp%NMODES) ) 
      allocate(bp%Ma   (bp%NMODES))
      
      bp%Ma = rho
      
      ! Read global data     
      do a = 1, bp%NMODES+RBM    
        do b = 1, bp%NMODES+RBM         
          if ((a.le.RBM).or.(b.le.RBM)) then
            read(mfid,*)
          else
            read(mfid,*) ((bp%Mabij(a-RBM,b-RBM,i,j), i=1,3),j=1,3)     
          endif
        enddo        
      enddo   
      bp%Mabij = rho*bp%Mabij 
     
      do a = 1, bp%NMODES+RBM 
        if (a.le.RBM) then 
          read(mfid,*)
        else  
          read(mfid,*) ((bp%Saij(a-RBM,i,j), i=1,3),j=1,3)          
        endif
      enddo  
      bp%Saij = rho*bp%Saij 
  
      do a = 1, bp%NMODES+RBM   
        if (a.le.RBM) then 
          read(mfid,*)
        else   
          read(mfid,*) bp%Ka(a-RBM)
        endif     
      enddo   
      bp%Ka = Emod*bp%Ka
    
      write(*,*) bp%Ka
    
      ! Read mode mesh
      read(mfid,*) bp%NDOFMODE,bp%NFACEMODE
      allocate(bp%xg_mode(bp%NDOFMODE,3))
      do i = 1, bp%NDOFMODE
        read(mfid,*) (bp%xg_mode(i,j), j = 1, 3)       
      end do
      if (bp%NDOFMODE > 0) then
        write(*,*) "   Body minimum :", minval(bp%xg_mode,1) 
        write(*,*) "   Body maximum :", maxval(bp%xg_mode,1) 
      endif        
              
      allocate(bp%IEN_MODE(bp%NFACEMODE,3))
      do i = 1, bp%NFACEMODE
        read(mfid,*) (bp%IEN_MODE(i,j), j = 1, 3) ! Hardcoded triangle
      end do
    
      ! Read modes
      allocate(bp%modes  (bp%NMODES,bp%NDOFMODE,3))  
      do a = 1, bp%NMODES+RBM
        if (a.le.RBM) then 
          do i = 1, bp%NDOFMODE
            read(mfid,*) 
          enddo     
        else 
          do i = 1, bp%NDOFMODE
            read(mfid,*) (bp%modes(a-RBM,i,j), j=1,3) 
          enddo
        endif
      enddo 
      bp%modes = bp%modes
      close(mfid)
    endif    
  		 
  ! Communicate
#if USEMPI    
    call MPI_BCAST(bp%xcg,3,MPI_DOUBLE_PRECISION, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)  
    call MPI_BCAST(bp%NBOUND,1,MPI_INTEGER, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)    
    call MPI_BCAST(bp%NMODES,1,MPI_INTEGER, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)    
    call MPI_BCAST(bp%NDOFMODE,1,MPI_INTEGER, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)
  
    if (.not.ismaster) then      
      allocate(bp%Bound_ID(bp%NBOUND))                 
      allocate(bp%modes   (bp%NMODES,bp%NDOFMODE,3))  
      allocate(bp%Ma      (bp%NMODES))   
      allocate(bp%xg_mode (bp%NDOFMODE,3))   
    endif 
  
    call MPI_BCAST(bp%Bound_ID,bp%NBOUND, MPI_INTEGER, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)
                 
    call MPI_BCAST(bp%modes,bp%NMODES*bp%NDOFMODE*3, MPI_DOUBLE_PRECISION, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)
		  
    call MPI_BCAST(bp%xg_mode,bp%NDOFMODE*3, MPI_DOUBLE_PRECISION, &
                   mpi_master,MPI_COMM_WORLD,mpi_err)

#endif 

  end subroutine readBodyParams 

  !----------------------------------------------------------------------
  !  Function meshContains
  !   Checks whether globaldof is assigned to local proc
  !----------------------------------------------------------------------   
  subroutine  FindModeMap(bp, mesh)
          
    use mpi
    use meshData
    implicit none

    type(bodyParams) bp   
    type(meshStruct) mesh 
  
    integer :: b,i,ii,j,dd
    real(8) :: x(3),d(3),dn,mn,mx
    logical :: contains
    
    real(8) :: xmin(3),xmax(3)
  
    xmin= 9d9
    xmax= -9d9
  
    allocate(bp%modemap(mesh%NNODE))
    bp%modemap = -1
          
    do b = 1, mesh%nbound

      if (contains(bp%NBOUND, bp%Bound_ID, mesh%bound(b)%FACE_ID)) then     
  
        mx = 0d0
        mn = 9d9
        do i= 1, mesh%bound(b)%NNODE
          ii =  mesh%bound(b)%BNODES(i)
          x  = mesh%xg(ii,:)
          
          !do dd =1, 3
          !  xmin(dd) = min(x(dd),xmin(dd))         
          !  xmax(dd) = max(x(dd),xmax(dd)) 
          !enddo 
          
          mn = 9d9
          do j = 1, bp%NDOFMODE
            d  = x - bp%xg_mode(j,:) 
            dn = sum(d*d)
            if (dn < mn) then
              mn = dn
              bp%modemap(ii) = j
            endif
          enddo      
          !write(*,*) myid, i,ii,mn 
          mx = max(mx,mn)
        enddo
      
#if USEMPI
        call MPI_Reduce(mx, mn, 1, MPI_DOUBLE_PRECISION,  &
                        MPI_MAX, mpi_master,MPI_COMM_WORLD,mpi_err)  
#endif
    
        if (ismaster) write(*,*) b,"Maximum mismatch in matching of eigenmode = ", mn 
      endif
    enddo  
    
      !write(*,*) " Modemap minimum :", xmin
      !write(*,*) " Modemap maximum :", xmax  
       
  end subroutine FindModeMap
  
  !----------------------------------------------------------------------   
  ! Use commu to find out which nodes are shared between procs
  ! Only use dofs that are unique to master
  !----------------------------------------------------------------------
  subroutine  FindDofMap(bp, mesh, offset)
          
    use mpi
    use meshData
    implicit none

    type(bodyParams), intent(inout) :: bp   
    type(meshStruct), intent(in)    :: mesh 
    integer, intent(inout) :: offset

    real(8) :: tmp(mesh%NNODE)
    integer :: NBDOF,i,j
    tmp = 1d0
    
#if USEMPI
    call commu(tmp, 1, 'in ')
              
    call MPI_BARRIER (MPI_COMM_WORLD,mpi_err)  

    call commu(tmp, 1, 'out')
#endif 
    if (ismaster) then

    NBDOF = 2*mesh%NSD + bp%NMODES
    allocate(bp%DofMap(NBDOF))  
      
    j = 0
    do i = offset, mesh%NNODE    
      if (tmp(i) < 1.1d0) then      
        j = j + 1
        bp%DofMap(j) = i
        if (j == NBDOF ) exit
      endif
    enddo

    if (j /= NBDOF ) then
      write(*,*) "Not enough free DOFS to map elastic body"
      stop
    endif
    
    offset = i+1
    
    endif
    
  end subroutine FindDofMap
  
       
  !--------------------------------------------------------------------
  subroutine predictBodyState(new, old) 
    use params
    implicit none
    
    type(bodyState), intent(inout) :: new
    type(bodyState), intent(in)    :: old 

    new%vb = old%vb
    new%db = old%db + Delt*old%vb
  
    new%wb = old%wb
    call integrateBodyRotation(new%Qb,old%Qb,new%wb,Dtgl)        
        
    if (new%NMODES /= 0 ) then
      new%dlambda = old%dlambda    
      new%lambda  = old%lambda + Delt*old%dlambda
    endif

  end subroutine 
  !--------------------------------------------------------------------
  subroutine incrementBodyState(new, old, param, vec, fact) 
    use params
    implicit none  
    type(bodyState),  intent(inout) :: new,old
    type(bodyParams), intent(inout) :: param
    real(8),          intent(in)    :: vec(*),fact  
    integer :: i,j            
    integer, parameter :: NSD = 3  
           
    i = 0
    do j = 1, NSD
      i = i + 1
      new%vb(j) = new%vb(j) + fact*vec(param%DofMap(i))
    enddo
    
    do j = 1, NSD
      i = i + 1            
      new%wb(j) = new%wb(j) + fact*vec(param%DofMap(i))
    enddo
                
    do j = 1, param%NMODES                   
      i = i + 1          
      new%dlambda(j) = new%dlambda(j) + fact*vec(param%DofMap(i))
    enddo  
    
    ! Integrate velocity
    new%db = old%db  + Delt*0.5d0*(new%vb + old%vb)
    
    call integrateBodyRotation(new%Qb,old%Qb, 0.5d0*(new%wb+old%wb),Dtgl)
    
    do j = 1, param%NMODES   
      new%lambda(j) = old%lambda(j)  + Delt*0.5d0*(new%dlambda(j) + old%dlambda(j))
    end do
    
    ! Update inertia data
    call compInertia  (param,new) 
    
  end subroutine incrementBodyState   
    
  !--------------------------------------------------------------------  
  ! Hughes-Winget integration ==> preserves proper rotation
  !-------------------------------------------------------------------- 
  subroutine integrateBodyRotation(Qb1,Qb0,wb,Dtgl)

    implicit none
  
    real(8), intent(out) :: Qb1(3,3)
    real(8), intent(in)  :: Qb0(3,3),wb(3),Dtgl
  
    real(8) :: RHS(3,3),LHS(3,3),Ainv(3,3),det
  
    integer :: i,j,k
    real(8) :: Eijk

    RHS = Qb0*Dtgl
    LHS = 0d0   
    do i = 1, 3  
      LHS(i,i) = Dtgl
      do j = 1, 3   
        do k = 1, 3
          RHS(i,:) = RHS(i,:) + 0.5d0*Eijk(i,j,k)*wb(j)*Qb0(k,:)
          LHS(i,k) = LHS(i,k) - 0.5d0*Eijk(i,j,k)*wb(j)
        enddo
      enddo     
    enddo  
          
    call get_inverse_3x3(LHS,Ainv,det)
  
    do i = 1, 3
      do j = 1, 3
        Qb1(i,j) = sum(Ainv(i,:)*RHS(:,j))
      enddo
    enddo
      
  end subroutine integrateBodyRotation
        
  !--------------------------------------------------------------------
  !  Subroutine for getting a nodal displacement
  !   x - node coordinate
  !   i - node index
  !   f - interpolation factor (0 ==> n0,  1 ==> n1)
  !   d - displacement 
  !--------------------------------------------------------------------
  subroutine getBodyDisp(d,bst,bp,x,i)

    use mpi

    implicit none
    real(8),          intent(out) :: d(3) 
    type(bodyState),  intent(in)  :: bst
    type(bodyParams), intent(in)  :: bp
    real(8),          intent(in)  :: x(3)     
    integer,          intent(in)  :: i      
    
    integer :: a
           
    d = bst%db    
    d = d + matmul(bst%Qb, (x-bp%xcg)) - (x-bp%xcg)

    do a = 1, bst%NMODES      
      d = d + bst%lambda(a)*bp%modes(a,bp%modemap(i),:)
    enddo

  end subroutine getBodyDisp         

  !--------------------------------------------------------------------
  ! This is actually the mesh velocity due to the body displacement
  !  - bodyState bst1 at n+1 
  !  - state    st0  at n
  !--------------------------------------------------------------------
  subroutine getBodyVel(u1,bst1,bp,st0,x,i)

    use mpi
    use solution
    use params
    implicit none    
    real(8),          intent(out) :: u1(3)  
    type(bodyState),  intent(in)  :: bst1
    type(bodyParams), intent(in)  :: bp    
    type(state),      intent(in)  :: st0
    real(8),          intent(in)  :: x(3)     
    integer,          intent(in)  :: i      
            
    real(8) :: d1(3),ac1(3)
      
    call getBodyDisp(d1,bst1,bp,x,i)
      
    ac1 = obet*(dtgl**2*(d1-st0%dg(i,:))- dtgl*st0%ugm(i,:) - mbet*st0%acgm(i,:)) 
    u1  = st0%ugm(i,:) + Delt*(mgam*st0%acgm(i,:) + gami*ac1)
                    
  end subroutine getBodyVel         
        	              
  !----------------------------------------------------------------------
  !  Subroutine for computing inertia parameters 
  !  that depend on rotation/deformation of body
  !----------------------------------------------------------------------
  subroutine compInertia(pb,st)
    use mpi
    implicit none
    
    type(bodyParams), intent(in)     :: pb       
    type(bodyState),  intent(inout)  :: st
     
    real(8) :: dIda(pb%NMODES,3,3)
    real(8) :: tmp(3,3), tmp1
    integer :: a,b,i,j,k,l
    real(8) :: Eijk
 
    if (.not.ismaster) return
 
    st%Iij = pb%Jij
    return
     
    ! Compute I
    tmp = pb%Jij
    do a = 1, pb%NMODES 
      do b = 1, pb%NMODES
        tmp = tmp + st%lambda(a)*st%lambda(b)*pb%Mabij(a,b,:,:)
      enddo      
      do i = 1, 3
        do j = 1, 3
          tmp(i,j) = tmp(i,j) + st%lambda(a)*(pb%Saij(a,i,j)+pb%Saij(a,j,i))
        enddo 
      enddo         
    enddo
   
    tmp1 = 0d0 
    do i = 1, 3
      tmp1 = tmp1 + tmp(i,i)
    enddo 
        
    st%Iij = 0d0  
    do i = 1, 3
      st%Iij(i,i) = st%Iij(i,i) + tmp1
      do j = 1, 3   
        do k = 1, 3
          do l = 1, 3   
            st%Iij(i,j) = st%Iij(i,j) -  st%Qb(i,k)*st%Qb(j,l)*tmp(k,l)
          enddo 
        enddo  
      enddo 
    enddo

    ! Compute dI/da  
    do a = 1, 0 !!pb%NMODES
      tmp = 0d0
      do b = 1, pb%NMODES 
        tmp = tmp + 2d0*st%lambda(b)*pb%Mabij(a,b,:,:) 
      enddo 
      
      do i = 1, 3
        do j = 1, 3
          tmp(i,j) = tmp(i,j) + pb%Saij(a,i,j)+pb%Saij(a,j,i)
        enddo       
      enddo 
            
      tmp1 = 0d0 
      do i = 1, 3
        tmp1 = tmp1 + tmp(i,i)
      enddo 
  
      dIda(a,:,:) = 0d0  
      do i = 1, 3
        dIda(a,i,i) = dIda(a,i,i) + tmp1
        do j = 1, 3   
          do k = 1, 3
            do l = 1, 3   
              dIda(a,i,j) = dIda(a,i,j) - st%Qb(i,k)*st%Qb(j,l)*tmp(k,l)
            enddo 
          enddo  
        enddo 
      enddo    
    
    enddo    
   
    ! Compute N 
    do a = 1, pb%NMODES 
      st%Nai(a,:) = 0d0 
      tmp = 0d0        
      do i = 1, 3        
        do j = 1, 3   
          do k = 1, 3
            do l = 1, 3          
              tmp(i,j) = tmp(i,j) + st%Qb(i,k)*st%Qb(j,l)*pb%Saij(a,k,l)
            enddo                                   
          enddo        
        enddo 
      enddo
          
      do i = 1, 3  
        do j = 1, 3   
          do k = 1, 3
            st%Nai(a,i) = st%Nai(a,i) + Eijk(i,j,k)*tmp(j,k)
          enddo
        enddo     
      enddo  
      
    enddo
    
  end subroutine compInertia
  
end module elasticBody
