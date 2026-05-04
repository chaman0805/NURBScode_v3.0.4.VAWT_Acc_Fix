
      
      subroutine SparseCG_BDIAG_mp(col, row, soln, kSpace)
 (col, row,
     &     IBC, IPER, D_FLAG, P_FLAG, rhsGm,
     &     solm, lhsK22, icntu, Utol, Kspaceu,Kspaceu_mn,
     &     NNODZu, NSHLu, NSD)
          
      use mpi
      implicit none   
      
      integer :: n, i, j, k, iK, iKs, jK, lK, Kspaceu,Kspaceu_mn, 
     &           count,icntu, NNODZu, NSHLu, NSD, itask, is, lenseg
      integer :: col(NNODZu+1), row(icntu) 
      
            
      real*8 :: lhsKBdiagm(NNODZu,NSD*NSD)
      

      integer n, i, j, k, iter, kSpace
      real*8  rhstmp(mNDZ,3), prodtmp(mNDZ,3), rr, rr1, rr2, rr3, 
     &  pp, alpha, soln(mNDZ,3), tmprr, beta, pv(mNDZ,3), rr0,
     &  tauk, taukm1, zs(mNDZ,3), Binv(3,3), tmp
      
      
      allocate(lhsKBdiag(3,3,mNDZ))
      
c...  extract block-diagonal matrix vector

      
      lhsKBdiag = 0d+0
      
      do i = 1, NNODZu
        do j = colvec(i), colvec(i+1)-1
          n = rowvec(j)
          if (n.eq.i) then
            lhsKBdiagm(i,:) = LHSK22(:,j)         
          endif
        enddo
      enddo
      
      
c...  communicate the block-diagonal
      
      if (numnodes .gt. 1) then         
         call commu(lhsKBdiagm, NSD*NSD, 'in ')
         call commu(lhsKBdiagm, NSD*NSD, 'out')         
      endif      
 

c...  invert block-diagonal

      Binv = 0d+0
      
      do i = 1, mNDZ
        
        Binv(1,1) =   lhsKBdiag(2,2,i) * lhsKBdiag(3,3,i) 
     &    - lhsKBdiag(3,2,i) * lhsKBdiag(2,3,i)
        Binv(1,2) =   lhsKBdiag(3,2,i) * lhsKBdiag(1,3,i) 
     &    - lhsKBdiag(1,2,i) * lhsKBdiag(3,3,i)
        Binv(1,3) =  lhsKBdiag(1,2,i) * lhsKBdiag(2,3,i) 
     &    - lhsKBdiag(1,3,i) * lhsKBdiag(2,2,i)
        tmp          = 1d+0 / ( Binv(1,1) * lhsKBdiag(1,1,i) 
     &    + Binv(1,2) * lhsKBdiag(2,1,i)  
     &    + Binv(1,3) * lhsKBdiag(3,1,i) )
        Binv(1,1) = Binv(1,1) * tmp
        Binv(1,2) = Binv(1,2) * tmp
        Binv(1,3) = Binv(1,3) * tmp
        Binv(2,1) = (lhsKBdiag(2,3,i) * lhsKBdiag(3,1,i) 
     &    - lhsKBdiag(2,1,i) * lhsKBdiag(3,3,i)) * tmp
        Binv(2,2) = (lhsKBdiag(1,1,i) * lhsKBdiag(3,3,i) 
     &    - lhsKBdiag(3,1,i) * lhsKBdiag(1,3,i)) * tmp
        Binv(2,3) = (lhsKBdiag(2,1,i) * lhsKBdiag(1,3,i) 
     &    - lhsKBdiag(1,1,i) * lhsKBdiag(2,3,i)) * tmp
        Binv(3,1) = (lhsKBdiag(2,1,i) * lhsKBdiag(3,2,i) 
     &    - lhsKBdiag(2,2,i) * lhsKBdiag(3,1,i)) * tmp
        Binv(3,2) = (lhsKBdiag(3,1,i) * lhsKBdiag(1,2,i) 
     &    - lhsKBdiag(1,1,i) * lhsKBdiag(3,2,i)) * tmp
        Binv(3,3) = (lhsKBdiag(1,1,i) * lhsKBdiag(2,2,i) 
     &    - lhsKBdiag(1,2,i) * lhsKBdiag(2,1,i)) * tmp
        
        lhsKBdiagm(i,1) = Binv(1,1)
        lhsKBdiagm(i,2) = Binv(1,2)
        lhsKBdiagm(i,3) = Binv(1,3)
        lhsKBdiagm(i,4) = Binv(2,1)
        lhsKBdiagm(i,5) = Binv(2,2)
        lhsKBdiagm(i,6) = Binv(2,3)
        lhsKBdiagm(i,7) = Binv(3,1)
        lhsKBdiagm(i,8) = Binv(3,2)
        lhsKBdiagm(i,9) = Binv(3,3)        
      enddo


      

      

      rhstmp=0d+0
      
      rhstmp(:,1) = RHSG(:,1)
      rhstmp(:,2) = RHSG(:,2)
      rhstmp(:,3) = RHSG(:,3)
      
      soln = 0d+0
      rr = 0d+0
      rr1 = 0d+0
      rr2 = 0d+0
      rr3 = 0d+0
      
      do n = 1, NNODZu
        rr1 = rr1 + rhstmp(n,1)*rhstmp(n,1)
        rr2 = rr2 + rhstmp(n,2)*rhstmp(n,2)
        rr3 = rr3 + rhstmp(n,3)*rhstmp(n,3)
      enddo
      rr = rr1+rr2+rr3
      
      rrglob=rr
      if (numnodes .gt. 1) then         
         call MPI_ALLREDUCE (rrglob, rr, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)         
      endif      
      
      rr0 = rr
      
      do iter = 1, kSpace
        
        
c...    Premultiply residual by inverse of Bdiag
        
        do i = 1, NNODZu
          zs(i,1) = lhsKBdiagm(i,1)*rhstmpm(i,1) +
     &              lhsKBdiagm(i,2)*rhstmpm(i,2) +
     &              lhsKBdiagm(i,3)*rhstmpm(i,3)
          zs(i,2) = lhsKBdiagm(i,4)*rhstmpm(i,1) +
     &              lhsKBdiagm(i,5)*rhstmpm(i,2) +
     &              lhsKBdiagm(i,6)*rhstmpm(i,3)
          zs(i,3) = lhsKBdiagm(i,7)*rhstmpm(i,1) +
     &              lhsKBdiagm(i,8)*rhstmpm(i,2) +
     &              lhsKBdiagm(i,9)*rhstmpm(i,3)
        end do
                
        tauk = 0d+0

        do i = 1, mNDZ
          do j = 1, 3
            tauk = tauk + zs(i,j)*rhstmp(i,j)
          enddo
        enddo
        
        if(iter.eq.1) then
          beta = 0d+0
          pv = zs
        else
          beta = tauk/taukm1
          pv(:,:) = zs(:,:) + beta*pv(:,:)
        endif

        

        
c       *********************************************************
c...    Product
        if (numnodes .gt. 1) call commu(pv,NSD,'out')
               
        call SparseProdM_3D(col, row, lhsK22, pv, prodtmp,
     &                      D_FLAG, P_FLAG, NNODZu, NSHLu, icntu, NSD)
     
        if (numnodes .gt. 1) call commu(prodtmp,NSD,'in ')   
        
c       *********************************************************

			
        
c...            
        pp = 0d+0
        rr1 = 0d+0
        rr2 = 0d+0
        rr3 = 0d+0
        
        do n = 1, mNDZ
          rr1 = rr1 + pv(n,1)*prodtmp(n,1)
          rr2 = rr2 + pv(n,2)*prodtmp(n,2)
          rr3 = rr3 + pv(n,3)*prodtmp(n,3)
        enddo
        
        pp = rr1 + rr2 + rr3
        ppglob=pp
        if (numnodes .gt. 1) then         
           call MPI_ALLREDUCE (ppglob, pp, 1,
     &          MPI_DOUBLE_PRECISION,MPI_SUM, 
     &          MPI_COMM_WORLD,mpi_err)         
        endif
              
        alpha = tauk/pp
        
c....   calculate the next guess
        
        soln(:,1) = soln(:,1) + alpha*pv(:,1)
        soln(:,2) = soln(:,2) + alpha*pv(:,2)
        soln(:,3) = soln(:,3) + alpha*pv(:,3)
        
c....   calculate new res. and dot prod. 
        
        rhstmp(:,1) = rhstmp(:,1) - alpha*prodtmp(:,1)
        rhstmp(:,2) = rhstmp(:,2) - alpha*prodtmp(:,2)
        rhstmp(:,3) = rhstmp(:,3) - alpha*prodtmp(:,3)
        
        tmprr = rr
        
        rr = 0
        rr1 = 0
        rr2 = 0
        rr3 = 0
        
        do n = 1, mNDZ
          rr1 = rr1 + rhstmp(n,1)*rhstmp(n,1)
          rr2 = rr2 + rhstmp(n,2)*rhstmp(n,2)
          rr3 = rr3 + rhstmp(n,3)*rhstmp(n,3)
        enddo
        
        rr = rr1 + rr2 + rr3
        
        rrglob=rr
        if (numnodes .gt. 1) then         
           call MPI_ALLREDUCE (rrglob, rr, 1,
     &          MPI_DOUBLE_PRECISION,MPI_SUM, 
     &          MPI_COMM_WORLD,mpi_err)         
        endif  
              
c....   check for convergence 

        taukm1 = tauk
        ercheck=sqrt(rr/rr0)
        if (mod(iter,50).eq.0) then
          write(*,*) "Iteration =  ", iter
          write(*,*) "Residual reduction =  ", ercheck
        endif 
        
        ercheck=sqrt(rr/rr0)
        if ((ercheck <= epsnrm).and.(iKs.ge.Kspaceu_mn)) exit 
       
      enddo
      
c...  communicate solution      

      if (numnodes .gt. 1) then 
         call commu(soln, NSD, 'out')
      endif      
      
      deallocate(lhsKBdiag)
      
      write(*,*) "Iteration =  ", iter
      write(*,*) "Residual reduction =  ", ercheck
                
      return      
      end 
