      subroutine SingleRigidBody3D(massb,Ibhat,Rn,Delt,Fb,Mb,FbJ,MbJ,
     &     vbn,dbn,wbn,VBC,MBC,Nnewt,vbnp1,dbnp1,wbnp1,Rnp1,
     &     normFb,normMb,relfac)

      implicit none

!
! Input (CM - center of mass):
!
! massb - mass of the rigid body 
! Ibhat(3,3) - inertia tensor in the reference configuration
! 
! Ibhat = \int_{\Omega_0}
!                       (\rho (X-X_0) \cdot (X-X_0)) \B{I}
!                       -\rho (X-X_0) \otimes (X-X_0)
!               d \Omega_0
!
! Rn(3,3) - rotation matrix at time step n
! Delt - time step size
! Fb(3) - vector of global forces at the n+1/2 config.
! Mb(3) - vector of global moments (about CM) at the n+1/2 config.
! FbJ(3,3) - jacobian of global forces at the n+1/2 config.
! MbJ(3,3) - jacobian of global moments (about CM) at the n+1/2 config.
! vbn(3) - CM translational velocity vector at time step n
! dbn(3) - CM displacement vector at time step n
! wbn(3) - rotational velocity vector at time step n
! VBC(3) - array of transl. vel. constraints (1 - constr, 0 - not)
! MBC(3) - array of rot. vel. constraints (1 - constr, 0 - not)
! Nnewt - number of "Newton" iterations for the rotation solve (use 10)
! relfac - relaxation factor  
!
! Output:
! 
! vbnp1(3) - CM translational velocity vector at time step n+1
! dbnp1(3) - CM displacement vector at time step n+1
! wbnp1(3) - rotational velocity vector at time step n+1
! Rnp1(3,3) - rotation matrix at time step n+1
! normFb - Norm of initial linear  momentum residual
! normMb - Norm of initial angular momentum residual

      integer i,j,k,inewt, Nnewt
      integer VBC(3), MBC(3)
      real*8 massb, Ibhat(3,3), Rn(3,3), Delt, Fb(3), Mb(3),   
     &     FbJ(3,3), MbJ(3,3),vbn(3),dbn(3), wbn(3)
      real*8 Ibn(3,3), Ibnp1(3,3), Rnp1(3,3), vbnp1(3), dbnp1(3), 
     &     wbnp1(3), relfac(2),
     &     Wmatn(3,3), Wmatnp1(3,3), dvnp1(3), ddnp1(3), dwnp1(3), 
     &     dRnp1(3,3),
     &     Resv(3), Resw(3), ResR(3,3), Ainv(3,3),  Wmatnp5(3,3),
     &     Rnp5(3,3), RmatLHS(3,3),normFb,normMb

! Get the CM velocity at n+1

      Resv(:) = Fb(:) - massb/Delt*(vbnp1-vbn)
      
! Set constraints on CM velocity
      do i = 1, 3
         if(VBC(i).eq.1) Resv(i) = 0d0
      enddo
      
! Compute linear momentum residual norm
      normFb = sqrt(sum(Resv*Resv))
      
! Solve for CM velocity increment and update
      dvnp1 = Delt/massb*Resv
      vbnp1 = vbnp1+ relfac(1)*dvnp1

! Use velocity solution to update displacement solution
      Resv(:) = 0.5d0*(vbn+vbnp1) - (dbnp1-dbn)/Delt

! Set constraints on CM displacement
      do i = 1, 3
         if(VBC(i).eq.1) Resv(i) = 0d0
      enddo

! Solve for CM displacement increment and update
      ddnp1(:) = Delt*Resv(:)
      dbnp1 = dbnp1 + ddnp1

      
          
!------------------------------
! Angular momentum balance
!------------------------------

! Initialization
      Wmatn = 0d0
      Wmatn(1,2) = -wbn(3)
      Wmatn(1,3) =  wbn(2)
      Wmatn(2,1) =  wbn(3)
      Wmatn(2,3) = -wbn(1)
      Wmatn(3,1) = -wbn(2)
      Wmatn(3,2) =  wbn(1)
 
! Use R_{n} to compute I_{n}
      Ibn = 0d0
      do i = 1,3
         do j = 1,3

            do k = 1,3
                  Ibn(i,j) = Ibn(i,j) + 
     &                 Rn(i,k)*sum(Ibhat(k,:)*Rn(j,:))
            enddo

         enddo
      enddo

! ITERATE to get angular velocity and rotation matrix at n+1
      do inewt = 1, Nnewt

         ! Use R_{n+1} to compute I_{n+1}
         Ibnp1 = 0d0
         do i = 1,3
            do j = 1,3
               
               do k = 1, 3
                  Ibnp1(i,j) = Ibnp1(i,j) + 
     &                 Rnp1(i,k)*sum(Ibhat(k,:)*Rnp1(j,:))
               enddo
               
            enddo
         enddo

         do i = 1, 3
         Resw(i) = Mb(i) - 
     &        (sum(Ibnp1(i,:)*wbnp1(:)) - sum(Ibn(i,:)*wbn(:)))/Delt
         enddo

! Set constraints on rotations
         do i = 1, 3
           if (MBC(i).eq.1) then
              Resw(i) = 0d0
              Ibnp1(i,:) = 0d0
              Ibnp1(:,i) = 0d0
              Ibnp1(i,i) = 1d0
           endif
         enddo
	 
	 ! Compute angular momentum residual norm
         if (inewt.eq.1) normMb = sqrt(sum(Resw*Resw))         

! Solve for rotational velocity increment and update

         call Inv3by3(Ibnp1/Delt,Ainv)
         
         do i = 1, 3
            dwnp1(i) = sum(Ainv(i,:)*Resw(:))
         enddo

         wbnp1 = wbnp1+ relfac(2)*dwnp1

! Solve for the increment of the rotation matrix

! Obtain the spin tensor and rotation matrix at 1/2 step

         Wmatnp1 = 0d0
         Wmatnp1(1,2) = -wbnp1(3)
         Wmatnp1(1,3) =  wbnp1(2)
         Wmatnp1(2,1) =  wbnp1(3)
         Wmatnp1(2,3) = -wbnp1(1)
         Wmatnp1(3,1) = -wbnp1(2)
         Wmatnp1(3,2) =  wbnp1(1)

         Wmatnp5 = 0.5d0*(Wmatn+Wmatnp1)
         Rnp5 = 0.5d0*(Rn+Rnp1)

         do i = 1, 3
            do j = 1, 3
               ResR(i,j) = sum(Wmatnp5(i,:)*Rnp5(:,j))
     &                   - (Rnp1(i,j)-Rn(i,j))/Delt
            enddo
         enddo

         RmatLHS = 0d0
         do i = 1,3
            RmatLHS(i,i) = 1d0/Delt
         enddo
         RmatLHS = RmatLHS - 0.5d0*Wmatnp5


! Solve for rotation matrix increment         
         call Inv3by3(RmatLHS,Ainv)
         
         do i = 1, 3
            do j = 1, 3
               dRnp1(i,j) = sum(Ainv(i,:)*ResR(:,j))
            enddo
         enddo

         Rnp1 = Rnp1 + dRnp1

      enddo ! end "Newton" iteration loop

      return
      end
 
!========================================      
      subroutine integrateRotation(Rnp1,Rnp0,wbn0,Delt)

!  Integrate rotation matrix using Hughes-Winget formulation

      implicit none

      real*8 Rnp1(3,3),Rnp0(3,3),wbn0(3),Delt      
      real*8 Wmat(3,3), ResR(3,3),RmatLHS(3,3), Ainv(3,3)
      integer i,j
      
      Wmat = 0d0
      Wmat(1,2) = -wbn0(3)
      Wmat(1,3) =  wbn0(2)
      Wmat(2,1) =  wbn0(3)
      Wmat(2,3) = -wbn0(1)
      Wmat(3,1) = -wbn0(2)
      Wmat(3,2) =  wbn0(1)

      do i = 1, 3
         do j = 1, 3
            ResR(i,j) = 0.5d0*sum(Wmat(i,:)*Rnp0(:,j))+Rnp0(i,j)/Delt
         enddo
      enddo

      RmatLHS = - 0.5d0*Wmat
      do i = 1,3
         RmatLHS(i,i) = RmatLHS(i,i)+ 1d0/Delt
      enddo
        
      call Inv3by3(RmatLHS,Ainv)
      
      do i = 1, 3
         do j = 1, 3
            Rnp1(i,j) = sum(Ainv(i,:)*ResR(:,j))
         enddo
      enddo

      return
      end

!===========================================
      subroutine Inv3by3(Amat,Ainv)

      implicit none

      real*8 Amat(3,3), Ainv(3,3), tmp

      Ainv = 0d+0
            
      Ainv(1,1) =   Amat(2,2) * Amat(3,3) 
     &  - Amat(3,2) * Amat(2,3)
      Ainv(1,2) =   Amat(3,2) * Amat(1,3) 
     &  - Amat(1,2) * Amat(3,3)
 
      Ainv(1,3) =  Amat(1,2) * Amat(2,3) 
     &  - Amat(1,3) * Amat(2,2)
      tmp          = 1d+0 / ( Ainv(1,1) * Amat(1,1) 
     &  + Ainv(1,2) * Amat(2,1)  
     &  + Ainv(1,3) * Amat(3,1) )
      Ainv(1,1) = Ainv(1,1) * tmp
      Ainv(1,2) = Ainv(1,2) * tmp
      Ainv(1,3) = Ainv(1,3) * tmp
      Ainv(2,1) = (Amat(2,3) * Amat(3,1) 
     &  - Amat(2,1) * Amat(3,3)) * tmp
      Ainv(2,2) = (Amat(1,1) * Amat(3,3) 
     &  - Amat(3,1) * Amat(1,3)) * tmp
      Ainv(2,3) = (Amat(2,1) * Amat(1,3) 
     &  - Amat(1,1) * Amat(2,3)) * tmp
      Ainv(3,1) = (Amat(2,1) * Amat(3,2) 
     &  - Amat(2,2) * Amat(3,1)) * tmp
      Ainv(3,2) = (Amat(3,1) * Amat(1,2) 
     &  - Amat(1,1) * Amat(3,2)) * tmp
      Ainv(3,3) = (Amat(1,1) * Amat(2,2) 
     &  - Amat(1,2) * Amat(2,1)) * tmp

      return
      end 
