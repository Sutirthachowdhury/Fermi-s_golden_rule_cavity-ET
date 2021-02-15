!module global

!integer, parameter :: nstate = 100

!end module global
!=================================
program marcus 

!use global

implicit none


integer, parameter :: nstate = 300
integer, parameter :: nfock=100
integer, parameter :: nstep = 240
real*8, parameter:: pi=3.1415926  

real*8, allocatable :: sij(:,:)  ! these are the overlaps
real*8, allocatable:: EVECTOR(:,:),eigen(:)  ! eigenvalues and eigenvectores

real*8, allocatable :: k(:,:,:)    ! this will be the rate const for each channel
real*8, allocatable :: vda(:,:)  ! this is calculating total coupling
real*8, allocatable :: A(:,:) !this is the prefactor part
real*8, allocatable :: total_rate(:) ! this is the total rate 

real*8 :: xi,mu_eg,lam,omega_c,beta,part,mu_dd,mu_aa,R_0,diab,arg,gc

integer :: n,m,g,i,j 

real*8 :: dg

real*8 :: overlap,rho

allocate(k(nfock,nfock,nstep),vda(nfock,nfock),A(nfock,nfock),total_rate(nstep),sij(nstate,nstate))
allocate(eigen(nstate),EVECTOR(nstate,nstate))

!------ parameters(in a.u.) -------------------
diab = 0.005/27.2114
xi = 0.0/27.2114  
mu_eg = 1.0 
lam = 0.65/27.2114
omega_c = 0.04/27.2114
beta = 1052.85
mu_dd = 0.0
mu_aa = 0.0
R_0 = sqrt(2.0/(omega_c)**3)*xi*(mu_dd-mu_aa)
gc = 0.005/27.2114
arg = 2.0*(omega_c)*(gc**2)*(1.0/(beta*(omega_c**2)))
!-----------------------------------------

call getoverlap(R_0,omega_c,sij)

!-----------------  prefactor calculations-----------------------------
vda(:,:) = 0.0
A(:,:) = 0.0


 do n = 1,nfock
    do m = 1,nfock

       ! The general formula for vda  (computing in a.u.)                                                                                
         
       vda(n,m) = (diab)*(sij(n,m)) + xi*(mu_eg)*sqrt(real(m-1)+1)*(sij(n,m+1)) &
            + xi*(mu_eg)*sqrt(real(m-1))*(sij(n,m-1))
         
       !write(*,*) "n=",real(n),"m=",real(m), vda(n+1,m+1)

       ! prefactor term corresponding to each channel                                                                                    

       A(n,m) = 2.0*pi*((vda(n,m))**2 + arg*sij(n,m))*(1.0/sqrt(4.0*pi*lam*(1.0/beta)))

       !--- for wigner--------
      !A(n,m) = 2.0*pi*((vda(n,m))**2 + 1.35169E-07*sij(n,m))*(1.0/sqrt(4.0*pi*lam*(1.0/beta)))

       A(n,m) = A(n,m)*41341.37


    enddo
 enddo
!---------------------------------------------------------------------------------------------
!------------------------  rate for each channel with particular dg--------------------

 k(:,:,:) = 0.0

 do g = 1,nstep

    dg = -g*0.0005

    do n = 1,nfock
       do m = 1,nfock
        
          k(n,m,g) =  A(n,m)*exp(-(lam + dg - (real(n)*omega_c) + (real(m)*omega_c))**2/(4.0*lam*(1.0/beta)))  ! in (ps)^-1 

         !write(*,*) "n=",real(n),"m=",real(m), A(n+1,m+1)

       enddo
    enddo

   enddo

   ! total partiton function----

   part = 0.
   do n = 1,nfock

      part = part + exp(-beta*real(n-1)*omega_c)

!      write(*,*) "n=",real(n),  exp(-beta*real(n)*omega_c)

   enddo

   !--------------------------

   !---- calculating net rate-------------------------------  
   total_rate(:) = 0.
 
   do g = 1,nstep
       dg = -g*0.0005
   
      ! total_rate(g) = 0.

      do n = 1,nfock
          do m = 1,nfock

             !write(102,222) -dg*27.2114,(k(1,i+1,g),i=0,nfock-1),(k(2,j+1,g),j=0,nfock-1),(k(3,j+1,g),j=0,nfock-1)
            total_rate(g) = total_rate(g) + k(n,m,g)*(exp(-beta*real(n-1)*omega_c)/part)

            !write(*,*) "n=",real(n),"m=",real(m), k(n+1,m+1,g)*(exp(-beta*real(n)*omega_c)/part),k(n+1,m+1,g) 

         enddo
      enddo
      
      open(unit=110, file="beratan_wc0.04ev_g0.005_diab0.005.out")
      write(110,*) -dg*27.2114,total_rate(g)

   enddo
 !----------------------------------------------------------

222 format(50(e13.6,2x))

stop
end program marcus
!-----------------------------------------------
subroutine getoverlap(R_0,omega_c,sij)

!use global

implicit none


integer, parameter :: nstate = 300

integer :: i,j
real*8 :: del

real*8, intent(in) :: R_0,omega_c
real*8, intent(out) :: sij(nstate,nstate)


real*8 :: EVECTOR(nstate,nstate), eigen(nstate)
 

  do i = 1,nstate
      do j = 1,nstate

              sij(i,j) = ((i-1)+0.5)*(omega_c)*del(i,j) + 0.5*((omega_c)**2)*(R_0**2)*del(i,j) &
                        -(((omega_c**1.5)*R_0)/(sqrt(2.0)))*(sqrt(real(j-1)+1)*del(i,j+1) + sqrt(real(j-1))*del(i,j-1))   

         !write(*,*) i,j, sij(i+1,j+1)

        enddo
    enddo
           !  write(*,*) del(0,-1)
   ! write(*,*) sij(1,1)

   call DIAG(eigen,EVECTOR,sij)

   sij = EVECTOR

   ! write(*,*) "after diag=", sij(1,1)
Return
end subroutine getoverlap
!------------------------------------------------
SUBROUTINE DIAG(EVALUES,EVECT,CRV)
 
! this subroutine only for IBMSP
       
!c
!c     CRV: HERMITIAN MATRIX (INPUT)
!c     EVECT: EIGENVECTORS (OUTPUT)
!c     EVALUES: EIGENVALUES (OUTPUT)
!c
!c     THIS ROUTINE SOLVES THE EIGENVALUE PROBLEM BY CALLING
!c     THE FOLLOWING IBM ESSL FUNCTION
!c
!c     SSPEV, DSPEV, CHPEV, and ZHPEV--Eigenvalues and, Optionally,
!c     the Eigenvectors of a Real Symmetric Matrix or a Complex
!c     Hermitian Matrix
!c

      !use global

      IMPLICIT NONE

      integer, parameter :: nstate = 300
      

       CHARACTER JOBZ,UPLO
       INTEGER INFO
       INTEGER N,NAP,LDZ

       INTEGER I,J,IND
    
       real*8  AP(nstate*(nstate+1)/2),WORK(3*nstate)

       !real*8,allocatable:: EVALUES(:)
       !real*8,allocatable:: CRV(:,:),EVECT(:,:)
      
       !allocate(EVALUES(nstate))
       !allocate(CRV(nstate,nstate),EVECT(nstate,nstate))
       
       real*8 EVALUES(nstate)
       real*8 CRV(nstate,nstate),EVECT(nstate,nstate)

       N=nstate  ! electronic part of the ham

       NAP=N*(N+1)/2
       LDZ=N
      
       EVALUES=0.
       EVECT=0.

       JOBZ='V' ! calculate both eigenvalue and eogenvector
       UPLO='L' ! lower diagonal matrix

       IND=0

       DO J=1,N
          DO I=J,N
             IND=IND+1
             AP(IND)=CRV(I,J)
          END DO
       END DO

       CALL DSPEV(JOBZ,UPLO,N,AP,EVALUES,EVECT,LDZ,WORK,INFO)
 
       RETURN
     END SUBROUTINE DIAG
!---------------------------------------
!---------------------------------------
Real*8 function del(i,j)

implicit none

Integer i,j

if(i.eq.j)then
 del = 1.0
else
 del = 0.0
end if

Return
End function del
!----------------------------
