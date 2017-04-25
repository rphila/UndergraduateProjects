       PROGRAM displacement

C=====================================================
C Rita Philavanh
C 06/23/08
C This Program finds the displacement of a particle after an oscillation
C (Option M=1) and the average displacement over time (Option M=2)
C 07/07/08
C Finds the trajectory of a specified atom (Option M=3)
C======================================================


       IMPLICIT NONE
       INTEGER :: I,J,K,L,N,P,M,cross,MaxAtom,MaxMolecule
       INTEGER :: startT,endT,initial,num,prob,maxR, test
       PARAMETER(N=1200,MaxAtom=8000,MaxMolecule=1000)
       REAL :: x(0:N,1:MaxAtom),y(0:N,1:MaxAtom),z(0:N,1:MaxAtom)
       REAL :: r2(2000),xnew(0:N,1:MaxAtom),ynew(0:N,1:MaxAtom) 
       REAL :: r,startpt,endpt,step,r2Tot,r2Avg

C-------------------------------------------------------------
	!Output for check
       OPEN(unit=90,file='check.dat',status='unknown')

	!Output for r^2 vs # of atom
       OPEN(unit=60,file='displacement600-700.dat',status='unknown')
       OPEN(unit=61,file='displacement700-800.dat',status='unknown')
       OPEN(unit=62,file='displacement800-900.dat',status='unknown')
       OPEN(unit=63,file='displacement900-1000.dat',status='unknown')
       OPEN(unit=64,file='displacement1000-1100.dat',status='unknown')
       OPEN(unit=65,file='displacement1100-1200.dat',status='unknown')

	!Output for <r^2> over time
       OPEN(unit=70,file='displacementAvg.dat',status='unknown')


	!Output for trajectories
       OPEN(unit=80,file='trajectory0.dat',status='unknown')       
       OPEN(unit=81,file='trajectory1.dat',status='unknown')
       OPEN(unit=82,file='trajectory2.dat',status='unknown')
       OPEN(unit=83,file='trajectory3.dat',status='unknown')
       OPEN(unit=84,file='trajectory4.dat',status='unknown')
       OPEN(unit=85,file='trajectory5.dat',status='unknown')
C--------------------------------------------------------------

C=========================================================================
C These DO loops read in the trajectories of each atom into a matrix
C for every time step, and then accounts for boundary conditions
C==========================================================================
       PRINT *,'Check the values of which atom? If none, enter 0'
       READ *,P
         
       DO 15 J=1,N !0,N 
       OPEN(unit=50,file='pict.oin',status='unknown')
       OPEN(unit=51,file='pict.o305',status='unknown')
       OPEN(unit=52,file='pict.305b',status='unknown')
       OPEN(unit=53,file='pict.305c',status='unknown')
       OPEN(unit=54,file='pict.303d',status='unknown')
       OPEN(unit=55,file='pict.303e',status='unknown')

         DO 10 I=1,MaxAtom

         IF(I.EQ.1)THEN
	    IF(J.EQ.0)READ(50,100)x(J,I),y(J,I),z(J,I)
	    IF(J.EQ.1)READ(51,100)x(J,I),y(J,I),z(J,I)
	    IF(J.GT.1.AND.J.LE.500)READ(51,200)x(J,I),y(J,I),z(J,I)
	    IF(J.EQ.501)READ(52,100)x(J,I),y(J,I),z(J,I)
	    IF(J.GT.501.AND.J.LE.1000)READ(52,200)x(J,I),y(J,I),z(J,I)
	    IF(J.EQ.1001)READ(53,100)x(J,I),y(J,I),z(J,I)
	    IF(J.GT.1001.AND.J.LE.1500)READ(53,200)x(J,I),y(J,I),z(J,I)
	    IF(J.EQ.1501)READ(54,100)x(J,I),y(J,I),z(J,I)
	    IF(J.GT.1501.AND.J.LE.2000)READ(54,200)x(J,I),y(J,I),z(J,I)
	    IF(J.EQ.2001)READ(55,100)x(J,I),y(J,I),z(J,I)
	    IF(J.GT.2001.AND.J.LE.2500)READ(55,200)x(J,I),y(J,I),z(J,I)
          END IF
          IF(I.NE.1)THEN
	    IF(J.EQ.0)READ(50,333)x(J,I),y(J,I),z(J,I)
	    IF(J.GE.1.AND.J.LE.500)READ(51,333)x(J,I),y(J,I),z(J,I)
	    IF(J.GE.501.AND.J.LE.1000)READ(52,333)x(J,I),y(J,I),z(J,I)
	    IF(J.GE.1001.AND.J.LE.1500)READ(53,333)x(J,I),y(J,I),z(J,I)
	    IF(J.GE.1501.AND.J.LE.2000)READ(54,333)x(J,I),y(J,I),z(J,I)
	    IF(J.GE.2001.AND.J.LE.2500)READ(55,333)x(J,I),y(J,I),z(J,I)
          END IF

100    FORMAT(4(/),F10.4,F7.4,F8.3)
200    FORMAT(1600(/),3(/),F10.4,F7.4,F8.3)
333    FORMAT(F10.4,F7.4,F8.3)

C----------------------------------------------------------------           
C----------------------------------------------------------------
         ! checks if particle crosses x boundary
               cross=0
               L=1
               DO L=1,J
                 IF(x(L,I)-x(L-1,I).GT.0.48)cross=cross-1  ! left of boundary 
                 IF(x(L,I)-x(L-1,I).LT.-0.48)cross=cross+1 ! right of boundary
               END DO
               xnew(J,I)=(x(J,I)*23.6989)+((FLOAT(cross))*23.6989)
C----------------------------------------------------------------
C-----------------------------------------------------------------
         ! checks if particle crosses y boundary
               cross=0
               L=1
               DO L=1,J
                 IF(y(L,I)-y(L-1,I).GT.0.48)cross=cross-1  ! left of boundary 
                 IF(y(L,I)-y(L-1,I).LT.-0.48)cross=cross+1 ! right of boundary
               END DO
               ynew(J,I)=(y(J,I)*20.5239)+((FLOAT(cross))*20.5239)
C-----------------------------------------------------------
C-----------------------------------------------------------

C-----------------------------------------
       ! Writes the input trajectories of specified atom to check values
         IF(P.NE.0.AND.I.EQ.P)THEN
          WRITE(90,*)I,J,xnew(J,I),ynew(J,I),z(J,I)  
         END IF
C----------------------------------------

10       CONTINUE
15     CONTINUE
       CLOSE(50)
       CLOSE(55)
C==========================================================================


C======================================================================
C The rest of this code either calculates r^2 or the trajectories, depending
C on what the user chooses
C=====================================================================
      PRINT *,'1.) r^2 vs # of atom, 2.) time vs <r^2> 3.)trajectory'
      READ *,M

C==========================
C 3.) Trajectory
C---------------------------
      IF(M.EQ.3)THEN
      PRINT *,'How many atom (<7)?'
      READ *,P
      PRINT *,'start time, end time, and time step?'
      READ *,startT,endT,step
      DO K=1,P
        PRINT *,'Which atom?'
        READ *,I
!        WRITE(80,*)''
!        WRITE(80,*)I
          DO J=startT,endT,step
!            IF(J.EQ.1000.OR.J.EQ.1001.OR.J.EQ.1100)THEN
              IF(K.EQ.1)WRITE(80,*)xnew(J,I),ynew(J,I),z(J,I)
!          IF(J.NE.1000)THEN
!           r=(((xnew(J,I)-xnew(1000,I))**2)+((ynew(J,I)
!     &        -ynew(1000,I))**2) + ((z(J,I)-z(1000,I))**2))
!            WRITE(80,*)'                         ',J,'-1000: r^2= ',r
!         END IF
             IF(K.EQ.2)WRITE(81,*)xnew(J,I),ynew(J,I),z(J,I)
             IF(K.EQ.3)WRITE(82,*)xnew(J,I),ynew(J,I),z(J,I)
             IF(K.EQ.4)WRITE(83,*)xnew(J,I),ynew(J,I),z(J,I)
             IF(K.EQ.5)WRITE(84,*)xnew(J,I),ynew(J,I),z(J,I)
             IF(K.EQ.6)WRITE(85,*)xnew(J,I),ynew(J,I),z(J,I)
          END DO 
      END DO   
      GOTO 345
      END IF
C--------------------------------
C=================================

C======================================================================
C 1.) probability and 2.) time  : vs.     r^2
C These next loops calculate the displacement for inputed time step(s)
C--------------------------------------------------

       IF(M.EQ.1)THEN
         PRINT *,'How many times to runs (<7)?'
         READ *,num
         initial=1 
         step=1      
       END IF
       IF(M.EQ.2)THEN
         startT=1 !0
         endT=100 !
         initial=1 !0
         num=N
         step=100 !
       END IF
         DO 40 J=initial,num,step
           
           IF(M.EQ.1)THEN
              PRINT *,'Enter (start,end) time to calculate r^2:'
              READ *,startT,endT
              maxR=0.0D0
           END IF          
C---------------------------------------
           I=1
           P=1
           DO 20 K=1,1000

           r2(P)=(((xnew(endT,I)-xnew(startT,I))**2)+((ynew(endT,I)
     &        -ynew(startT,I))**2) + ((z(endT,I)-z(startT,I))**2))

             IF(M.EQ.1.AND.r2(P).GT.maxR)maxR=r2(P)
             IF(M.EQ.2)r2Tot=r2Tot+r2(P)

           r2(P+1)=(((xnew(endT,I+7)-xnew(startT,I+7))**2)+((ynew(endT,
     &     I+7)-ynew(startT,I+7))**2)+((z(endT,I+7)-z(startT,I+7))**2))

             IF(M.EQ.1.AND.r2(P+1).GT.maxR)maxR=r2(P+1)
             IF(M.EQ.2)r2Tot=r2Tot+r2(P)    

             I=I+8
             P=P+2
20         CONTINUE
           IF(M.EQ.1)GOTO 123 
C--------------------------------
!	PRINT *,endT,r2Tot,r2Avg
C----------------------------------
	! 2.) time vs. <r^2>
   !     IF(mod(endT,100).EQ.0)THEN
          r2Avg=r2Tot/(2000.0D0)
          WRITE(70,*)endT,r2Avg
          r2Avg=0.0D0
          r2Tot=0.0D0
    !    END IF
          startT=startT+100 !
          endT=endT+100 !
          GOTO 234
C----------------------------------

C-----------------------------------------------------
	! 1.)r^2 vs. probability for each displacement
123     CONTINUE
         r=.1
         maxR=maxR/r
         startpt=0
         endpt=.1
         DO 35 K=0,maxR

           prob=0
           DO 30 P=1,2000
             IF(r2(P).GE.startpt.AND.r2(P).LT.endpt)prob=prob+1
C-------! A Check
!         IF(r.GT.6.AND.prob.GT.0.AND.prob.NE.test)THEN
!          PRINT *,K,prob,P,r2(P)
!        END IF
!        test=prob
C--------!!

30         CONTINUE        
           IF(J.EQ.1)WRITE(60,*)r,prob
           IF(J.EQ.2)WRITE(61,*)r,prob
           IF(J.EQ.3)WRITE(62,*)r,prob
           IF(J.EQ.4)WRITE(63,*)r,prob
           IF(J.EQ.5)WRITE(64,*)r,prob
           IF(J.EQ.6)WRITE(65,*)r,prob
           
           startpt=startpt+.1
           endpt=endpt+.1
           r=r+.1
35       CONTINUE
C---------------------------------------------

234    CONTINUE  

40     CONTINUE
C-----------------------------------------------
C=============================================================
345    CONTINUE
       CLOSE(60)
       CLOSE(61)
       CLOSE(62)
       CLOSE(63)
       CLOSE(64)
       CLOSE(65)
       CLOSE(70)
       CLOSE(80)
       CLOSE(81)       
       CLOSE(82)
       CLOSE(83)
       CLOSE(84)
       CLOSE(85)
       CLOSE(90)
       STOP
       END PROGRAM displacement
C=============================================================
C=============================================================
