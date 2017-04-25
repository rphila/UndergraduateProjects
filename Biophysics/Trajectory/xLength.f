       PROGRAM xLength

C=====================================================
C Rita Philavanh
C 06/17/08
C This Program finds the average length in the x direction
C======================================================

       IMPLICIT NONE
!       LOGICAL :: bcross,left
       INTEGER :: I,J,K,L,N,T,cross,MaxAtom,MaxMolecule
       REAL :: xnew,xdiff,length,lengthTot,avgLength !,diffzero,diffbdry
       PARAMETER(N=100,MaxAtom=8000,MaxMolecule=1000)
       REAL :: x(1:N,1:MaxAtom),z(1:N,1:MaxAtom),xorig(N)
       REAL :: startpoint(N),endpoint(N)
   
       OPEN(unit=60,file='avgxLengthall-r2b.dat',status='unknown')

C==========================================================================
C These first 2 DO loops read in the trajectories of each atom into a matrix
C==========================================================================
       DO 15 J=1,N !N=100 iterations=2000 timestep/every 20 steps
       OPEN(unit=50,file='pict.r2b',status='unknown')

         DO 10 I=1,MaxAtom

	  IF(I.EQ.1)THEN
	    IF(J.EQ.1)READ(50,100)x(J,I),z(J,I)
	    IF(J.NE.1)READ(50,200)x(J,I),z(J,I)
	  END IF
	    IF(I.NE.1)READ(50,333)x(J,I),z(J,I)

100    FORMAT(4(/),F10.4,7X,F8.3)
200    FORMAT(1603(/),F10.4,7X,F8.3)
333    FORMAT(F10.4,7X,F8.3)

10       CONTINUE
15     CONTINUE
       CLOSE(50)

C============================================================
C These next Do loops takes into account the periodic boundary conditions 
C It also finds the x length of each molecule and averages them
C============================================================
       T=0
       DO 25 J=1,N
       lengthTot=0
       avgLength=0
       I=1
         DO 20 K=1,MaxMolecule
           length=0
           endpoint(J)=x(J,I+7)
           startpoint(J)=x(J,I)
           length=abs(xnew(endpoint,J)-xnew(startpoint,J))*23.6989
           lengthTot=lengthTot+length
           I=I+8
20       CONTINUE

	avgLength=lengthTot/FLOAT(MaxMolecule)
	WRITE(60,*)T,avgLength
	T=T+20

25      CONTINUE

       CLOSE(60)
       STOP
       END PROGRAM xLength
C============================================================
C============================================================
       REAL FUNCTION xnew(xorig,J)
       IMPLICIT NONE
       INTEGER :: I,J,cross
       REAL :: xorig(J)

       cross=0
       DO I=2,J
         IF(xorig(I)-xorig(I-1).GT.0.48)cross=cross-1  ! left of boundary 
         IF(xorig(I)-xorig(I-1).LT.-0.48)cross=cross+1 ! right of boundary
       END DO
       xnew=xorig(J)+FLOAT(cross)
       END FUNCTION xnew
C==============================================================

