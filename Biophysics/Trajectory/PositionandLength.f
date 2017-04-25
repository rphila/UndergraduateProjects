	PROGRAM PositionandLength
C========================================
C Rita Philavanh
C 6/16/08
C This program computes the trajectory of a certain atom 
C throughout time, and also the length of a molecule.
C Things that need to be changed for a different run is 
C marked with a !!! next to it.
C=========================================

	IMPLICIT NONE
	INTEGER :: I,J,crs,time
	LOGICAL :: bcross,left 
	REAL :: x,z,xdiff,diffzero,diffbdry,xnew(100),xorig(100) 
	REAL :: startpoint(100), endpoint, length

!	OPEN(unit=60,file='xz161-r2.dat',status='unknown') !!! file name
!	OPEN(unit=61,file='xz168-r2.dat',status='unknown') !!! file name
	OPEN(unit=65,file='xlength9and16-r2b.dat',status='unknown') !!! file name

	DO 20 J=1,2 ! Does loop for start and end of molecule
	bcross=.false.
	left=.false.
	crs=0
	  time=0
	  DO 10 I=1,100 ! 2000 timesteps/every 20 steps=100 iterations
	   OPEN(unit=50,file='pict.r2b',status='unknown') !!! .r2 file input

	   IF(I.EQ.1)THEN
		IF(J.EQ.1)READ(50,100)x,z
		IF(J.EQ.2)READ(50,200)x,z
	   END IF
	   IF(I.NE.1)THEN 
		IF(J.EQ.1)READ(50,101)x,z
		IF(J.EQ.2)READ(50,201)x,z
		xorig(I)=x ! Records original x w/o considering boundary condition 
		xdiff=abs(xorig(I)-xorig(I-1))

		! Used to check if the particle really crosses boundary or just zero:
			diffzero=abs(xorig(I))+abs(xorig(I-1))
			diffbdry=(0.5-abs(xorig(I)))+(0.5-abs(xorig(I-1)))

	      IF(diffzero.GT.diffbdry.AND.xdiff.GT.0.485)THEN !particle crosses boundary
		bcross=.true.
		IF((xorig(I)-xorig(I-1))>0.48)THEN
		  left=.true.
		  crs=crs-1
		END IF
		IF((xorig(I)-xorig(I-1))<-0.48)THEN
		  left=.false.
		  crs=crs+1
		END IF
		x=x+(FLOAT(crs))
		GOTO 123
		END IF
	      IF(bcross.AND.xdiff.LT.abs(xorig(I)-xnew(I-1)))THEN
		!particle still across boundary 
		IF(.NOT.left)x=x+(FLOAT(crs))
		IF(left)x=x-(FLOAT(crs))
		GOTO 123
	      END IF
	   END IF

123	   CONTINUE
	   xnew(I)=x ! The actual position considering boundary condition
	   x=x*23.6989 ! scale factor
	   IF(J.EQ.1)THEN
		WRITE(60,*)x,z
		startpoint(I)=x
	   END IF
 	   IF(J.EQ.2)THEN
		WRITE(61,*)x,z
		endpoint=x
		length=abs(endpoint-startpoint(I))
	!	WRITE(65,*)time,length
	PRINT *,time,length
	   END IF
	   time=time+20

C===================== atom 1 and 8 =====================
 100        FORMAT(4(/),8(/),F10.4,7X,F8.3,9591(/)) !!! position of start atom 
 101        FORMAT(19(9603(/)),3(/),8(/),F10.4,7X,F8.3,9591(/)) !!! every 20 steps
 200        FORMAT(4(/),15(/),F10.4,7X,F8.3,9584(/)) !!! position of end atom
 201        FORMAT(19(9603(/)),3(/),15(/),F10.4,7X,F8.3,9584(/)) !!! every 20 steps
C=================================================
C 100	   FORMAT(4(/),80(/),F10.4,7X,F8.3,9519(/)) !!! position of start atom 
C 101	   FORMAT(19(9603(/)),3(/),80(/),F10.4,7X,F8.3,9519(/)) !!! every 20 steps
C 200	   FORMAT(4(/),87(/),F10.4,7X,F8.3,9512(/)) !!! position of end atom
C 201	   FORMAT(19(9603(/)),3(/),87(/),F10.4,7X,F8.3,9512(/)) !!! every 20 steps
C==================================================
C 100      FORMAT(4(/),160(/),F10.4,7X,F8.3,9439(/)) !!! position of start atom 
C 101      FORMAT(19(9603(/)),3(/),160(/),F10.4,7X,F8.3,9439(/)) !!! every 20 steps
C 200      FORMAT(4(/),167(/),F10.4,7X,F8.3,9432(/)) !!! position of end atom
C 201      FORMAT(19(9603(/)),3(/),167(/),F10.4,7X,F8.3,9432(/)) !!! every 20 steps
C==================================================
C 100      FORMAT(4(/),240(/),F10.4,7X,F8.3,9359(/)) !!! position of start atom 
C 101      FORMAT(19(9603(/)),3(/),240(/),F10.4,7X,F8.3,9359(/)) !!! every 20 steps
C 200      FORMAT(4(/),247(/),F10.4,7X,F8.3,9352(/)) !!! position of end atom
C 201      FORMAT(19(9603(/)),3(/),247(/),F10.4,7X,F8.3,9352(/)) !!! every 20 steps
C==================================================
C 100      FORMAT(4(/),320(/),F10.4,7X,F8.3,9279(/)) !!! position of start atom 
C 101      FORMAT(19(9603(/)),3(/),320(/),F10.4,7X,F8.3,9279(/)) !!! every 20 steps
C 200      FORMAT(4(/),327(/),F10.4,7X,F8.3,9272(/)) !!! position of end atom
C 201      FORMAT(19(9603(/)),3(/),327(/),F10.4,7X,F8.3,9272(/)) !!! every 20 steps
C==================================================


10	  CONTINUE
	CLOSE(50)
20	CONTINUE
	CLOSE(60)
	CLOSE(61)
	CLOSE(65)
	STOP
	END
C====================================================
