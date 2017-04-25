	PROGRAM xz8
	IMPLICIT NONE
	INTEGER :: I,T,cross
	REAL x,y,xdiff,x2,xnew(100),xorig(100)
 	OPEN(unit=60,file='xz1-100.dat',status='unknown')
	
	  DO I=1,100 ! 100 iterations*every 10 steps=1000 time steps total
	   OPEN(unit=50,file='pict.89aa',status='unknown')
	   IF(I.EQ.1)READ(50,100)x,y
	   IF(I.NE.1)THEN 
		READ(50,101)x,y
	  	xarray(I)=x 
		xdiff=abs(xorig(I)-xorig(I-1))
	      IF(xdiff.GT.0.485)THEN !particle crosses boundary
		x=x2(xorig,I)
	      ELSE IF(xdiff.LT.abs(xorig(I)-xnew(I-1)))THEN
		IF(xorig(I).LT.0)x=x+1 !particle still across boundary
		IF(xorig(I).GT.0)x=x
	      END IF
	   END IF
	   xnew(I)=x
	   x=x*23.6989 ! scale factor
	   WRITE(60,*)x,y
 	   
100	   FORMAT(4(/),7(/),F10.4,7X,F8.3,9592(/)) ! position of 8th atom for every 10 iterations
101	   FORMAT(9(9603(/)),3(/),7(/),F10.4,7X,F8.3,9592(/))
		
	  END DO
	CLOSE(50)
	CLOSE(60)
	STOP
	END

C==================================================
C==================================================
! Finds where particle is when it crosses boundary
	REAL FUNCTION x2(xarray,I)
        IMPLICIT NONE
        INTEGER :: I,T,cross
        REAL :: xarray(I) 
	cross=0 
	 DO T=2,I
            IF (xarray(T)-xarray(T-1)>0.48) THEN !particle moved left 
               cross=cross+-1
            ELSE IF (xarray(T)-xarray(T-1)<-0.48) THEN !particle moved right
               cross=cross+1
            ELSE
               cross=cross
            END IF
          END DO
	x2=xarray(I)+cross
	END

