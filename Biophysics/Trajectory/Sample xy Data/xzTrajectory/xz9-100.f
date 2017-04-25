	PROGRAM xz9
	IMPLICIT NONE
	INTEGER :: I,T,cross
	REAL x,y,xdiff,x2,x3(100),xarray(100)
 	OPEN(unit=60,file='xz9-100.dat',status='unknown')
	
	  DO I=1,100 ! 100 iterations*every 10 steps=1000 time steps total
	   OPEN(unit=50,file='pict.89aa',status='unknown')
	   IF(I.EQ.1)READ(50,100)x,y
	   IF(I.NE.1)THEN 
		READ(50,101)x,y
	  	xarray(I)=x 
		xdiff=xarray(I)-xarray(I-1)
	      IF(xdiff.GT.0.48.OR.xdiff.LT.-0.48)THEN !particle crosses boundary
		x=x2(xarray,I)
	      ELSE IF(xdiff.LT.abs(xarray(I)-x3(I-1)))THEN
		IF(xarray(I).LT.0)x=x+1 !particle still across boundary
		IF(xarray(I).GT.0)x=x
	      END IF
	   END IF
	   x3(I)=x
	   x=x*23.6989 ! scale factor
	   WRITE(60,*)x,y
 	   
100	   FORMAT(4(/),8(/),F10.4,7X,F8.3,9591(/)) ! position of 9th atom for every 10 iterations
101	   FORMAT(9(9603(/)),3(/),8(/),F10.4,7X,F8.3,9591(/))
		
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

