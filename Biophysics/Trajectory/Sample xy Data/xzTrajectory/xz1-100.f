	PROGRAM xz1
	IMPLICIT NONE
	INTEGER :: I,T
	REAL x,y,xdiff,x2,cross,xnew(100),xorig(100)
 	OPEN(unit=60,file='xz1-100.dat',status='unknown')
	
	  DO I=1,100 !1000 timesteps total / every 10 points =100 iterations
	   OPEN(unit=50,file='pict.89aa',status='unknown')
	   IF(I.EQ.1)READ(50,100)x,y
	   IF(I.NE.1)THEN 
		READ(50,101)x,y
	  	xorig(I)=x 
		xdiff=abs(xorig(I)-xorig(I-1))
	      IF(xdiff.GT.0.485)THEN !particle crosses boundary
		x=x2(xorig,I) ! actual x position
	      ELSE IF(xdiff.LT.abs(xorig(I)-xnew(I-1)))THEN
		IF(xorig(I).LT.0)x=x+1.0D0 !particle still across boundary
		IF(xorig(I).GT.0)x=x
	      END IF
	   END IF
	   xnew(I)=x
	   x=x*23.6989 ! scale factor
	   WRITE(60,*)x,y
 	   
100	   FORMAT(4(/),F10.4,7X,F8.3,9599(/)) ! position of 1st point for every 10 timesteps
101	   FORMAT(9(9603(/)),3(/),F10.4,7X,F8.3,9599(/))
		
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
        INTEGER :: I,T
        REAL :: xarray(I), cross
	cross=0.0D0
	 DO T=2,I
            IF (xarray(T)-xarray(T-1)>0.48) THEN !left of boundary 
               cross=cross+-1.0D0
            ELSE IF (xarray(T)-xarray(T-1)<-0.48) THEN !right of boundary
               cross=cross+1.0D0
            ELSE
               cross=cross
            END IF
          END DO
	x2=xarray(I)+cross
	END

