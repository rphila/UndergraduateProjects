	PROGRAM xz9
	IMPLICIT NONE
	INTEGER :: I
	REAL x,y
 	OPEN(unit=60,file='xz9.dat',status='unknown')
	
	  DO I=1,100 ! 100 iterations*every 10 steps=1000 time steps total
	   OPEN(unit=50,file='pict.89aa',status='unknown')
	   IF(I.EQ.1)READ(50,100)x,y
	   IF(I.NE.1)READ(50,101)x,y
	  	 
	   x=x*23.6989 ! scale factor
	   WRITE(60,*)x,y
 	   
100	   FORMAT(4(/),8(/),F10.4,7X,F8.3,9591(/)) ! position of 9th atom for every 10 iterations
101	   FORMAT(9(9603(/)),3(/),8(/),F10.4,7X,F8.3,9591(/))
		
	  END DO
	CLOSE(50)
	CLOSE(60)
	STOP
	END

