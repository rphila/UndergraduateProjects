	PROGRAM xz1
	INTEGER I
	REAL x,y
	OPEN(unit=60,file='xz1.dat',status='unknown')

	  DO I=1,1000
	   OPEN(unit=50, file='pict.89aa',status='unknown')
	   IF(I.EQ.1)READ(50,100)x,y
	   IF(I.NE.1)READ(50,101)x,y
	  
	   call bdrycros(x,xcros, )
	   x=x*23.6989
	   WRITE(60,*)x,y
	
100	   FORMAT(4(/),F10.4,7X,F8.3,9599(/))
101	   FORMAT(3(/),F10.4,7X,F8.3,9599(/))
		
	  END DO
	CLOSE(50)
	CLOSE(60)
	STOP
	END
C==================================================
C Subroutine taken from Joris and edited to fit my code
C==================================================
	SUBROUTINE bdrycros(X,Y,Xcros,Ycros,timesteps,pictamount)
!I want for all particles arrays that show at all timesteps how many times the boundary was crossed in both x and y -directions
        IMPLICIT NONE
        REAL, DIMENSION(1:40,1:100,1:pictamount*timesteps), INTENT(IN) :: X,Y
        REAL, DIMENSION(1:40,1:100,1:pictamount*timesteps), INTENT(INOUT) :: Xcros,Ycros
        INTEGER ::  I,J,T
        INTEGER, INTENT(IN) ::timesteps,pictamount
	  DO T=2,timesteps*pictamount
            DO I=1,40
                DO J=1,100
                        IF (X(I,J,T)-X(I,J,T-1)>0.5) THEN !particle moved left 
                                Xcros(I,J,T)=Xcros(I,J,T-1)-1
                        ELSE IF (X(I,J,T)-X(I,J,T-1)<-0.5) THEN !particle moved right
                                Xcros(I,J,T)=Xcros(I,J,T-1)+1
                        ELSE
                                Xcros(I,J,T)=Xcros(I,J,T-1)
                        END IF
                        
                        IF (Y(I,J,T)-Y(I,J,T-1)>0.5) THEN
                                Ycros(I,J,T)=Ycros(I,J,T-1)-1
                        ELSE IF (Y(I,J,T)-Y(I,J,T-1)<-0.5) THEN
                                Ycros(I,J,T)=Ycros(I,J,T-1)+1
                        ELSE
                                Ycros(I,J,T)=Ycros(I,J,T-1)
                        END IF
                END DO
            END DO
	  END DO
	END SUBROUTINE bdrycros

