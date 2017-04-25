       PROGRAM clusters

C================================
C Rita Philavanh
C 7/16/08
C This program reads in a struc file and and finds prob of cluster of
C certain size within specific time range
C or plots the cluster break (#/size) through time
C================================

       INTEGER :: I,J,K,L,P,N,M,T,MaxAtom,MaxSizeTot,a,b,c,d,x,f
       INTEGER :: xx, temp,z
       PARAMETER(N=500,MaxAtom=8000,f=50,z=999)
 !! N=timestep, f=0 file, z=# timestep before last step
       INTEGER :: ClustSize(10000000),time(0:N),MaxClust(0:N)
       INTEGER :: MaxSize(0:N,1:MaxAtom)
       INTEGER :: ClustMatrix(1:MaxAtom,1:MaxAtom)
       INTEGER :: OldClustMatrix(1:MaxAtom,1:MaxAtom)

       OPEN(unit=50,file='hbonds.90x_struc',status='unknown')   !! 0 file
       OPEN(unit=51,file='hbonds.40_struc',status='unknown')
!       OPEN(unit=52,file='hbonds.90m_struc',status='unknown')
!       OPEN(unit=53,file='hbonds.90n_struc',status='unknown')
!       OPEN(unit=54,file='hbonds.90p_struc',status='unknown')
!       OPEN(unit=55,file='hbonds.90r_struc',status='unknown')
!       OPEN(unit=56,file='hbonds.90s_struc',status='unknown')
!       OPEN(unit=57,file='hbonds.90t_struc',status='unknown')
!       OPEN(unit=58,file='hbonds.90u_struc',status='unknown')
!       OPEN(unit=59,file='hbonds.90v_struc',status='unknown')
!       OPEN(unit=60,file='hbonds.90w_struc',status='unknown')
!       OPEN(unit=61,file='hbonds.90x_struc',status='unknown')
!       OPEN(unit=62,file='hbonds.89y_struc',status='unknown')

!       OPEN(unit=63,file='hbonds.89z_struc',status='unknown')
!       OPEN(unit=64,file='hbonds.89aa_struc',status='unknown')
!       OPEN(unit=65,file='hbonds.89ab_struc',status='unknown')
!       OPEN(unit=66,file='hbonds.89ac_struc',status='unknown')
!       OPEN(unit=67,file='hbonds.89ad_struc',status='unknown')
!       OPEN(unit=68,file='hbonds.89ae_struc',status='unknown')
!       OPEN(unit=69,file='hbonds.89af_struc',status='unknown')
!       OPEN(unit=70,file='hbonds.89ag_struc',status='unknown')
!       OPEN(unit=71,file='hbonds.89ah_struc',status='unknown')
!       OPEN(unit=72,file='hbonds.89ai_struc',status='unknown')
!       OPEN(unit=73,file='hbonds.89aj_struc',status='unknown')
!       OPEN(unit=74,file='hbonds.89ak_struc',status='unknown')
!       OPEN(unit=75,file='hbonds.89al_struc',status='unknown')
!       OPEN(unit=76,file='hbonds.89am_struc',status='unknown')
!       OPEN(unit=77,file='hbonds.89an_struc',status='unknown')
!       OPEN(unit=78,file='hbonds.89ap_struc',status='unknown')
!       OPEN(unit=79,file='hbonds.89ar_struc',status='unknown')
!       OPEN(unit=80,file='hbonds.89as_struc',status='unknown')
!       OPEN(unit=81,file='hbonds.89at_struc',status='unknown')
!       OPEN(unit=82,file='hbonds.89au_struc',status='unknown')
!       OPEN(unit=83,file='hbonds.89av_struc',status='unknown')
!       OPEN(unit=84,file='hbonds.89aw_struc',status='unknown')
!       OPEN(unit=85,file='hbonds.89ax_struc',status='unknown')
!       OPEN(unit=86,file='hbonds.89ay_struc',status='unknown')
!       OPEN(unit=87,file='hbonds.89az_struc',status='unknown')

C---------------------------------------------------------
!       OPEN(unit=40,file='ClustSizeAvg0m-90.dat',status='unknown')
!       OPEN(unit=41,file='ClustSizeTot0m-90.dat',status='unknown')
!       OPEN(unit=42,file='ClustNum100-40.dat',status='unknown')
       OPEN(unit=43,file='ClustSame0-40.dat',status='unknown')
!       OPEN(unit=44,file='ClustGrow100-40.dat',status='unknown')
!       OPEN(unit=45,file='ClustTot10m-90.dat',status='unknown')
!       OPEN(unit=46,file='ClustTot20m-90.dat',status='unknown')


!       OPEN(unit=90,file='cluster0.dat',status='unknown')
!       OPEN(unit=91,file='cluster1.dat',status='unknown')
!       OPEN(unit=92,file='cluster2.dat',status='unknown')
!       OPEN(unit=93,file='cluster3.dat',status='unknown')
!       OPEN(unit=94,file='cluster4.dat',status='unknown')
!       OPEN(unit=95,file='cluster5.dat',status='unknown')
C---------------------------------------------------------

       PRINT *,'1. (# or Size) Cluster Break vs. Time'
       PRINT *,'2. Cluster Size vs. # of Cluster'
       READ *,M
       PRINT *,'How many timesteps in each file?'
       READ *,d
       IF(M.EQ.1)PRINT *,'1.) Compare to 0? 2.) Every 10? 3.)100?'
       READ *,x
       IF(x.EQ.2)T=10
       IF(x.EQ.3)T=100
       IF(x.EQ.1)THEN
          PRINT *,'Every ? step?'
          READ *,T       
       END IF
  
C======================================================
C Reads in struc file
C======================================================

!-------initialization-------------!
       DO J=0,N
         DO K=1,MaxAtom
            MaxSize(J,K)=0
         END DO
       END DO
       MaxSizeTot=0
       P=0
       time(0)=0

C---------- Runs through 1st file until the last time step----------------
       DO J=1,z 
         DO K=1,MaxAtom
           IF(K.EQ.1)READ(f,100)temp
           IF(K.NE.1)READ(f,101)temp
            IF(temp.EQ.0)THEN ! End of all clusters in time step
             GOTO 321 ! Go to next time step
            END IF
           DO I=1,temp
             READ(f,200)xx
           END DO
         END DO
321   CONTINUE
      END DO
C---------------------------------------------------

!-----------------------------------!
C-------------------------------------------------------------
         a=1
         b=d
         c=f+1            
       DO 18 J=0,N
         DO 14 K=1,MaxAtom
            IF(K.EQ.1.AND.J.EQ.0)READ(f,100)ClustSize(P)
            IF(K.EQ.1.AND.J.GE.a.AND.J.LE.b)READ(c,100)ClustSize(P)
!!
	IF((J.EQ.490).AND.K.EQ.1)PRINT *,J,ClustSize(P)
!!
            IF(K.NE.1.AND.J.EQ.0)READ(f,101)ClustSize(P)
            IF(K.NE.1.AND.J.GE.a.AND.J.LE.b)READ(c,101)ClustSize(P)
!!
	IF((J.EQ.490))PRINT *,J,ClustSize(P)
!!

100        FORMAT((/),9X,I3)
101        FORMAT(9X,I3)
C-----------------------------------------------------------------
!----Records size of all clusters in a time step, and finds the end----!

           IF(ClustSize(P).GT.MaxSize(J,K))MaxSize(J,K)=ClustSize(P)
           IF(MaxSize(J,K).GT.MaxSizeTot)MaxSizeTot=MaxSize(J,K)
           IF(ClustSize(P).EQ.0)THEN ! End of all clusters in time step
             MaxClust(J)=K
             time(J+1)=P+1 ! used to convert what time step P is equivalent to
             GOTO 123 ! Go to next time step
           END IF
!----------------------------------------------------------------------!

C-----------------------------------------------------------
C! Writes new time step's clusters into matrix
C-----------------------------------------------------------
           DO I=1,MaxSize(J,K)
              ClustMatrix(K,I)=0
           END DO

           DO 10 I=1,ClustSize(P)
             IF(J.EQ.0)READ(f,200)ClustMatrix(K,I)
             IF(J.GE.a.AND.J.LE.b)READ(c,200)ClustMatrix(K,I)
!!
	IF(J.EQ.0.AND.K.LE.2)PRINT *,ClustMatrix(K,I)
!!

200          FORMAT(8X,I4)
10         CONTINUE
C------------------------------------------------
           P=P+1
14      CONTINUE
123    CONTINUE

C---------------------------------------------------
C	!Call subroutines
C-----------------------------------------------
       IF((M.EQ.1).AND.(J.EQ.0.OR.MOD(J,T).EQ.0))THEN
        IF(J.GT.0)THEN
	PRINT *,J
          CALL ClusterBreak(OldClustMatrix,ClustMatrix,MaxSize,N,
     &    MaxAtom,M,MaxClust,J,x)
        END IF
C-----------------------------------------------------

	IF(x.EQ.1.AND.J.NE.0)GOTO 111 ! OldClustMatrix will remain 0
        ! Records previous time step's clusters into another matrix
         DO K=1,MaxClust(J)
           DO I=1,MaxAtom
             OldClustMatrix(K,I)=0
             OldClustMatrix(K,I)=ClustMatrix(K,I)
           END DO
         END DO
111     CONTINUE  !
        END IF
C-------------------------------------------------
C-------------------------------------------------

! If all files have same timesteps:
        IF(mod(J,d).EQ.0.AND.J.NE.0)THEN
           a=a+d
           b=b+d
           c=c+1
        END IF

! If not:
!         IF(c.LT.53.AND.mod(J,500).EQ.0.AND.J.NE.0)THEN
!           a=a+500
!           b=b+500
!           c=c+1
!        PRINT *,J,c,a,b,'mod J'
!         END IF
!         IF(c.GE.53.AND.mod(J,1000).EQ.0.AND.J.NE.0)THEN
!           IF(c.EQ.53)a=b+1
!           IF(c.NE.53)a=a+1000
!           b=b+1000
!           c=c+1
!        PRINT *,J,c,a,b,'mod J-500'
!         END IF

18    CONTINUE

      IF(M.EQ.2)CALL ClusterSize(MaxSizeTot,time,ClustSize,N)

      STOP
      END program clusters
C============================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C==============================================================
C SUBROUTINES:
C=============================================================

C===================================================================
        SUBROUTINE ClusterBreak(OldClustMatrix,ClustMatrix,MaxSize,N,
     &    MaxAtom,M,MaxClust,J,x)
C--------------------------------------------------------------
C Finds #/size of cluster breaking off at diffferent times
C----------------------------------------------------------------
        INTEGER :: N,MaxAtom,M,J,ClustRow,x,start
        INTEGER :: MaxClust(0:N),MaxSize(0:N,1:MaxAtom)
        INTEGER :: ClustMatrix(1:MaxAtom,1:MaxAtom)
        INTEGER :: OldClustMatrix(1:MaxAtom,1:MaxAtom)

        INTEGER :: SimCount,SimTot,CountTot,diffCount,K,I,K2,I2
        INTEGER :: sim(1:MaxAtom), AtomCount, CountGrow

        INTEGER :: ClustSize,MaxClustNum,CountSame,NewClustSize
        INTEGER :: SizeBreak,BreakTot,BreakCount,LargestBreak
        REAL :: BreakAvg

C----------------------------------------
      IF(x.EQ.1)start=0
      IF(x.EQ.2)start=J-10
      IF(x.EQ.3)start=J-100
      MaxClustNum=MaxClust(start)  !
      IF(MaxClust(J).GT.MaxClust(start))MaxClustNum=MaxClust(J)  !

	! Initialization
      CountTot=0
      CountSame=0
      CountGrow=0
      ClustRow=0
      BreakTot=0 
      BreakCount=0
      BreakAvg=0.0D0 
C------------------------------------------
C-----------------------------------------------------
	! Compares 2 clusters of different times
      DO 39 K=1,MaxClust(start) !
      AtomCount=0
      NewClustSize=0
      DO I = 1,MaxClustNum
            sim(I) = 0
      END DO
      LargestBreak=0
         DO 36 I=1,MaxSize(start,K) !
             IF(OldClustMatrix(K,I).NE.0)AtomCount=AtomCount+1
             DO 33 K2=1,MaxClust(J)
                 DO 30 I2= 1,MaxSize(J,K2)

                     IF((OldClustMatrix(K,I).EQ.ClustMatrix(K2,I2)).AND.
     &                  (OldClustMatrix(K,I).NE.0))THEN
                          sim(K2) = sim(K2)+1

                     END IF
30               CONTINUE
 
                 IF(sim(K2).GT.NewClustSize)THEN
                    NewClustSize=sim(K2)
                    ClustRow=K2
                 END IF
               
33           CONTINUE      
36       CONTINUE
C------------------------------------------------
C------------------------------------
	! Calculates and records the number and size of clusters breaking
        SimCount=0
        SimTot=0
        diffCount=0
        SizeBreak=0 
        DO L=1,MaxClustNum
          SimTot=SimTot + sim(L)
          IF(sim(L).NE.0)SimCount=SimCount+1
        END DO

        IF(AtomCount.NE.2)SimCount=SimCount-1
        IF(AtomCount.EQ.2.AND.SimTot.EQ.2)SimCount=SimCount-1
C----------------------------------------------------------
      !A cluster broke but disappeared:
        IF(AtomCount.EQ.2.AND.SimTot.NE.2)SimCount=1
        diffCount=AtomCount-SimTot
        IF(AtomCount.NE.2.AND.diffCount.GT.0)SimCount=SimCount+1
C-------------------------------------------------------------

	!If a cluster breaks into several pieces count as 1 break:
        IF(SimCount.GT.1)SimCount=1 

        CountTot=CountTot+SimCount
        IF(NewClustSize.EQ.AtomCount.AND.
     &           MaxSize(J,ClustRow).EQ.AtomCount)CountSame=CountSame+1
        IF(NewClustSize.EQ.AtomCount.AND.
     &           MaxSize(J,ClustRow).GT.AtomCount)CountGrow=CountGrow+1

C-------------------------------------------------
C------------------------------------------------
        SizeBreak=AtomCount-NewClustSize

        IF(SizeBreak.GT.LargestBreak)LargestBreak=SizeBreak
        BreakTot=BreakTot+SizeBreak
!        IF(SizeBreak.GT.0)BreakCount=BreakCount+1
C-------------------------------------------------
C-------------------------------------------------
39     CONTINUE

        IF(CountTot.NE.0)BreakAvg=FLOAT(BreakTot)/FLOAT(CountTot)
!        WRITE(40,*)J,BreakAvg
!        WRITE(41,*)J,BreakTot
!        WRITE(42,*)J,CountTot
        WRITE(43,*)J,CountSame   !!!
!        WRITE(44,*)J,CountGrow
!        WRITE(45,*)J,CountTot+CountSame+CountGrow
!        WRITE(46,*)J,MaxClust(J)
        RETURN
        END
C============================================================

C============================================================
       SUBROUTINE ClusterSize(MaxSizeTot,time,ClustSize,N)
C----------------------------------------------------------
C Finds # of clusters with same size
C----------------------------------------------------------
       INTEGER num,starttime,endtime,Clust,prob,K,N

      PRINT *,'How many times to run?'
      DO 28 K=1,num
        PRINT *,'Enter start and end time'
        READ *,starttime,endtime

       DO 24 Clust=2,MaxSizeTot     
         prob=0
           DO 20 P=time(starttime),time(endtime)
             IF(ClustSize(P).EQ.Clust)prob=prob+1
20         CONTINUE
           IF(K.EQ.1)WRITE(90,*)Clust,prob
           IF(K.EQ.2)WRITE(91,*)Clust,prob
           IF(K.EQ.3)WRITE(92,*)Clust,prob
           IF(K.EQ.4)WRITE(93,*)Clust,prob
           IF(K.EQ.5)WRITE(94,*)Clust,prob
           IF(K.EQ.6)WRITE(95,*)Clust,prob
24     CONTINUE
28    CONTINUE
 
      RETURN
      END
C==========================================================
