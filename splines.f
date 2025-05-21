c
c These routines are taken from Numerical Recipes
c
      SUBROUTINE SPLINE (R,X,Y,S,N,MXN,V)
c	this is the routine which actually uses cubic splines for
c	interpolation. SUBR. ONESPL prepares the data which are used
c	by this routine to do its work.
c	if the point at which interpolation is desired lies outside 
c	the range covered by X, then a linear extrapolation is
c	performed.
c 
c history:
c 95/03/20 (mjw): allow one to only use part of a passed array
c                   by letting N be the number of valid indices and
c		    MXN the dimensions of the arrays
c

      implicit none

      integer N,MXN,NHI,NLO,MS,IR,IL
      real*8 X(MXN),Y(MXN),S(MXN),R,V,H,RS1,RS2,SLOPE

      IF (R.LT.X(1).OR.R.GT.X(N)) GO TO 80
      NHI=N
      NLO=1
      MS=N/2
    5 CONTINUE
      IF (R-X(MS)) 10,20,30
   20 CONTINUE
      V=Y(MS)
      RETURN
   10 CONTINUE
      NHI=MS
   15 CONTINUE
      IF (NHI.EQ.(NLO+1)) GO TO 50
      MS=(NHI+NLO)/2
      GO TO 5
   30 CONTINUE
      NLO=MS
      GO TO 15
   50 CONTINUE
      H=X(NHI)-X(NLO)
      RS1=X(NHI)-R
      RS2=R-X(NLO)
      V=(S(NLO)*RS1**3/6.+S(NHI)*RS2**3/6.+(Y(NLO)-S(NLO)*H*H/6.)
     .*RS1+(Y(NHI)-S(NHI)*H*H/6.)*RS2)/H
      RETURN
   80 CONTINUE
      IR=2
      IL=1
      IF (R.GT.X(N)) IR=N
      IF (R.GT.X(N)) IL=N-1
      SLOPE=(Y(IR)-Y(IL))/(X(IR)-X(IL))
      V=Y(IR)+SLOPE*(R-X(IR))
      RETURN
      END

      SUBROUTINE ONESPL (X,Y,S,Q,N,MXN)
c	this routine prepares some numbers used in cubic spline
c	interpolation. this routine is only called once for each
c	function X,Y. the actual interpolation is performed by
c	SUBROUTINE SPLINE, which can be called repetitively without
c	intervening calls to ONESPL.
c 
c history:
c 95/03/20 (mjw): allow one to only use part of a passed array
c                   by letting N be the number of valid indices and
c		    MXN the dimensions of the arrays
c
      implicit none
      integer N,MXN,N1,I,II
      real*8 X(MXN),Y(MXN),S(MXN),Q(MXN),H,YP1,HH,P,H1,YP

      Q(1)=0.0
      S(1)=0.0
      H=X(2)-X(1)
      YP=(Y(2)-Y(1))/H
      N1=N-1
      DO 1 I=2,N1
       H1=X(I+1)-X(I)
       YP1=(Y(I+1)-Y(I))/H1
       HH=H+H1
       P=2.0*HH+H*Q(I-1)
       Q(I)=-H1/P
       S(I)=(6.0*(YP1-YP)-H*S(I-1))/P
       H=H1
       YP=YP1
    1 CONTINUE
      S(N)=0.0
      DO 2 II=2,N1,1
       I=N+1-II
       S(I)=S(I)+Q(I)*S(I+1)
    2 CONTINUE
      RETURN
      END
