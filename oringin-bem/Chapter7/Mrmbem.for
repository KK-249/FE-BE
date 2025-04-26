C MRMBEM.FOR
C
C A
C MULTIPLE RECIPROCITY BEM PROGRAM
C FOR EIGENVALUE DETERMINATION OF THE 
C TWO DIMENSIONAL HELMHOLTZ EQUATION
C USING CONSTANT ELEMENTS
C
C THIS PROGRAM IS A FORTRAN IMPLEMENTATION OF THE EIGENVALUE
C DETERMINATION SCHEME PRESENTED IN EQUATIONS 51 THROUGH 61
C OF CHAPTER 7 OF THE BOOK ENTITLED:
C 
C BOUNDARY ELEMENT ACOUSTICS - FUNDAMENTALS AND COMPUTER CODES
C EDITED BY T.W. WU. 
C THE WIT PRESS.
C
C
C LAST REVISION: NOVEMBER 1999
C
C REMARK: NO ATTENTION HAS BEEN MADE TO OPTIMIZE THE CODE IN ANY
C         SENSE.
C
C DESCRIPTION OF THE PARAMETERS AND VARIABLES
C
C NELMAX:		MAXIMUM NUMBER OF (CONSTANT) BOUNDARY ELEMENTS	
C LSERMAX:	MAXIMUM NUMBER OF ELEMENTS IN THE MATRIX SERIES
C             SEE EQUATION (52)	
C MV:			MAXIMUM OF THE VECTOR CONTAINING THE ELEMENTS OF 
C             THE MATRIX SERIES
C NWAVEMAX:	MAXIMUM NUMBER OF FREQUENCY POINTS (WAVE NUMBERS)
C
C
	PROGRAM MRMBEM
C
      implicit double precision (a-h,o-y)
c
      parameter(nelmax=66,lsermax=80,mv=(lsermax+1)*nelmax*nelmax)
	parameter(nwavemax=4098)
c
      dimension g(nelmax,nelmax),h(nelmax,nelmax)
	dimension vg(mv),vh(mv),vfg(mv),vfh(mv)
	dimension ca(0:mv),cb(0:mv)
      dimension x(nelmax),y(nelmax),xm(nelmax),ym(nelmax)
	dimension right(nelmax)
	dimension bcond(nelmax)
	dimension vwave(nwavemax),dg(nwavemax),dgmod(nwavemax)
	dimension etadet(nwavemax)
c
      integer kodigo(nelmax)
      common /bloc1/ inp,ipr
      common /bloc2/ lser
      common /bloc3/ pi
      common /bloc4/ wave
	common /bloc5/ ncoord
	common /bloc6/ npgauss
c
	character*100 arq1,arq2
c
c names for the input and output files
c
      call sarqnat(arq1,arq2)
c
c reading input data from input file
c
      call sdatinp(nelc,nelmax,bcond,kodigo,nwave,wavefrst,wavelast,x,y)
c
c generating the mid-points of the boundary elements
c
      call smidpt(nelc,x,y,xm,ym)
c
c calculating the coefficients fo the series of fundamental solutions
c see equation (56)
c ca coefficient A(j)
c cb coefficient B(j)
c
      call scoefs(mv,ca,cb)
c
c determination of the matrices H(j) and G(j) of equation (52)
c which are stored in in the vectors vh and vg.
c
	call svetgh(nelc,mv,x,y,xm,ym,ca,cb,vh,vg)	 
c
c determination of the frequencies (wave numbers) for the loop
c
      call sfreqs(nwavemax,nwave,wavefrst,wavelast,
     *           vwave)
c
c Loop over the number of frequencies (nwave)	 
c
 	do 2001 ik=1,nwave
c
      write(*,*)"ik=",ik
c
c extracting the eigenfrequencies from n(k) 
c
      call sfilter(nelc,nelmax,mv,ik,kodigo,vg,vh,vfg,vfh,
     *            g,h,right,bcond,
     *            vwave,nwavemax,etadet,dg,dgmod) 
c
 2001 continue
c
      stop
	end
c
c**********************************************************************
      subroutine sarqnat(arq1,arq2)
c
c this subroutine asks for the name of input and output files
c
      common /bloc1/ inp,ipr
	character*100 arq1,arq2
c
c unit numbers for the input and output files
c
      inp=8
      ipr=9
c
c reading names for the input and output files
c
      write(*,40)
  40  format(t10,'write name of input file:',/)
      read(*,45)arq1
  45  format(a100)
      open(inp,file=arq1,status='old')
c
      write(*,50)
  50  format(t10,'write name of output file:',/)
      read(*,55)arq2
  55  format(a100)
      open(ipr,file=arq2,status='unknown')
c
      return
	end
c
c**********************************************************************
      subroutine sdatinp(nelc,nelmax,bcond,kodigo,
     *                   nwave,wavefrst,wavelast,x,y)
c
c this subroutine reads the input data file and writes these data to
c the output data file for checking purposes
c
      implicit double precision (a-h,o-y)
      character*100 title
c
	dimension bcond(nelmax),x(nelmax),y(nelmax)
c
      integer kodigo(nelmax)
c      
      common /bloc1/ inp,ipr
      common /bloc2/ lser
      common /bloc3/ pi
	common /bloc5/ ncoord
	common /bloc6/ npgauss
c
c  reading basic parameters
c
c nelc = number of constant boundary elements
c lser = number of elements in the mrm matrix series
c ncoord = number of coordinate used to normalize the det
c npgauss= number of gauss points for all integrations
c nwave = total of wave numbers k
c wavefrst = first wave number
c wavelast = last wave number
c
      pi=4.d0*datan(1.d0)
c         
      read (inp,'(a)') title
      write(ipr,'(a)') title
c
      read (inp,'(a)') title
      read(inp,*)      nelc,lser
      write(ipr,100)   nelc,lser
  100 format(/2x, 'number of constant boundary elements=',i4,
     */2x,'number of elements in the mrm matrix series=',i4)
c
      read (inp,'(a)') title
      read (inp,*)     ncoord,npgauss
      write(ipr,110)   ncoord,npgauss
  110 format(/ 2x, 'number of coordinate used to normalize the det=',i4,
     */2x,'number of gauss points for all integrations=',i4)
c
      read(inp,'(a)')title
      read(inp,*)nwave
      write(ipr,120)nwave
  120 format(/2x, 'total of wave numbers k=', i5)   
c
c  first and last wave numbers
c
      read(inp,'(a)') title
      read(inp,*)     wavefrst,wavelast
      write(ipr,130)  wavefrst,wavelast
  130 format(/2x, 'first wave number=', e14.7
     *       /2x, 'last wave number =', e14.7)   
c
c coordinates	of the elements x,y
c
      read(inp,'(a)')title
      write(ipr,140)
  140 format(/2x, 'boundary nodes coordinates'/2x, 'coord-number',6x,
     * 'x-coord', 19x,'y-coord')
c
	do 10 icord=1,nelc
	     read (inp,*) ic, x(icord), y(icord)
	     write(ipr,*) ic, x(icord), y(icord)
   10 continue
c
c boundary conditions
c stored in vector bcond(i).
c kodigo(i)=0 , u known
c kodigo(i)=1 , du/dn known
c
      read(inp,'(a)')title
      write(ipr,150)
  150 format(/2x, 'boundary conditions-bc'/9x, 'node',6x, 'bc-code', 7x, 
     *'bc-value')
c
      do 20 i=1,nelc
           read (inp,*)ic,kodigo(i),bcond(i)
           write(ipr,*)ic,kodigo(i),bcond(i)
   20 continue
c
      write(ipr,160)
  160 format(/2x, 'Determined Eigenfrequencies (wave number):')
c
      return
      end
c
c**********************************************************************
        subroutine smidpt(nelc,x,y,xm,ym)
c
c mid points of the boundary elements
c
        implicit double precision (a-h,o-y)
        dimension x(1),y(1),xm(1),ym(1)
c
        x(nelc+1)=x(1)
        y(nelc+1)=y(1)
        do 10 i=1,nelc
        xm(i)=(x(i)+x(i+1))/2.d0
  10    ym(i)=(y(i)+y(i+1))/2.d0
        return
        end

c
c***********************************************
	subroutine svetgh(nelc,mv,x,y,xm,ym,ca,cb,vh,vg)
c***********************************************
c
c the subroutine builds all the matrices H(j) and G(j) of the
c polynomial equation (52) for j=0,lser and stores them in the
c vectors vh and vg. They are not pre-multiplied by k**j. 
c
c the first nelc*nelc positions of the vector vg are occupied by
c the matriz g0, and so on.
c 
c non-diagonal elements are determined by subroutine sghgaus.
c diagonal elements are calculated (analytically) in sub sgijj
c
      implicit double precision (a-h,o-y)
	dimension x(1),y(1),xm(1),ym(1)
	dimension vg(mv),vh(mv),ca(0:mv),cb(0:mv)
	common /bloc2/ lser
	common /bloc1/ inp,ipr
	common /bloc3/ pi
c
      do 100 l=1,nelc*nelc*(lser+1)
	 vg(l)=0.d0
	 vh(l)=0.d0
  100 continue
c
      do 10 j=0,lser
	 do 20 i=1,nelc
	  do 30 k=1,nelc
	   ik=(k-1)*nelc+i+j*nelc*nelc
	   kk=k+1
	  if (i-k)40,45,40
   40    call sghgaus(j,mv,xm(i),ym(i),x(k),y(k),x(kk),y(kk),
     *                ca,cb,hout,gout)
         vg(ik)=gout
	   vh(ik)=hout	 
         go to 30
   45    call sgjii(j,mv,x(k),y(k),x(kk),y(kk),ca,cb,gdiag)
         vg(ik)=gdiag
	   vh(ik)=0.d0
   30   continue
   20  continue
   10 continue
c
      return
	end
c
c**************************************
      subroutine scoefs(mv,ca,cb)
c**************************************
c
c determination of the coefficients for the series of fundamental
c solutions. see equation (56). 
c a(j)=ca , b(j)=cb 
c
	implicit double precision (a-h,o-y)
      dimension ca(0:mv),cb(0:mv)
	common /bloc2/ lser 
c
	ca(0)=1.d0
	cb(0)=0.d0
      do 100 j=0,lser
	aux=ca(j)/(4*(j+1)**2)
	ca(j+1)=aux
	cb(j+1)=aux*(1.d0/(j+1)+cb(j)/ca(j))
  100 continue
c
      return
	end
c	  	   			 
c****************************************************
        subroutine sgjii(j,mv,x1,y1,x2,y2,ca,cb,gdiag)
c****************************************************
c
	implicit double precision (a-h,o-y)
      dimension ca(0:mv),cb(0:mv)        
      common /bloc3/ pi
c
c     diagonal elements of g(j)
c     see equation (58)
c 
        ax=(x2-x1)/2.d0
        ay=(y2-y1)/2.d0
        sr=dsqrt(ax**2+ay**2)
        gdiag=(sr*(sr)**(2*j))/(pi*(2*j+1))*
     *  (ca(j)*(1/(2*j+1)-dlog(sr))+cb(j))
c
        return
        end 
c
c    
c****************************************************************
        subroutine sghgaus(j,mv,xp,yp,x1,y1,x2,y2,ca,cb,hout,gout)
c****************************************************************
c
c  calculates the non-diagonal elements of all matrices H(j) and
c  G(j) of the polynomial equation (52) using the Gauss integration
c  rule. (j=0,lser)
c
	implicit double precision (a-h,o-y)
c        
      dimension xg(100),w(100)
	dimension ca(0:mv),cb(0:mv)
c
      common /bloc3/ pi
	common /bloc2/ lser
	common /bloc1/ inp,ipr
	common /bloc6/ npgauss
c
        ax=(x2-x1)/2.d0
        bx=(x2+x1)/2.d0
        ay=(y2-y1)/2.d0
        by=(y2+y1)/2.d0
c
        sr1=dsqrt(ax*ax+ay*ay)
c
c     abscissa and weights for the gauss integration
c
        call sptgaus(npgauss,xg,w)
c
c     coefficients of the h and g matrices
c     r is the normal, ra is the distance in the fund. sol.
c
        gout=0.d0
        hout=0.d0
c     
        do 15 i=1,npgauss 
        xco=bx+ax*xg(i)
        yco=by+ay*xg(i)
        r =((xco-xp)*ay-(yco-yp)*ax)/sr1
        ra=dsqrt((xco-xp)*(xco-xp)+(yco-yp)*(yco-yp))
c
        gout=gout+(-1/(2*pi)*((ra)**(2*j))*(ca(j)*
     *  dlog(ra)-cb(j)))*w(i)*sr1
c
        hout=hout+(-1/(2*pi)*((2*j*dlog(ra)+1)*ca(j)-
     *  2*j*cb(j))*(ra)**(2*j-1))*(r*w(i)*sr1)/(ra)
c 
  15    continue
c
        return
        end
c
c*****************************************
      subroutine series(nelc,mv,vg,vh,vfg,vfh)
c*****************************************
c
	implicit double precision (a-h,o-y)
      dimension vg(mv),vh(mv),vfg(mv),vfh(mv)
	common /bloc2/ lser
	common /bloc4/ wave
c
c multiplies each matrix (g0,g1,...;h0,h1,...) by k**(2j) 
c and performing the sum, storing the results in the vectors
c vfg and vfh.
c
      do 10 j=1,nelc*nelc
	vfg(j)=0.d0
 	vfh(j)=0.d0
	do 20 im=0,lser
	ck=((-1)**(im+1))*wave**(2*im)
	i=im*nelc*nelc+j
	vfg(j)=vfg(j)-ck*vg(i)
	vfh(j)=vfh(j)-ck*vh(i)
   20 continue
   10 continue
c
      return
	end
c
c*****************************************************
      subroutine smatgh(nelc,nelmax,mv,kodigo,vfg,vfh,g,h,
     *right,bcond)
c****************************************************
c
	implicit double precision (a-h,o-y)
      dimension vfg(mv),vfh(mv)
 	dimension g(nelmax,nelmax),h(nelmax,nelmax)
	dimension right(nelmax)
	dimension bcond(nelmax)
	common /bloc1/ inp,ipr
	integer kodigo(nelmax)
	common /bloc3/ pi
c
c this subroutine transforms the vectors vfg and vfh determined
c by the subroutine series, into the system matrices G and H, as
c defined in equation (53). Further it adds the free term Ci to
c the diagonal od matrix H. It also introduces the proper boundary
c conditions, forming the final system matrices A and B as defined
c in equation (54).
c
      do 1000 i=1,nelc
	do 2000 j=1,nelc
	g(i,j)=0.d0
	h(i,j)=0.d0
 2000 continue
 1000 continue
c
      do 10 i=1,nelc
 	do 20 j=1,nelc
	ij=(j-1)*nelc+i
	g(i,j)=vfg(ij)
	h(i,j)=vfh(ij)
   20 continue
   10 continue
c
c adding the free term
c
      do 30 i=1,nelc
	h(i,i)=.5d0
   30 continue
c
c introducing the boundary conditions
c
        do 45 j=1,nelc
        if (kodigo(j))45,45,40
   40   do 150 i=1,nelc
        ch=g(i,j)
        g(i,j)=-h(i,j)
        h(i,j)=-ch
  150   continue
   45   continue
c
c building the right side of the equation
c
      do 200 j=1,nelc
      right(j)=0.d0	 
    	do 250 k=1,nelc
	right(j)=right(j)+h(j,k)*bcond(k)
  250 continue
  200 continue
c
      return
	end
c	   		 
c**********************************************************************
        subroutine sptgaus(npg,xg,w)
c
	implicit double precision (a-h,o-y)
	dimension  xg(16),w(16)
c
	if(npg.eq.4) then
	xg(1) = - 0.861136311594053  
	xg(2) = - 0.339981043584856  
	xg(3) = - xg(2)
	xg(4) = - xg(1)
	w(1) = 0.347854845137454  
	w(2) = 0.652145154862546  
	w(3) = w(2)
	w(4) = w(1)
	else if(npg.eq.6) then
	xg(1) = - 0.238619186083197  
	xg(2) = - 0.661209386466256  
	xg(3) = - 0.932469514203152  
	xg(4) = - xg(1)
	xg(5) = - xg(2)
	xg(6) = - xg(3)
	w(1) = 0.467913934572691  
	w(2) = 0.360761573048139  
	w(3) = 0.171324492379170  
	w(4) = w(1)
	w(5) = w(2)
	w(6) = w(3)
	else if(npg.eq.8) then
	xg(1) = - 0.183434642495650  
	xg(2) = - 0.525532409916329  
	xg(3) = - 0.796666477413627  
	xg(4) = - 0.960289856497536  
	xg(5) = - xg(1)
	xg(6) = - xg(2)
	xg(7) = - xg(3)
	xg(8) = - xg(4)
	w(1) = 0.362683783378362  
	w(2) = 0.313706645877887  
	w(3) = 0.222381034453374  
	w(4) = 0.101228536290376  
	w(5) = w(1)
	w(6) = w(2)
	w(7) = w(3)
	w(8) = w(4)
	else if(npg.eq.16) then 
	xg(1) = - 0.0950125098376374  
	xg(2) = - 0.2816035507792589  
	xg(3) = - 0.4580167776572273  
	xg(4) = - 0.6178762444026437  
	xg(5) = - 0.7554044083550030  
	xg(6) = - 0.8656312023878317  
	xg(7) = - 0.9445750230732325  
	xg(8) = - 0.9894009349916499  
	xg(9) = - xg(1)
	xg(10) = - xg(2)
	xg(11) = - xg(3)
	xg(12) = - xg(4)
	xg(13) = - xg(5)
	xg(14) = - xg(6)
	xg(15) = - xg(7)
	xg(16) = - xg(8)
	w(1) = 0.18945061045506  
	w(2) = 0.18260341504492  
	w(3) = 0.16915651939500  
	w(4) = 0.14959598881657  
	w(5) = 0.12462897125553  
	w(6) = 0.09515851168249  
	w(7) = 0.06225352393864  
	w(8) = 0.02715245941175  
	w(9) = w(1)
	w(10) = w(2)
	w(11) = w(3)
	w(12) = w(4)
	w(13) = w(5)
	w(14) = w(6)
	w(15) = w(7)
	w(16) = w(8)
c
	end if
c
	return
	end
c
c
c************************************
      SUBROUTINE SOLGAUS(A,B,DET,N,NX)
c************************************	
C
C SOLUTION OF LINEAR SYSTEMS BY GAUSS ELIMINATION
C
C A = SYSTEM MATRIX
C B = IN THE INPUT PHASE CONTAINS THE RIGHT SIDE
C     IN THE OUTPUT PHASE CONTAINS THE SOLUTION
C DET= DETERMINANT OF THE SYSTEM, OUTPUT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
	DIMENSION A(NX,NX),B(NX)
	common /bloc1/ inp,ipr
C
      TOL=1.D-12
	N1=N-1
	DO 100 K=1,N1
	K1=K+1
	C=A(K,K)
	IF (DABS(C)-TOL)1,1,3
   1  DO 7 J=K1,N
C
      IF (DABS((A(J,K)))-TOL)7,7,5
   5  DO 6 L=K,N
      C=A(K,L)
	A(K,L)=A(J,L)
   6  A(J,L)=C
      C=B(K)
	B(K)=B(J)
	B(J)=C
	C=A(K,K)
c
	GO TO 3
   7  CONTINUE
      GO TO 8
C
   3  DO 4 J=K1,N
      A(K,J)=A(K,J)/C
   4	continue
C
      DO 10 I=K1,N
	C=A(I,K)
	DO 9 J=K1,N
   9  A(I,J)=A(I,J)-C*A(K,J)
  10  B(I)=B(I)-C*B(K)
  100 CONTINUE
C
      IF (DABS((A(N,N)))-TOL)8,8,101
C  101 B(N)=B(N)/A(N,N)
  101 B(N)=1.D0
C
      DO 200 L=1,N1
	K=N-L
	K1=K+1
	DO 200 J=K1,N
  200 B(K)=B(K)-A(K,J)*B(J)
C
      DET=1.D0
	DO 250 I=1,N
  250 DET=DET*A(I,I)
      GO TO 300
   8  WRITE(IPR,2)K
   2  FORMAT('**** SINGULARITY IN LINE', I5)
      DET=0.D0
  300 RETURN
      END
c
c**********************************************************************
c
	subroutine subdet(nelc,nelmax,mv,ik,kodigo,vg,vh,vfg,vfh,
     *                  g,h,right,bcond,gmod,nwavemax,
     *                  etadet,dg,dgmod) 
c
c this subroutine takes matrix [A] from the subroutine smatgh, builds
c the matrix [R] using equation (61) and the normalizing coordinate
c ncoord (defined in the input file), calculates det[A] and det[R]
c using subroutine solgaus. In the sequence it builds the function
c n(k)=(etadet(k)) as described in equation (60)
c
      implicit  double precision (a-h,o-y)
	integer   kodigo(nelmax)
	dimension vg(mv),vh(mv),vfg(mv),vfh(mv)
	dimension bcond(nelmax),right(nelmax)
     	dimension g(nelmax,nelmax),h(nelmax,nelmax),gmod(nelmax,nelmax)
	dimension dg(nwavemax),dgmod(nwavemax),etadet(nwavemax)
	common    /bloc4/ wave
	common    /bloc5/ ncoord
c 	
      call series(nelc,mv,vg,vh,vfg,vfh)
c
      call smatgh(nelc,nelmax,mv,kodigo,vfg,vfh,g,h,right,bcond)
c 
      do 1020 igmod=1,nelc
	do 1020 jgmod=1,nelc
	gmod(igmod,jgmod)=0.0d0
	gmod(igmod,jgmod)=g(igmod,jgmod)
 1020 continue
c
c normalizing coordinate, equation (59)
c
	do 1030 i1mod=1,nelc
	gmod(i1mod,ncoord)=-h(i1mod,ncoord)
 1030 continue
c
      call solgaus(g,right,detg,nelc,nelmax)
	call solgaus(gmod,right,detgmod,nelc,nelmax)
c
      dg(ik)=detg
	dgmod(ik)=detgmod
	etadet(ik)=dg(ik)/dgmod(ik)
c
c
	return
	end
c**********************************************************************
      subroutine soutfnt(raiz)
c
c this subroutine prints to the output file  the eigenfrequencies
c that have been found in the specified frequency range 
c
	implicit double precision (a-h,o-y)
	common /bloc1/ inp,ipr
c
      write(ipr,335)raiz
  335 format(e14.7)
c
	return
	end
c**********************************************************************
      subroutine sfilter(nelc,nelmax,mv,ik,kodigo,vg,vh,vfg,vfh,
     *            g,h,right,bcond,
     *            vwave,nwavemax,etadet,dg,dgmod)  
c
c this subroutine takes the value of n(k)=etadet(k) from equation (60)
c calculated in subroutine subdet and extracts the eigenfrequencies.
c 
c the strategy is based on the change of sign of n(k)=etadet(k).
c once there is a change of sign it must be checked if it comes from
c a sign change in det[A] or det[R]. Only zero changes in det[A] are
c related to eigenfrequencies 
c
	implicit double precision (a-h,o-y)
c
	integer kodigo(nelmax)
	dimension vg(mv),vh(mv),vfg(mv),vfh(mv)
	dimension bcond(nelmax)
	dimension right(nelmax)
     	dimension g(nelmax,nelmax),h(nelmax,nelmax),gmod(nelmax,nelmax)
	dimension vwave(nwavemax)
	dimension dg(nwavemax),dgmod(nwavemax),etadet(nwavemax)
c
	common /bloc4/ wave
	common /bloc1/ inp,ipr
c
	etadet(ik)=0.d0
	dg(ik)=0.d0
	dgmod(ik)=0.d0
c 	
      wave1=vwave(ik)
	wave=wave1
c
	call subdet(nelc,nelmax,mv,ik,kodigo,vg,vh,vfg,vfh,
     *            g,h,right,bcond,gmod,nwavemax,
     *            etadet,dg,dgmod)
c
c
      if (ik .eq. 1) then
	return
	else
	prod=etadet(ik-1)*etadet(ik)
	endif
c
        if (prod .gt. 0.d0) then
           return
	  else
	  detg1=dabs(dg(ik-1))
	  detgmod1=dabs(dgmod(ik-1))
	  detg2=dabs(dg(ik))
	  detgmod2=dabs(dgmod(ik))
	     if (detg1 .gt. detgmod1 .and. detg2 .gt. detgmod2) then
		    return
		 else if (detg1 .lt. detgmod1 .and. detg2 .lt. detgmod2) then 
		 etadet1=dabs(etadet(ik-1))
		 etadet2=dabs(etadet(ik))
		    if (etadet1 .lt. etadet2) then
			 raiz=vwave(ik-1)
			else
			 raiz=vwave(ik)
			endif
	        call soutfnt(raiz)
		 endif
        endif
c	  		 	  		 	 
	return
	end
c**********************************************************************
	subroutine sfreqs(nwavemax,nwave,wavefrst,wavelast,
     *           vwave) 
c
c this subroutine determines the wave numbers for which the 
c determinants will be calculated
c
      implicit double precision (a-h,o-y)
	dimension vwave(nwavemax)
c
      if(nwave.eq.1)then
	wstep=0
	else   
      wstep=(wavelast-wavefrst)/(nwave-1)
	endif
      do 20 ifreq=1,nwave
 20   vwave(ifreq)=wavefrst+wstep*(ifreq-1)
c
	return
	end
c**********************************************************************
