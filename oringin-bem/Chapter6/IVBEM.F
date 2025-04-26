C234567890
      implicit double precision (a-h,o-z)
      complex*16 vn
c nodal coordinates
      dimension xnd(500), ynd(500), znd(500)
c element connectivity
      dimension nel1(500), nel2(500), nel3(500)
c coordinates for field points
      dimension xf(100), yf(100), zf(100)
c complex element velocity
      dimension vn(500)
c nnodes = number of nodes in bem model
c nelem = number of elements in bem model
c nfp = number of field points
c
c bem.model contains the bem model in NASTRAN format
c bem.vel.bc contains the velocity boundary conditions assigned
c            on each element
c field.model contains the coordinates of the field points in NASTRAN format
      open (unit=15, file='bem.model', status='UNKNOWN')
      open (unit=17, file='bem.vel.bc', status='UNKNOWN')
      open (unit=18, file='field.model', status='UNKNOWN')
      open (unit=20, file='output', status='UNKNOWN')
c intitialize arrays
      n1=500
      n2=100
      call intiz(n1,nel1)
      call intiz(n1,nel2)
      call intiz(n1,nel3)
      call realiz(n1,xnd)
      call realiz(n1,ynd)
      call realiz(n1,znd)
      call realiz(n2,xf)
      call realiz(n2,yf)
      call realiz(n2,zf)
      call compliz(n1,vn)
c
      call rdbem(nnodes,nelem,xnd,ynd,znd,nel1,nel2,nel3)
      call rdfield(nfield,xf,yf,zf)
      call rdvel(nelem,vn)
      call indbem(nnodes,nelem,nfield,xnd,ynd,znd,nel1,nel2,nel3
     *,xf,yf,zf,vn)
      stop 
      end
c
c
c
      subroutine indbem(nnodes,nelem,nfield,xnd,ynd,znd
     *,nel1,nel2,nel3,xf,yf,zf,vn)
      implicit double precision (a-h,o-z)
      complex*16 vn, fv, gr,  cro, prod1, prod2, prod3, amat
      complex*16 prod4
      complex*16 sc1, sc2, p, press, prs, tmat
      dimension fv(500), amat(500,500), tmat(500,500)
      dimension xnd(500), ynd(500), znd(500)
      dimension nel1(500), nel2(500), nel3(500), iv(500)
      dimension xf(100), yf(100), zf(100)
      dimension vn(500)
c
      dimension iel1(3),iel2(3)
      dimension xx1(3), yy1(3), zz1(3), rr1(3,3), fn1(3), g1(3)
      dimension xx2(3), yy2(3), zz2(3), rr2(3,3), fn2(3), g2(3)
      dimension xd1(3), yd1(3), zd1(3), iint1(3)
      dimension p(3)
c
      write(*,*)'insert f'
      read(*,*)f
      write(*,*)'insert rho'
      read(*,*)rho
      write(*,*)'insert csound'
      read(*,*)csound
      pi=3.141d0
      fk=2.0d0*pi*f/csound
      cro=dcmplx(0.0d0,rho*2.0d0*pi*f)
c
      do 10 i=1,nelem
       iel1(1)=nel1(i)
       iel1(2)=nel2(i)
       iel1(3)=nel3(i)
       do 20 ii=1,3
        xx1(ii)=xnd(iel1(ii))
        yy1(ii)=ynd(iel1(ii))
        zz1(ii)=znd(iel1(ii))
 20    continue
       call elem(rr1,xx1,yy1,zz1,fn1,ar1)
       do 30 int=1,7
        call integ(g1,w1,int)
        do 40 j=1,3
         fv(iel1(j))=-vn(i)*ar1*w1*g1(j)*cro+fv(iel1(j))
 40     continue
 30    continue
       do 100 j=1,nelem
        if(i.ne.j)then
         iel2(1)=nel1(j)
         iel2(2)=nel2(j)
         iel2(3)=nel3(j)
         do 110 ii=1,3
          xx2(ii)=xnd(iel2(ii))
          yy2(ii)=ynd(iel2(ii))
          zz2(ii)=znd(iel2(ii))
 110     continue
         call elem(rr2,xx2,yy2,zz2,fn2,ar2)
         prd=fn1(1)*fn2(1)+fn1(2)*fn2(2)+fn1(3)*fn2(3)
         do 120 int1=1,7
          call integ(g1,w1,int1)
          x1=0.0d0
          y1=0.0d0
          z1=0.0d0
          do 130 ii=1,3
           x1=x1+xx1(ii)*g1(ii)
           y1=y1+yy1(ii)*g1(ii)
           z1=z1+zz1(ii)*g1(ii)
 130      continue
          s1=w1*ar1
          do 140 int2=1,7
           call integ(g2,w2,int2)
           x2=0.0d0
           y2=0.0d0
           z2=0.0d0
           do 150 ii=1,3
            x2=x2+xx2(ii)*g2(ii)
            y2=y2+yy2(ii)*g2(ii)
            z2=z2+zz2(ii)*g2(ii)
 150       continue
           r=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
           gr=dcmplx(dcos(r*fk),-dsin(r*fk))/(4.0d0*pi*r)
           prod1=s1*w2*ar2*gr
           prod4=prod1*prd
           do 160 ip=1,3
            do 170 jp=1,3
             apr=0.0d0
             do 180 ijp=1,3
              apr=apr+rr1(ijp,ip)*rr2(ijp,jp)
 180         continue
             prod2=prod4*g1(ip)*g2(jp)*fk*fk
             prod3=(prod1*apr)
c             if(i.eq.1)then
c             endif
             amat(iel1(ip),iel2(jp))=amat(iel1(ip),iel2(jp))
     *+prod2-prod3
 170        continue
 160       continue
 140      continue
 120     continue
        else
         call sgint1(xx1,yy1,zz1,iel1)
         call sgint2(xx1,yy1,zz1,xxx,yyy,zzz)
         call elem(rr1,xx1,yy1,zz1,fn1,ar1)
         do 200 m1=1,2
          if(m1.eq.1)then
           xd1(1)=xx1(1)
           yd1(1)=yy1(1)
           zd1(1)=zz1(1)
           xd1(2)=xxx
           yd1(2)=yyy
           zd1(2)=zzz
           xd1(3)=xx1(3)
           yd1(3)=yy1(3)
           zd1(3)=zz1(3)
          endif
          if(m1.eq.2)then
           xd1(1)=xxx
           yd1(1)=yyy
           zd1(1)=zzz
           xd1(2)=xx1(2)
           yd1(2)=yy1(2)
           zd1(2)=zz1(2)
           xd1(3)=xx1(3)
           yd1(3)=yy1(3)
           zd1(3)=zz1(3)
          endif
          call elem(rr2,xd1,yd1,zd1,fn2,ar2)
          do 210 int1=1,7
           call sgint3(g1,w1,int1)
           pr=w1*ar1
           x1=0.0d0
           y1=0.0d0
           z1=0.0d0
           do 220 ind=1,3
            x1=x1+xx1(ind)*g1(ind)
            y1=y1+yy1(ind)*g1(ind)
            z1=z1+zz1(ind)*g1(ind)
 220       continue
           do 230 int2=1,7
            call integ(g2,w2,int2)
            x2=0.0d0
            y2=0.0d0
            z2=0.0d0
            do 240 ind=1,3
             x2=x2+xd1(ind)*g2(ind)
             y2=y2+yd1(ind)*g2(ind)
             z2=z2+zd1(ind)*g2(ind)
 240        continue
            r=dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            gr=dcmplx(dcos(fk*r),-dsin(fk*r))/(4.0d0*pi*r)
            prod1=pr*w2*ar2*gr
            do 250 ip=1,3
             do 260 jp=1,3
              apr=0.0d0
              do 270 ijp=1,3
               apr=apr+rr1(ijp,ip)*rr2(ijp,jp)
 270          continue
              prod2=prod1*g1(ip)*g2(jp)*fk*fk
              prod3=prod1*apr
              if(m1.eq.1)then
               if(jp.eq.2)then
                amat(iel1(ip),iel1(1))=amat(iel1(ip),iel1(1))+
     *prod2*0.50d0-prod3*0.50d0
                amat(iel1(ip),iel1(2))=amat(iel1(ip),iel1(2))+
     *prod2*0.50d0-prod3*0.50d0
               else
                amat(iel1(ip),iel1(jp))=amat(iel1(ip),iel1(jp))+
     *prod2-prod3
               endif
              else
               if(jp.eq.1)then
                amat(iel1(ip),iel1(1))=amat(iel1(ip),iel1(1))+
     *prod2*0.50d0-prod3*0.50d0
                amat(iel1(ip),iel1(2))=amat(iel1(ip),iel1(2))+
     *prod2*0.50d0-prod3*0.50d0
               else
                amat(iel1(ip),iel1(jp))=amat(iel1(ip),iel1(jp))+
     *prod2-prod3
               endif
              endif
 260         continue
 250        continue
 230       continue
 210      continue
 200     continue
        endif
 100   continue
 10   continue
c
      ndm=500
c
c
      call cdlud(nnodes,ndm,amat,ndm,tmat,iv)
      call cdbs(nnodes,ndm,tmat,iv,fv)
c
      do 300 i=1,nfield
       press=dcmplx(0.0d0,0.0d0)
       fdx=xf(i)
       fdy=yf(i)
       fdz=zf(i)
       do 310 j=1,nelem
        iel1(1)=nel1(j)
        iel1(2)=nel2(j)
        iel1(3)=nel3(j)
        do 320 ii=1,3
         xx1(ii)=xnd(iel1(ii))
         yy1(ii)=ynd(iel1(ii))
         zz1(ii)=znd(iel1(ii))
         p(ii)=fv(iel1(ii))
 320    continue
        call elem(rr1,xx1,yy1,zz1,fn1,ar1)
        do 330 int1=1,7
         call integ(g1,w1,int1)
         x1=0.0d0
         y1=0.0d0
         z1=0.0d0
          do 340 ii=1,3
           x1=x1+xx1(ii)*g1(ii)
           y1=y1+yy1(ii)*g1(ii)
           z1=z1+zz1(ii)*g1(ii)
 340      continue
         dx=x1-fdx
         dy=y1-fdy
         dz=z1-fdz 
         r=dsqrt(dx**2+dy**2+dz**2)
         gr=dcmplx(dcos(r*fk),-dsin(r*fk))/(4.0d0*pi*r)
         fnor=(fn1(1)*dx+fn1(2)*dy+fn1(3)*dz)/r
         sc1=dcmplx(-1.0d0,-fk*r)
         sc2=dcmplx(dcos(fk*r),-dsin(fk*r))
         dd=pi*r*r*4.0d0
         do 350 ii=1,3
          prs=(sc1*sc2*fnor*ar1*p(ii)*w1*g1(ii))/dd
          press=press+prs
 350     continue
 330    continue
 310   continue
       pmag=cdabs(press)
       db=20.0d0*dlog10(pmag/(dsqrt(2.0d0)*2.0d0*1.0e-5))
       write(20,1000)i,db,pmag,press
 1000  format(i5,4f10.4)
 300  continue
      return
      end
c
      subroutine elem(rr,x,y,z,f,ar)
      implicit double precision (a-h,o-z)
      dimension rr(3,3), x(3), y(3), z(3), f(3)
      f(1)=((y(1)-y(3))*(z(2)-z(3)))-((y(2)-y(3))*(z(1)-z(3)))
      f(2)=((z(1)-z(3))*(x(2)-x(3)))-((z(2)-z(3))*(x(1)-x(3)))
      f(3)=((x(1)-x(3))*(y(2)-y(3)))-((x(2)-x(3))*(y(1)-y(3)))
      ar=0.50d0*dsqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
      ar2=2.0d0*ar
      f(1)=f(1)/ar2
      f(2)=f(2)/ar2
      f(3)=f(3)/ar2
      rr(1,1)=(x(2)-x(3))/ar2
      rr(1,2)=(x(3)-x(1))/ar2
      rr(1,3)=(x(1)-x(2))/ar2
      rr(2,1)=(y(2)-y(3))/ar2
      rr(2,2)=(y(3)-y(1))/ar2
      rr(2,3)=(y(1)-y(2))/ar2
      rr(3,1)=(z(2)-z(3))/ar2
      rr(3,2)=(z(3)-z(1))/ar2
      rr(3,3)=(z(1)-z(2))/ar2
      return
      end
c
c
c
      subroutine integ(s,w,n)
      implicit double precision (a-h,o-z)
      dimension s(3)
      go to (1,2,3,4,5,6,7)n
 1    s(1)=1.0d0/3.0d0
      s(2)=1.0d0/3.0d0
      w=0.2250d0
      go to 8
 2    s(1)=0.79742699d0
      s(2)=0.10128651d0
      w=0.12593918d0
      go to 8
 3    s(1)=0.10128651d0
      s(2)=0.79742699d0
      w=0.12593918d0
      go to 8
 4    s(1)=0.10128651d0
      s(2)=0.10128651d0
      w=0.12593918d0
      go to 8
 5    s(1)=0.05971587d0
      s(2)=0.47014206d0
      w=0.13239415d0
      go to 8
 6    s(1)=0.47014206d0
      s(2)=0.05971587d0
      w=0.13239415d0
      go to 8
 7    s(1)=0.47014206d0
      s(2)=0.47014206d0
      w=0.13239415d0
 8    s(3)=1.0d0-s(1)-s(2)
      return
      end
c
c
c
      subroutine sgint1(x,y,z,iel)
      implicit double precision (a-h,o-z)
      dimension x(3), y(3), z(3), iel(3)
      dimension x2(3), y2(3), z2(3), iel2(3)
      s1=dsqrt((x(2)-x(3))**2+(y(2)-y(3))**2+(z(2)-z(3))**2)
      s2=dsqrt((x(1)-x(3))**2+(y(1)-y(3))**2+(z(1)-z(3))**2)
      s3=dsqrt((x(1)-x(2))**2+(y(1)-y(2))**2+(z(1)-z(2))**2)
      smax=max(s1,s2,s3)
      if(dabs(s1-smax).lt.0.00001*s1)then
       x2(1)=x(2)
       y2(1)=y(2)
       z2(1)=z(2)
       iel2(1)=iel(2)
       x2(2)=x(3)
       y2(2)=y(3)
       z2(2)=z(3)
       iel2(2)=iel(3)
       x2(3)=x(1)
       y2(3)=y(1)
       z2(3)=z(1)
       iel2(3)=iel(1)
       go to 10
      endif
      if(dabs(s2-smax).lt.0.00001*s2)then
       x2(1)=x(3)
       y2(1)=y(3)
       z2(1)=z(3)
       iel2(1)=iel(3)
       x2(2)=x(1)
       y2(2)=y(1)
       z2(2)=z(1)
       iel2(2)=iel(1)
       x2(3)=x(2)
       y2(3)=y(2)
       z2(3)=z(2)
       iel2(3)=iel(2)
       go to 10
      endif
      return
 10   continue
      do 20 i=1,3
      x(i)=x2(i)
      y(i)=y2(i)
      z(i)=z2(i)
      iel(i)=iel2(i)
 20   continue
      return
      end
c
c
c
      subroutine sgint2(x,y,z,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      dimension x(3), y(3), z(3)
      xxx=x(1)+(x(2)-x(1))/2.0d0
      yyy=y(1)+(y(2)-y(1))/2.0d0
      zzz=z(1)+(z(2)-z(1))/2.0d0
      return
      end
c
c
c
      subroutine sgint3(s,w,n)
      implicit double precision (a-h,o-z)
      dimension s(3)
      go to (1,2,3,4,5,6,7)n
 1    s(1)=1.0d0/3.0d0
      s(2)=1.0d0/3.0d0
      w=27.0d0/60.0d0
      go to 8
 2    s(1)=0.5d0
      s(2)=0.5d0
      w=8.0d0/60.0d0
      go to 8
 3    s(1)=0.0d0
      s(2)=0.50d0
      w=8.0d0/60.0d0
      go to 8
 4    s(1)=0.5d0
      s(2)=0.0d0
      w=8.0d0/60.0d0
      go to 8
 5    s(1)=1.0d0
      s(2)=0.0d0
      w=3.0d0/60.0d0
      go to 8
 6    s(1)=0.0d0
      s(2)=1.0d0
      w=3.0d0/60.0d0
      go to 8
 7    s(1)=0.0d0
      s(2)=0.0d0
      w=3.0d0/60.0d0
 8    s(3)=1.0d0-s(1)-s(2)
      return
      end
c
c
c
      subroutine rdvel(nelem,vn)
      implicit double precision (a-h,o-z)
      complex*16 vn
      dimension vn(500)
      in=0
      do 10 i=1,nelem
      read(17,'(6x,2f10.4)')vn(i)
 10   continue
      return
      end
c
c
c
      subroutine rdfield(nfield,xf,yf,zf)
      implicit double precision (a-h,o-z)
      dimension xf(100), yf(100), zf(100)
      character*80 line
      character*4 grid
      nfield=0
 10   continue
      read(18,'(a80)',end=100)line
      read(line,'(a4)')grid
      if(grid.eq.'GRID')then
         nfield=nfield+1
         read(line,'(16x,8x,3f8.3)')xf(nfield),yf(nfield),
     *zf(nfield)
      endif
      go to 10
 100  continue
      return
      end
c
c
c
      subroutine rdbem(nnodes,nelem,xnd,ynd,znd,nel1,nel2,nel3)
      implicit double precision (a-h,o-z)
      dimension xnd(500), ynd(500), znd(500)
      dimension nel1(500), nel2(500), nel3(500)
      character*80 line
      character*4 grid
      character*6 tria
c
      nnodes=0
      nelem=0
 10   continue
      read(15,'(a80)',end=100)line
      read(line,'(a4)')grid
      read(line,'(a6)')tria
      if(grid.eq.'GRID')then
        nnodes=nnodes+1
        read(line,'(16x,8x,3f8.3)')xnd(nnodes),ynd(nnodes),
     *znd(nnodes)
      endif
      if(tria.eq.'CTRIA3')then
        nelem=nelem+1
        read(line,'(16x,8x,3i8)')nel1(nelem),nel2(nelem),nel3(nelem)
      endif
      go to 10
 100  continue
      return
      end
c
c
c
      subroutine intiz(n,iv)
      dimension iv(n)
      do 10 i=1,n
      iv(i)=0
 10   continue
      return
      end
c
c
c
      subroutine realiz(n,v)
      implicit double precision (a-h,o-z)
      dimension v(n)
      do 10 i=1,n
      v(i)=0.0d0
 10   continue
      return
      end
c
c
c
      subroutine compliz(n,v)
      implicit double precision (a-h,o-z)
      complex*16 v
      dimension v(n)
      do 10 i=1,n
      v(i)=dcmplx(0.0d0,0.0d0)
 10   continue
      return
      end
c
c
c
      subroutine cdlud(n,adim,a,tdim,t,iv)
      integer i,j,k,kp1,l
      integer n,adim,tdim,iv(1)
      double precision piv
      complex*16 tmp
      complex*16 a(adim,n), t(tdim,n)
c
      if(adim.lt.n)go to 8110
       do 8100 j=1,n
       do 8100 i=1,n
 8100  t(i,j)=a(i,j)
 8110 continue
      iv(n)=1
      do 8260 k=1,n
       piv=cdabs(t(k,k))
       if(k.ge.n)go to 8260
        l=k
        kp1=k+1
        do 8200 i=kp1,n
         if(piv.ge.cdabs(t(i,k)))go to 8200
          piv=cdabs(t(i,k))
          l=i
 8200    continue
         iv(k)=l
         tmp=t(l,k)
         if(k.ne.l)go to 8210
          if(piv.gt.0.0d0)go to 8220
          iv(n)=0
          go to 8260
 8210    continue
         iv(n)=-iv(n)
         t(l,k)=t(k,k)
         t(k,k)=tmp
 8220    continue
         tmp=-tmp
         do 8230 i=kp1,n
 8230    t(i,k)=t(i,k)/tmp
         do 8250 j=kp1,n
          tmp=t(l,j)
          t(l,j)=t(k,j)
          t(k,j)=tmp
          do 8240 i=kp1,n
 8240     t(i,j)=t(i,j)+t(i,k)*tmp
 8250     continue
 8260     continue
          if(piv.eq.0.0d0)iv(n)=0
          return
          end
c
c
c
          subroutine cdbs(n,tdim,t,iv,b)
          integer n, tdim, iv(1)
          complex*16 t(tdim,n), b(1)
          integer i,k,kp1,l
          complex*16 tmp
          do 8110 k=1,n
           if(k.ge.n)go to 8110
           l=iv(k)
           tmp=b(l)
           b(l)=b(k)
           b(k)=tmp
           kp1=k+1
           do 8100 i=kp1,n
 8100      b(i)=b(i)+t(i,k)*tmp
 8110      continue
           k=n
 8200      b(k)=b(k)/t(k,k)
           if(k.le.1)return
           tmp=-b(k)
           kp1=k
           k=k-1
           do 8210 i=1,k
 8210      b(i)=b(i)+t(i,kp1)*tmp
           go to 8200
           end
