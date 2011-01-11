ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ... SPLINE INTERPOLATION
c ... THESE INTERPOLATES BY FINDING CUBIC POLYNOMIALS THAT 
c ... REPRODUCE THE VALUES OF Y AT EACH X VALUE
c ... AND HAVE CONTINUOUS DERIVATIVES AT THESE POINTS 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ... THIS ROUTINE PREPARES THE INTERPOLATION VECTORS
c ... x, y, n are the n initial values
c ... yp1 and ypn choose the values of the endpoints
c ... y2 stores the vector of derivatives
      SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)

      PARAMETER (NMAX=5000)

      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)

      if(yp1.gt..99e30) then
         y2(1)=0.d0
         u(1)=0.d0
      else
         y2(1)=-0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do 11 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
 11   enddo

      if(ypn.gt..99e30) then
         qn=0.d0
         un=0.d0
      else
         qn=0.5d0
         un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

      do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 12   enddo

      RETURN
      END

c ... THIS ROUTINE USES THE PREVIOUSLE PREPARED VECTORS
c ... TO ACTUALLY CARRY OUT THE INTERPOLATION
c ... xa, ya, y2a are the n initial values for x,y and derivatives
c ... x is the coordinate where the interpolation is done
c ... y is the interpolated value
      SUBROUTINE SPLINT(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x) then
            khi=k
         else
            klo=k
         endif
       go to 1
       endif

       h=xa(khi)-xa(klo)
       if(h.eq.0.d0) then
        write(*,*) 'bad xa input'
        read(*,*)
       endif
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+
     &  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
       RETURN
       END

c ... THIS ROUTINE USES THE PREVIOUSLE PREPARED VECTORS
c ... TO ACTUALLY CARRY OUT THE INTERPOLATION
c ... xa, ya, y2a are the n initial values for x,y and derivatives
c ... x is the coordinate where the interpolation is done
c ... y is the interpolated value
c ... The extra variable ki is used because sometimes the x values
c ... are consecutive and one does not need to explore the whole initial
c ... x values to find the closest one. ki should be 1 at the first 
c ... iteration. Then it is stored in the program and used to see if 
c ... the x is close to the previously used values. 
       SUBROUTINE SPLIN2(xa,ya,y2a,n,x,y,ki)
       INTEGER n,ki
       DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
       INTEGER k,khi,klo
       DOUBLE PRECISION a,b,h
       klo=1
       khi=n

       if(ki.eq.1) go to 1

       if(xa(ki).lt.x.and.xa(ki+1).gt.x) then
          klo=ki
          khi=klo+1
          go to 2
       else if(xa(ki+1).lt.x .and. xa(ki+2).gt.x) then
          klo=ki+1
          khi=klo+1
          go to 2
       endif

  1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x) then
             khi=k
          else
             klo=k
          endif
       go to 1
       endif

  2    continue
       ki=klo
       h=xa(khi)-xa(klo)
       if(h.eq.0.d0) then
        write(*,*) 'bad xa input'
        read(*,*)
       endif
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+
     &  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
       RETURN
       END

c ... USER FRIENDLY SPLINE
c ... One only enters the initial and final vectors here
      SUBROUTINE SPINT(xa,ya,na,xb,yb,nb)
       INTEGER na,nb,ib,iw
       DOUBLE PRECISION xa(na),ya(na),y2a(na),xb(nb),yb(nb),yn
       yn=5e30
       call spline(xa,ya,na,yn,yn,y2a)

       iw=1
       do ib=1,nb
         call splin2(xa,ya,y2a,na,xb(ib),yb(ib),iw)
       enddo

      RETURN
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ... LINEAR INTERPOLATION

      SUBROUTINE LIN_INT(xa,ya,n,x,y,ki)
      INTEGER n,ki
      DOUBLE PRECISION x,y,xa(n),ya(n)

      INTEGER k,khi,klo
      DOUBLE PRECISION h
      klo=1
      khi=n

      if(ki.eq.1) go to 1

      if(xa(ki).lt.x.and.xa(ki+1).gt.x) then
         klo=ki
         khi=klo+1
         go to 2
      else if(xa(ki+1).lt.x .and. xa(ki+2).gt.x) then
         klo=ki+1
         khi=klo+1
         go to 2
      endif

 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
          if(xa(k).gt.x) then
             khi=k
          else
             klo=k
          endif
       go to 1
       endif

 2     continue
       ki=klo

       h=xa(khi)-xa(klo)
       if(h.eq.0.d0) then
        write(*,*) 'bad xa input'
        read(*,*)
       endif

       y=ya(klo)+(ya(khi)-ya(klo))/(xa(khi)-xa(klo))*(x-xa(klo))

       RETURN
       END

c ... USER FRIENDLY LINEAR INTERPOLATION
c ... One only enters the initial and final vectors here
      SUBROUTINE LININT(xa,ya,na,xb,yb,nb)
       INTEGER na,nb,ib,iw
       DOUBLE PRECISION xa(na),ya(na),xb(nb),yb(nb)

       iw=1
       do ib=1,nb
         call lin_int(xa,ya,na,xb(ib),yb(ib),iw)
       enddo

      RETURN
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ... 2D INTERPOLATIONS

c ... LINEAR
      SUBROUTINE LIN_INT2D(xa,nx,ya,ny,za,x,y,z,kx,ky)
      INTEGER nx,ny,kx,ky
      DOUBLE PRECISION xa(nx),ya(ny),za(nx,ny),x,y,z
     
      INTEGER kxlo,kxhi,kylo,kyhi,ix,iy,k
      DOUBLE PRECISION z1,z2,z3,z4,t,u

      kxlo=1
      kxhi=nx

      kylo=1
      kyhi=ny

      ix=0
      iy=0

      if(kx.eq.1 .and. ky.eq.1) go to 100

c ... Try to find the index which bounds the x and y variable from the previous interpolations
      if(xa(kx).lt.x .and. xa(kx+1).gt.x) then
         kxlo=kx
         kxhi=kxlo+1
         ix=1
      else if(xa(kx+1).lt.x .and. xa(kx+2).gt.x) then
         kxlo=kx+1
         kxhi=kxlo+1
         ix=1
      endif

      if(ya(ky).lt.y .and. ya(ky+1).gt.y) then
         kylo=ky
         kyhi=kylo+1
         iy=1
      else if(ya(ky+1).lt.y .and. ya(ky+2).gt.y) then
         kylo=ky+1
         kyhi=kylo+1
         iy=1
      endif

 100  continue

c ... If x and y are not close to the previous indices, find new indices
      if(ix.eq.0) then
 1       if (kxhi-kxlo.gt.1) then
            k=(kxhi+kxlo)/2
            if(xa(k).gt.x) then
             kxhi=k
            else
             kxlo=k
            endif
            go to 1
         endif
      endif

      if(iy.eq.0) then
 2       if (kyhi-kylo.gt.1) then
            k=(kyhi+kylo)/2
            if(ya(k).gt.y) then
             kyhi=k
            else
             kylo=k
            endif
            go to 2
         endif
      endif

c ... The indices are known. Prepare to Interpolate
      z1=za(kxlo,kylo)
      z2=za(kxhi,kylo)
      z3=za(kxhi,kyhi)
      z4=za(kxlo,kyhi)

      t=(x-xa(kxlo))/(xa(kxhi)-xa(kxlo))
      u=(y-ya(kylo))/(ya(kyhi)-ya(kylo))

      z= (1.d0-t)*(1.d0-u)*z1 
     & + t*(1.d0-u)*z2
     & + t*u*z3
     & + (1.d0-t)*u*z4


      RETURN
      END

c ... SPLINE PREPARATION IN 2D
c ... PREPARES SPLINES ON THE ROWS OF THE za VARIABLE
      SUBROUTINE SPLINE2D(nx,ya,ny,za,z2a)
      INTEGER nx,ny !,NN
      DOUBLE PRECISION ya(ny),za(nx,ny),z2a(nx,ny)
      
      INTEGER j,k
      DOUBLE PRECISION yspl,y2tmp(ny),ytmp(ny)

      yspl=5e30

      do 13 j=1,nx

         do 11 k=1,ny
            ytmp(k)=za(j,k)
 11      enddo

         call spline(ya,ytmp,ny,yspl,yspl,y2tmp)

         do 12 k=1,ny
            z2a(j,k)=y2tmp(k)
 12      enddo

 13   enddo     

      RETURN
      END



c ... SPLINE INTERPOLATION IN 2D USING PREVIOUS DATA
      SUBROUTINE SPLINT2D(xa,nx,ya,ny,za,z2a,xf,yf,zf)
      INTEGER nx,ny !,NNm
      DOUBLE PRECISION xf,yf,zf,xa(nx),ya(ny),za(nx,ny),z2a(nx,ny)

      INTEGER j,k
      DOUBLE PRECISION yspl,yaux,ytmp(ny),y2tmp(ny),yytmp(nx),x2tmp(nx)

      yspl=5e30

      do 12 j=1,nx

         do 11 k=1,ny
            ytmp(k)=za(j,k)
            y2tmp(k)=z2a(j,k)
 11      enddo

         call splint(ya,ytmp,y2tmp,ny,yf,yaux)
         yytmp(j)=yaux
 12   enddo

      call spline(xa,yytmp,nx,yspl,yspl,x2tmp)
      call splint(xa,yytmp,x2tmp,nx,xf,zf)

      RETURN
      END




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ... POLYNOMIAL INTERPOLATION
c ... THESE INTERPOLATES BY FINDING POLYNOMIALS OF ORDER NPOL-1
c ... THAT REPRODUCE NPOL VALUES OF Y CLOSE TO THE DESIRED X VALUE
c ... THEN IT INTERPOLATES TO THE DESIRED X VALUE 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This routine finds the Npol-1 polynomial that goes through   c
c the Npol closest points to x in the array xa(n). It then     c
c interpolates the polynomial to x, yielding y and the error   c
c dy. The nch variable carries the index of the closer xa      c
c element to x                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE POL_INT(xa,ya,n,Npol,x,y,dy,nchng,ki)

      INTEGER n,nchng,Npol
      DOUBLE PRECISION xa(n),ya(n),x,y,dy

      INTEGER k,ki,khi,klo,Nmin,Nbel,ix,ixa
      DOUBLE PRECISION xb(Npol),yb(Npol)

c ... Choose the number of points above and below x needed
c ... according to the order of the polynomial, Npol
      if(mod(Npol,2).eq.0) then
         Nmin=Npol/2
         Nbel=-Npol/2+1
      else
         Nmin=int(Npol/2)+1
         Nbel=-int(Npol/2)
      endif

c ... Find the index of the two elements xa that bound x
      klo=1
      khi=n

      if(ki.eq.1) go to 1

      if(xa(ki).lt.x.and.xa(ki+1).gt.x) then
         klo=ki
         khi=klo+1
         go to 2
      else if(xa(ki+1).lt.x .and. xa(ki+2).gt.x) then
         klo=ki+1
         khi=klo+1
         go to 2
      endif

  1   if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x) then
            khi=k
         else
            klo=k
         endif
      go to 1
      endif

  2   continue
      ki=klo
c ... Choose the closer index
       nchng=khi
       if( abs(x-xa(khi)) .gt. abs(x-xa(klo)) ) nchng=klo

       if(nchng.lt.Nmin) nchng=Nmin
       if(nchng.gt.n-Nmin) nchng=n-Npol+Nmin

       do ix=1,Npol
          ixa=nchng+Nbel+ix-1
          xb(ix)=xa(ixa)
          yb(ix)=ya(ixa)
       enddo
       call polint(xb,yb,Npol,x,y,dy)

       RETURN 
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Given arrays xa and ya of length n and a given value of x     c
c this routine returns a value y and an error estimate dy of    c
c the polynomial interpolation of x at y.                       c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE POLINT(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER(NMAX=10) 
      
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif=abs(x-xa(1))
      
      do 11 i=1,n
         dift=abs(x-xa(i))
         if(dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
 11   enddo

      y=ya(ns)
      ns=ns-1
      
      do 13 m=1,n-1
         do 12 i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0) then
             write(*,*) 'failure in polint'
             read(*,*)
            endif
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
 12      enddo
         if(2*ns.lt.n-m) then
            dy=c(ns+1)
         else 
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
 13   enddo

      RETURN 
      END


