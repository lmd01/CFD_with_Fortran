program main
implicit none
integer:: i,j,k,l
double precision::T,X,pi
integer,parameter::m=320,n=390!10,20,40,80,160,320;n is given by T/cfl/dx^(5/3)
double precision,parameter::ep=1d-6
double precision::dx,dt
double precision::u(-7:2*m+7,0:n),f(-1:2*m+1,0:n)
double precision::f1(-1:2*m+1,0:n),f2(-1:2*m+1,0:n),u1(-6:2*m+6),u2(-6:2*m+6)
double precision::ua(-6:2*m+6,0:n),ur(0:1,-3:2*m+3,0:n,0:2)
double precision::alpha(0:1,-2:2*m+2,0:n,0:2),beta(-2:2*m+2,0:n,0:2),b(0:1,0:2)
double precision::omega(0:1,-2:2*m+2,0:n,0:2),omega1(0:1,-2:2*m+2,0:2)
double precision::up(-3:2*m+3,0:n),un(-3:2*m+3,0:n)
double precision::up1(-3:2*m+3,0:n),un1(-3:2*m+3,0:n),up2(-3:2*m+3,0:n),un2(-3:2*m+3,0:n)
double precision::u1r(0:1,-3:2*m+3,0:2),beta1(-2:2*m+2,0:2),alpha1(0:1,-2:2*m+2,0:2)
pi=4.d0*datan(1.d0)
X=2.d0*pi
T=0.5d0!1.5d0,0.5d0
dx=X/(2*m)
dt=T/n
!write(*,*) ep
!-----------------------------------------------
!IC
do i=-7,2*m+7
  u(i,0)=0.5d0+dsin(i*dx)
enddo 
!FV method
do i=-6,2*m+6,2
  ua(i,0)=(u(i+1,0)+2d0*u(i,0)+u(i-1,0))/4d0
enddo
	
!-----------------------------------------------
!5th order WENO, with k=2
!------------------------------------------------
!linear weights
b(1,0)=1.d0/10
b(1,1)=3.d0/5
b(1,2)=3.d0/10
b(0,0)=3.d0/10!dr(1,2)
b(0,1)=3.d0/5!dr(1,1)
b(0,2)=1.d0/10!dr(1,0)
!--------------------------------------------------
!time loop
do j=0,n-1
!----------------------------------------------------------
!the value of stencil
  do i=-2,2*m+2,2
    ur(1,i+1,j,0)=ua(i-4,j)/3d0-7d0*ua(i-2,j)/6+11d0*ua(i,j)/6
    ur(1,i+1,j,1)=-ua(i-2,j)/6d0+5d0*ua(i,j)/6+ua(i+2,j)/3d0
    ur(1,i+1,j,2)=ua(i,j)/3d0+5d0*ua(i+2,j)/6-ua(i+4,j)/6d0
    ur(0,i-1,j,0)=-ua(i-4,j)/6d0+5d0*ua(i-2,j)/6+ua(i,j)/3d0
    ur(0,i-1,j,1)=ua(i-2,j)/3d0+5d0*ua(i,j)/6-ua(i+2,j)/6d0
    ur(0,i-1,j,2)=11d0*ua(i,j)/6-7d0*ua(i+2,j)/6+ua(i+4,j)/3d0
!--------------------------------------------------------
!the smooth indicators
    beta(i,j,0)=(13d0/12)*(ua(i-4,j)-2*ua(i-2,j)+ua(i,j))**2+&
&	(1d0/4)*(ua(i-4,j)-4*ua(i-2,j)+3*ua(i,j))**2
    beta(i,j,1)=(13d0/12)*(ua(i-2,j)-2*ua(i,j)+ua(i+2,j))**2+&
&   (1d0/4)*(ua(i-2,j)-ua(i+2,j))**2
    beta(i,j,2)=(13d0/12)*(ua(i,j)-2*ua(i+2,j)+ua(i+4,j))**2+&
&	(1d0/4)*(3*ua(i,j)-4*ua(i+2,j)+ua(i+4,j))**2
!-------------------------------------------------------
!nonlinear weights
    do k=0,1
	  do l=0,2
        alpha(k,i,j,l)=b(k,l)/(beta(i,j,l)+ep)**2			
	  enddo	
	enddo
	
	do k=0,1
	  do l=0,2
	    omega(k,i,j,l)=alpha(k,i,j,l)/(alpha(k,i,j,0)+alpha(k,i,j,1)+&
&       alpha(k,i,j,2))
      enddo
	enddo	  
!---------------------------------------------------------------------
!get up&un
    do l=0,2
      un(i+1,j)=omega(1,i,j,l)*ur(1,i+1,j,l)+un(i+1,j)
	  up(i-1,j)=omega(0,i,j,l)*ur(0,i-1,j,l)+up(i-1,j)
    enddo
  enddo
!-----------------------------------------------------------------------
!get the numerical flux
  do i=-1,2*m+1,2
    if(un(i,j).lt.up(i,j)) then
      f(i,j)=min((un(i,j)**2)/2,(up(i,j)**2)/2)
	else
	  f(i,j)=max((un(i,j)**2)/2,(up(i,j)**2)/2)
    endif  
  enddo

!---------------------------------------------------------
!3rd-order R-K
    do i=0,2*m,2 !1st step
	  u1(i)=ua(i,j)-(dt/(2*dx))*(f(i+1,j)-f(i-1,j)) 
	enddo
	
	do k=-6,-2,2 
    u1(k)=u1(2*m+k)
    enddo
    do k=2*m+2,2*m+6,2
    u1(k)=u1(k-2*m)
    enddo
	
  do i=-2,2*m+2,2
    u1r(1,i+1,0)=u1(i-4)/3d0-7d0*u1(i-2)/6+11d0*u1(i)/6
    u1r(1,i+1,1)=-u1(i-2)/6d0+5d0*u1(i)/6+u1(i+2)/3d0
    u1r(1,i+1,2)=u1(i)/3d0+5d0*u1(i+2)/6-u1(i+4)/6d0
    u1r(0,i-1,0)=-u1(i-4)/6d0+5d0*u1(i-2)/6+u1(i)/3d0
    u1r(0,i-1,1)=u1(i-2)/3d0+5d0*u1(i)/6-u1(i+2)/6d0
    u1r(0,i-1,2)=11d0*u1(i)/6-7d0*u1(i+2)/6+u1(i+4)/3d0
	
    beta1(i,0)=(13d0/12)*(u1(i-4)-2*u1(i-2)+u1(i))**2+&
&	(1d0/4)*(u1(i-4)-4*u1(i-2)+3*u1(i))**2
    beta1(i,1)=(13d0/12)*(u1(i-2)-2*u1(i)+u1(i+2))**2+&
&   (1d0/4)*(u1(i-2)-u1(i+2))**2
    beta1(i,2)=(13d0/12)*(u1(i)-2*u1(i+2)+u1(i+4))**2+&
&	(1d0/4)*(3*u1(i)-4*u1(i+2)+u1(i+4))**2
 
    do k=0,1
	  do l=0,2
        alpha1(k,i,l)=b(k,l)/(beta1(i,l)+ep)**2	
	  enddo	
	enddo
	
	do k=0,1
	  do l=0,2
	    omega1(k,i,l)=alpha1(k,i,l)/(alpha1(k,i,0)+alpha1(k,i,1)+&
&       alpha1(k,i,2))
      enddo
	enddo	

    do l=0,2
      un1(i+1,j)=omega1(1,i,l)*u1r(1,i+1,l)+un1(i+1,j)
	  up1(i-1,j)=omega1(0,i,l)*u1r(0,i-1,l)+up1(i-1,j)
    enddo	
  enddo
!---------------------------------------------------------	
	do k=-1,2*m+1,2! 2nd step
       if(un1(k,j).lt.up1(k,j)) then
        f1(k,j)=min((un1(k,j)**2)/2,(up1(k,j)**2)/2)
	   else
	    f1(k,j)=max((un1(k,j)**2)/2,(up1(k,j)**2)/2)
       endif  
    enddo 
	
	do i=0,2*m,2  
	  u2(i)=3*ua(i,j)/4+(u1(i)-(dt/(2*dx))*(f1(i+1,j)-f1(i-1,j)))/4
	enddo
	
	do k=-6,-2,2 
    u2(k)=u2(2*m+k)
    enddo
    do k=2*m+2,2*m+6,2
    u2(k)=u2(k-2*m)
    enddo
	
  do i=-2,2*m+2,2
    u1r(1,i+1,0)=u2(i-4)/3d0-7d0*u2(i-2)/6+11d0*u2(i)/6
    u1r(1,i+1,1)=-u2(i-2)/6d0+5d0*u2(i)/6+u2(i+2)/3d0
    u1r(1,i+1,2)=u2(i)/3d0+5d0*u2(i+2)/6-u2(i+4)/6d0
    u1r(0,i-1,0)=-u2(i-4)/6d0+5d0*u2(i-2)/6+u2(i)/3d0
    u1r(0,i-1,1)=u2(i-2)/3d0+5d0*u2(i)/6-u2(i+2)/6d0
    u1r(0,i-1,2)=11d0*u2(i)/6-7d0*u2(i+2)/6+u2(i+4)/3d0
	
    beta1(i,0)=(13d0/12)*(u2(i-4)-2*u2(i-2)+u2(i))**2+&
&	(1d0/4)*(u2(i-4)-4*u2(i-2)+3*u2(i))**2
    beta1(i,1)=(13d0/12)*(u2(i-2)-2*u2(i)+u2(i+2))**2+&
&   (1d0/4)*(u2(i-2)-u2(i+2))**2
    beta1(i,2)=(13d0/12)*(u2(i)-2*u2(i+2)+u2(i+4))**2+&
&	(1d0/4)*(3*u2(i)-4*u2(i+2)+u2(i+4))**2
 
    do k=0,1
	  do l=0,2
        alpha1(k,i,l)=b(k,l)/(beta1(i,l)+ep)**2		
	  enddo	
	enddo
	
	do k=0,1
	  do l=0,2
	    omega1(k,i,l)=alpha1(k,i,l)/(alpha1(k,i,0)+alpha1(k,i,1)+&
&       alpha1(k,i,2))
      enddo
	enddo	

    do l=0,2
      un2(i+1,j)=omega1(1,i,l)*u1r(1,i+1,l)+un2(i+1,j)
	  up2(i-1,j)=omega1(0,i,l)*u1r(0,i-1,l)+up2(i-1,j)
    enddo	
  enddo
!---------------------------------------------------------------	
	do k=-1,2*m+1,2! 3rd step
       if(un2(k,j).lt.up2(k,j)) then
        f2(k,j)=min((un2(k,j)**2)/2,(up2(k,j)**2)/2)
	   else
	    f2(k,j)=max((un2(k,j)**2)/2,(up2(k,j)**2)/2)
       endif  
    enddo
	
	do i=0,2*m,2  
	  ua(i,j+1)=ua(i,j)/3+2*(u2(i)-(dt/(2*dx))*(f2(i+1,j)-f2(i-1,j)))/3
	enddo
!-------------------------------------------------------------------
!reset the ghost point to fit the period BC
  do k=-6,-2,2 
    ua(k,j+1)=ua(2*m+k,j+1)
  enddo
  do k=2*m+2,2*m+6,2
    ua(k,j+1)=ua(k-2*m,j+1)
  enddo
enddo 

!do i=-6,2*m+6,2
!j=1
!  write(*,*) ua(i,j),i,j
!enddo
!-----------------------------------------------
!output
open(111,file='WENO_m320.txt')
do i=0,2*m,2
  write(111,*) ua(i,n)
enddo
end
