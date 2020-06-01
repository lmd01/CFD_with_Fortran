program main
implicit none
integer:: i,j,k,l,h
double precision::T,X,pi
integer,parameter::m=160,n=1200!m=20,40,80,160
double precision,parameter::gamma=1.4d0,ep=1d-6
double precision::dx,dt!,error
double precision::v(-6:2*m+6,0:2),p(-6:2*m+6,0:2)
double precision::w(0:2),xk(0:2),u(3,-6:2*m+6,0:2),ua(3,-6:2*m+6,0:n)
double precision::f(3,-1:2*m+1,0:n),f1(3,-1:2*m+1,0:n),f2(3,-1:2*m+1,0:n)
double precision::xx(-4:2*m+6,0:2),up(3,-3:2*m+3,0:n),un(3,-3:2*m+3,0:n)
double precision::up1(3,-3:2*m+3,0:n),un1(3,-3:2*m+3,0:n),up2(3,-3:2*m+3,0:n),un2(3,-3:2*m+3,0:n)
double precision::vn(-1:2*m+3,0:n),vp(-3:2*m+1,0:n),pn(-1:2*m+3,0:n),pp(-3:2*m+1,0:n)
double precision::vn1(-1:2*m+3,0:n),vp1(-3:2*m+1,0:n),pn1(-1:2*m+3,0:n),pp1(-3:2*m+1,0:n)
double precision::vn2(-1:2*m+3,0:n),vp2(-3:2*m+1,0:n),pn2(-1:2*m+3,0:n),pp2(-3:2*m+1,0:n)
double precision::alpha(3,0:1,-2:2*m+2,0:n,0:2),beta(3,-2:2*m+2,0:n,0:2),b(0:1,0:2)
double precision::alpha1(3,0:1,-2:2*m+2,0:2),beta1(3,-2:2*m+2,0:2)
double precision::ur(3,0:1,-3:2*m+3,0:n,0:2),omega(3,0:1,-2:2*m+2,0:n,0:2),omega1(3,0:1,-2:2*m+2,0:2)
double precision::u1(3,-6:2*m+6),u2(3,-6:2*m+6),u1r(3,0:1,-3:2*m+3,0:2)
!------------------------------------------------------
pi=4.d0*datan(1.d0)
X=2.d0
T=0.6d0
dx=X/(2*m)
dt=T/n
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
!-----------------------------------------------
!set the Gauss-Legendre constant
w(2)=0.555555556d0
xk(2)=0.7745966692d0
w(0)=0.555555556d0
xk(0)=-0.7745966692d0
w(1)=0.8888888889d0
xk(1)=0.d0
do i=-4,2*m+6,2
  do k=0,2
   xx(i,k)=dx*xk(k)+i*dx
  enddo
enddo 
!-----------------------------------------------
!set the IC
!u1=/rou u2=/rou*v u3=/E, with E=u3/(gamma-1)+rou*v^2 /2
!-----------------------------------------------------
!left side
v(:,0)=0.d0
do i=-6,m-26,2
  p(i,0)=1.d0
  do k=0,2
    u(1,i,k)=1.d0
	u(2,i,k)=u(1,i,k)*v(i,0)
	u(3,i,k)=p(i,0)/(gamma-1.d0)+u(2,i,k)*v(i,0)/2.d0
  enddo
enddo
!----------------------------------------------------
!right side
do i=m-24,2*m+6,2
  p(i,0)=0.1d0
  do k=0,2
    u(1,i,k)=0.125d0
	u(2,i,k)=u(1,i,k)*v(i,0)
	u(3,i,k)=p(i,0)/(gamma-1.d0)+u(2,i,k)*v(i,0)/2.d0
  enddo
enddo
    
!-----------------------------------------------------
!FV method with Gauss-Legendre quad
do i=-6,2*m+6,2
  do l=1,3
    do k=0,2
     ua(l,i,0)=ua(l,i,0)+w(k)*u(l,i,k)/2
!	 write(*,*)ua(l,i,0),'l=',l
	enddo
  enddo
enddo
!------------------------------------------------------
!time loop
do j=0,n-1
!----------------------------------------------------------
!the value of stencil
  do i=-2,2*m+2,2
    do l=1,3
    ur(l,1,i+1,j,0)=ua(l,i-4,j)/3d0-7d0*ua(l,i-2,j)/6+11d0*ua(l,i,j)/6
    ur(l,1,i+1,j,1)=-ua(l,i-2,j)/6d0+5d0*ua(l,i,j)/6+ua(l,i+2,j)/3d0
    ur(l,1,i+1,j,2)=ua(l,i,j)/3d0+5d0*ua(l,i+2,j)/6-ua(l,i+4,j)/6d0
    ur(l,0,i-1,j,0)=-ua(l,i-4,j)/6d0+5d0*ua(l,i-2,j)/6+ua(l,i,j)/3d0
    ur(l,0,i-1,j,1)=ua(l,i-2,j)/3d0+5d0*ua(l,i,j)/6-ua(l,i+2,j)/6d0
    ur(l,0,i-1,j,2)=11d0*ua(l,i,j)/6-7d0*ua(l,i+2,j)/6+ua(l,i+4,j)/3d0
!--------------------------------------------------------
!the smooth indicators
    beta(l,i,j,0)=(13d0/12)*(ua(l,i-4,j)-2*ua(l,i-2,j)+ua(l,i,j))**2+&
&	(1d0/4)*(ua(l,i-4,j)-4*ua(l,i-2,j)+3*ua(l,i,j))**2
    beta(l,i,j,1)=(13d0/12)*(ua(l,i-2,j)-2*ua(l,i,j)+ua(l,i+2,j))**2+&
&   (1d0/4)*(ua(l,i-2,j)-ua(l,i+2,j))**2
    beta(l,i,j,2)=(13d0/12)*(ua(l,i,j)-2*ua(l,i+2,j)+ua(l,i+4,j))**2+&
&	(1d0/4)*(3*ua(l,i,j)-4*ua(l,i+2,j)+ua(l,i+4,j))**2
!-------------------------------------------------------
!nonlinear weights
    do k=0,1
	  do h=0,2
        alpha(l,k,i,j,h)=b(k,h)/(beta(l,i,j,h)+ep)**2			
	  enddo	
	enddo
	
	do k=0,1
	  do h=0,2
	    omega(l,k,i,j,h)=alpha(l,k,i,j,h)/(alpha(l,k,i,j,0)+alpha(l,k,i,j,1)+&
&       alpha(l,k,i,j,2))
      enddo
	enddo	  
!---------------------------------------------------------------------
!get up&un
    do k=0,2
      un(l,i+1,j)=omega(l,1,i,j,k)*ur(l,1,i+1,j,k)+un(l,i+1,j)
	  up(l,i-1,j)=omega(l,0,i,j,k)*ur(l,0,i-1,j,k)+up(l,i-1,j)
    enddo
  enddo	
  enddo
!---------------------------------------------------------------------------  
  do i=-2,2*m+2,2
    vn(i+1,j)=un(2,i+1,j)/un(1,i+1,j)
	vp(i-1,j)=up(2,i-1,j)/up(1,i-1,j)
	pn(i+1,j)=(gamma-1)*(un(3,i+1,j)-un(2,i+1,j)*vn(i+1,j)/2)
	pp(i-1,j)=(gamma-1)*(up(3,i-1,j)-up(2,i-1,j)*vp(i-1,j)/2)
  enddo
!-------------------------------------------------------
!L-F flux!α=1.2=max|λi|
  do i=-1,2*m+1,2
     f(1,i,j)=(1.d0/2)*(un(2,i,j)+up(2,i,j)-1.2d0*(up(1,i,j)-un(1,i,j)))
     f(2,i,j)=(1.d0/2)*(un(2,i,j)*vn(i,j)+pn(i,j)+up(2,i,j)*vp(i,j)+pp(i,j)&
&	          -1.2d0*(up(2,i,j)-un(2,i,j)))
     f(3,i,j)=(1.d0/2)*((un(3,i,j)+pn(i,j))*vn(i,j)+(up(3,i,j)+pp(i,j))*vp(i,j)&
&	          -1.2d0*(up(3,i,j)-un(3,i,j)))
  enddo
!-------------------------------------------------------
!3rd-order R-K step1
  do i=0,2*m,2 
    do l=1,3
    u1(l,i)=ua(l,i,j)-(dt/2/dx)*(f(l,i+1,j)-f(l,i-1,j))
	enddo
  enddo 
!--------------------------------------------------------
!set the ghost point with BC
  do l=1,3
    u1(l,-2)=ua(l,-2,0)
    u1(l,-4)=ua(l,-4,0)
	u1(l,-6)=ua(l,-6,0)
	u1(l,2*m+2)=ua(l,2*m+2,0)
	u1(l,2*m+4)=ua(l,2*m+4,0)
	u1(l,2*m+6)=ua(l,2*m+6,0)
  enddo
  
  do i=-2,2*m+2,2
    do l=1,3
    u1r(l,1,i+1,0)=u1(l,i-4)/3d0-7d0*u1(l,i-2)/6+11d0*u1(l,i)/6
    u1r(l,1,i+1,1)=-u1(l,i-2)/6d0+5d0*u1(l,i)/6+u1(l,i+2)/3d0
    u1r(l,1,i+1,2)=u1(l,i)/3d0+5d0*u1(l,i+2)/6-u1(l,i+4)/6d0
    u1r(l,0,i-1,0)=-u1(l,i-4)/6d0+5d0*u1(l,i-2)/6+u1(l,i)/3d0
    u1r(l,0,i-1,1)=u1(l,i-2)/3d0+5d0*u1(l,i)/6-u1(l,i+2)/6d0
    u1r(l,0,i-1,2)=11d0*u1(l,i)/6-7d0*u1(l,i+2)/6+u1(l,i+4)/3d0
 
    beta1(l,i,0)=(13d0/12)*(u1(l,i-4)-2*u1(l,i-2)+u1(l,i))**2+&
&	(1d0/4)*(u1(l,i-4)-4*u1(l,i-2)+3*u1(l,i))**2
    beta1(l,i,1)=(13d0/12)*(u1(l,i-2)-2*u1(l,i)+u1(l,i+2))**2+&
&   (1d0/4)*(u1(l,i-2)-u1(l,i+2))**2
    beta1(l,i,2)=(13d0/12)*(u1(l,i)-2*u1(l,i+2)+u1(l,i+4))**2+&
&	(1d0/4)*(3*u1(l,i)-4*u1(l,i+2)+u1(l,i+4))**2

    do k=0,1
	  do h=0,2
        alpha1(l,k,i,h)=b(k,h)/(beta1(l,i,h)+ep)**2			
	  enddo	
	enddo
	
	do k=0,1
	  do h=0,2
	    omega1(l,k,i,h)=alpha1(l,k,i,h)/(alpha1(l,k,i,0)+alpha1(l,k,i,1)+&
&       alpha1(l,k,i,2))
      enddo
	enddo	  
!---------------------------------------------------------------------
!get up1&un1
    do k=0,2
      un1(l,i+1,j)=omega1(l,1,i,k)*u1r(l,1,i+1,k)+un1(l,i+1,j)
	  up1(l,i-1,j)=omega1(l,0,i,k)*u1r(l,0,i-1,k)+up1(l,i-1,j)
    enddo
  enddo	
  enddo
  
  do i=-2,2*m+2,2
    vn1(i+1,j)=un1(2,i+1,j)/un1(1,i+1,j)
	vp1(i-1,j)=up1(2,i-1,j)/up1(1,i-1,j)
	pn1(i+1,j)=(gamma-1)*(un1(3,i+1,j)-un1(2,i+1,j)*vn1(i+1,j)/2)
	pp1(i-1,j)=(gamma-1)*(up1(3,i-1,j)-up1(2,i-1,j)*vp1(i-1,j)/2)
  enddo
  
  do i=-1,2*m+1,2
     f1(1,i,j)=(1.d0/2)*(un1(2,i,j)+up1(2,i,j)-1.2d0*(up1(1,i,j)-un1(1,i,j)))
     f1(2,i,j)=(1.d0/2)*(un1(2,i,j)*vn1(i,j)+pn1(i,j)+up1(2,i,j)*vp1(i,j)+pp1(i,j)&
&	          -1.2d0*(up1(2,i,j)-un1(2,i,j)))
     f1(3,i,j)=(1.d0/2)*((un1(3,i,j)+pn1(i,j))*vn1(i,j)+(up1(3,i,j)+pp1(i,j))*vp1(i,j)&
&	          -1.2d0*(up1(3,i,j)-un1(3,i,j)))
  enddo
!-------------------------------------------------------------------------
!step2  
  do i=0,2*m,2 
    do l=1,3
	  u2(l,i)=3*ua(l,i,j)/4+(u1(l,i)-(dt/(2*dx))*(f1(l,i+1,j)-f1(l,i-1,j)))/4
	enddo
  enddo
  
  do l=1,3
    u2(l,-2)=ua(l,-2,0)
    u2(l,-4)=ua(l,-4,0)
	u2(l,-6)=ua(l,-6,0)
	u2(l,2*m+2)=ua(l,2*m+2,0)
	u2(l,2*m+4)=ua(l,2*m+4,0)
	u2(l,2*m+6)=ua(l,2*m+6,0)
  enddo
  
  do i=-2,2*m+2,2
    do l=1,3
    u1r(l,1,i+1,0)=u2(l,i-4)/3d0-7d0*u2(l,i-2)/6+11d0*u2(l,i)/6
    u1r(l,1,i+1,1)=-u2(l,i-2)/6d0+5d0*u2(l,i)/6+u2(l,i+2)/3d0
    u1r(l,1,i+1,2)=u2(l,i)/3d0+5d0*u2(l,i+2)/6-u2(l,i+4)/6d0
    u1r(l,0,i-1,0)=-u2(l,i-4)/6d0+5d0*u2(l,i-2)/6+u2(l,i)/3d0
    u1r(l,0,i-1,1)=u2(l,i-2)/3d0+5d0*u2(l,i)/6-u2(l,i+2)/6d0
    u1r(l,0,i-1,2)=11d0*u2(l,i)/6-7d0*u2(l,i+2)/6+u2(l,i+4)/3d0
 
    beta1(l,i,0)=(13d0/12)*(u2(l,i-4)-2*u2(l,i-2)+u2(l,i))**2+&
&	(1d0/4)*(u2(l,i-4)-4*u2(l,i-2)+3*u2(l,i))**2
    beta1(l,i,1)=(13d0/12)*(u2(l,i-2)-2*u2(l,i)+u2(l,i+2))**2+&
&   (1d0/4)*(u2(l,i-2)-u2(l,i+2))**2
    beta1(l,i,2)=(13d0/12)*(u2(l,i)-2*u2(l,i+2)+u2(l,i+4))**2+&
&	(1d0/4)*(3*u2(l,i)-4*u2(l,i+2)+u2(l,i+4))**2

    do k=0,1
	  do h=0,2
        alpha1(l,k,i,h)=b(k,h)/(beta1(l,i,h)+ep)**2			
	  enddo	
	enddo
	
	do k=0,1
	  do h=0,2
	    omega1(l,k,i,h)=alpha1(l,k,i,h)/(alpha1(l,k,i,0)+alpha1(l,k,i,1)+&
&       alpha1(l,k,i,2))
      enddo
	enddo	  
!---------------------------------------------------------------------
!get up2&un2
    do k=0,2
      un2(l,i+1,j)=omega1(l,1,i,k)*u1r(l,1,i+1,k)+un2(l,i+1,j)
	  up2(l,i-1,j)=omega1(l,0,i,k)*u1r(l,0,i-1,k)+up2(l,i-1,j)
    enddo
  enddo	
  enddo
  
  do i=-2,2*m+2,2
    vn2(i+1,j)=un2(2,i+1,j)/un2(1,i+1,j)
	vp2(i-1,j)=up2(2,i-1,j)/up2(1,i-1,j)
	pn2(i+1,j)=(gamma-1)*(un2(3,i+1,j)-un2(2,i+1,j)*vn2(i+1,j)/2)
	pp2(i-1,j)=(gamma-1)*(up2(3,i-1,j)-up2(2,i-1,j)*vp2(i-1,j)/2)
  enddo
  
  do i=-1,2*m+1,2
     f2(1,i,j)=(1.d0/2)*(un2(2,i,j)+up2(2,i,j)-1.2d0*(up2(1,i,j)-un2(1,i,j)))
     f2(2,i,j)=(1.d0/2)*(un2(2,i,j)*vn2(i,j)+pn2(i,j)+up2(2,i,j)*vp2(i,j)+pp2(i,j)&
&	          -1.2d0*(up2(2,i,j)-un2(2,i,j)))
     f2(3,i,j)=(1.d0/2)*((un2(3,i,j)+pn2(i,j))*vn2(i,j)+(up2(3,i,j)+pp2(i,j))*vp2(i,j)&
&	          -1.2d0*(up2(3,i,j)-un2(3,i,j)))
  enddo
  
  do i=0,2*m,2 
    do l=1,3  
	  ua(l,i,j+1)=ua(l,i,j)/3+2*(u2(l,i)-(dt/(2*dx))*(f2(l,i+1,j)-f2(l,i-1,j)))/3
    enddo
  enddo
  
  do l=1,3
    ua(l,-2,j+1)=ua(l,-2,0)
    ua(l,-4,j+1)=ua(l,-4,0)
	ua(l,-6,j+1)=ua(l,-6,0)
	ua(l,2*m+2,j+1)=ua(l,2*m+2,0)
	ua(l,2*m+4,j+1)=ua(l,2*m+4,0)
	ua(l,2*m+6,j+1)=ua(l,2*m+6,0)
  enddo
  
enddo
!-------------------------------------------------------------
!output needed to be trans to rou v and p form
open(111,file='1DEulerm160shockrho0.txt')
do i=0,2*m,2
!  do l=1,3
  write(111,*) ua(1,i,0)!(gamma-1)*(ua(3,i,n)-ua(2,i,n)**2/2/ua(1,i,n))!ua(1,i,n)!ua(2,i,n)/ua(1,i,n)
!  enddo
enddo
end 
