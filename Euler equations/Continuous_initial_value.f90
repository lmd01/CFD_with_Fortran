program main
implicit none
integer:: i,j,k,l
double precision::T,X,pi
integer,parameter::m=160,n=32!m=20,40,80,160
double precision,parameter::gamma=1.4d0
double precision::dx,dt,error
double precision::v(-4:2*m+6,0:2),p(-4:2*m+6,0:2)
double precision::w(0:2),xk(0:2),u(3,-4:2*m+6,0:2),ua(3,-4:2*m+6,0:n),f(3,-1:2*m+1,0:n)
double precision::u1(3,-4:2*m+6),u2(3,-4:2*m+6),f1(3,-1:2*m+1,0:n)
double precision::xx(-4:2*m+6,0:2),uap(3,-3:2*m+1,0:n),uan(3,-1:2*m+3,0:n)
double precision::uap1(3,-3:2*m+1,0:n),uan1(3,-1:2*m+3,0:n)
double precision::vn(-1:2*m+3,0:n),vp(-3:2*m+1,0:n),pn(-1:2*m+3,0:n),pp(-3:2*m+1,0:n)
double precision::vn1(-1:2*m+3,0:n),vp1(-3:2*m+1,0:n),pn1(-1:2*m+3,0:n),pp1(-3:2*m+1,0:n)
double precision::uac(0:2*m,0:2),uacc(0:2*m)
!------------------------------------------------------
pi=4.d0*datan(1.d0)
X=2.d0*pi
T=0.6d0
dx=X/(2*m)
dt=T/n
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
v(:,0)=1.d0
p(:,0)=1.d0
do i=-4,2*m+6,2
  do k=0,2
    u(1,i,k)=1.d0+0.2*dsin(xx(i,k))!accurate:dsin(xx-T)
	u(2,i,k)=u(1,i,k)*v(i,0)
	u(3,i,k)=p(i,0)/(gamma-1.d0)+u(2,i,k)*v(i,0)/2.d0
  enddo
enddo
!-----------------------------------------------------
!FV method with Gauss-Legendre quad
do i=-4,2*m+6,2
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
  do i=-2,2*m+2,2
    do l=1,3
    uan(l,i+1,j)=-ua(l,i-2,j)/6+5*ua(l,i,j)/6+ua(l,i+2,j)/3
	uap(l,i-1,j)=ua(l,i,j)/3+5*ua(l,i+2,j)/6-ua(l,i+4,j)/6
	enddo
  enddo
  
  do i=-2,2*m+2,2
    vn(i+1,j)=uan(2,i+1,j)/uan(1,i+1,j)
	vp(i-1,j)=uap(2,i-1,j)/uap(1,i-1,j)
	pn(i+1,j)=(gamma-1)*(uan(3,i+1,j)-uan(2,i+1,j)*vn(i+1,j)/2)
	pp(i-1,j)=(gamma-1)*(uap(3,i-1,j)-uap(2,i-1,j)*vp(i-1,j)/2)
  enddo
!-------------------------------------------------------
!L-F flux!α=1.2=max|λi|
  do i=-1,2*m+1,2
     f(1,i,j)=(1.d0/2)*(uan1(2,i,j)+uap1(2,i,j)-1.2d0*(uap1(1,i,j)-uan1(1,i,j)))
     f(2,i,j)=(1.d0/2)*(uan1(2,i,j)*vn1(i,j)+pn1(i,j)+uap1(2,i,j)*vp1(i,j)+pp1(i,j)&
&	          -1.2d0*(uap1(2,i,j)-uan1(2,i,j)))
     f(3,i,j)=(1.d0/2)*((uan1(3,i,j)+pn1(i,j))*vn1(i,j)+(uap1(3,i,j)+pp1(i,j))*vp1(i,j)&
&	          -1.2d0*(uap1(3,i,j)-uan1(3,i,j)))
  enddo
!-------------------------------------------------------
!3rd order R-K,step 1
  do i=0,2*m,2 
    do l=1,3
    u1(l,i)=ua(l,i,j)-(dt/2/dx)*(f(l,i+1,j)-f(l,i-1,j))
!	write(*,*)ua(2,i,j)/ua(1,i,j),'l=',l,'i=',i,'j=',j
	enddo
  enddo
!---------------------------------------------------------
!set the ghost point with BC
  do l=1,3
    u1(l,-2)=u1(l,2*m-2)
    u1(l,-4)=u1(l,2*m-4)
	u1(l,2*m+2)=u1(l,2)
	u1(l,2*m+4)=u1(l,4)
	u1(l,2*m+6)=u1(l,6)
  enddo
  
  do i=-2,2*m+2,2
    do l=1,3
    uan1(l,i+1,j)=-u1(l,i-2)/6+5*u1(l,i)/6+u1(l,i+2)/3
	uap1(l,i-1,j)=u1(l,i)/3+5*u1(l,i+2)/6-u1(l,i+4)/6
	enddo
  enddo
  
  do i=-2,2*m+2,2
    vn1(i+1,j)=uan1(2,i+1,j)/uan1(1,i+1,j)
	vp1(i-1,j)=uap1(2,i-1,j)/uap1(1,i-1,j)
	pn1(i+1,j)=(gamma-1)*(uan1(3,i+1,j)-uan1(2,i+1,j)*vn1(i+1,j)/2)
	pp1(i-1,j)=(gamma-1)*(uap1(3,i-1,j)-uap1(2,i-1,j)*vp1(i-1,j)/2)
  enddo
  
  do i=-1,2*m+1,2
     f1(1,i,j)=(1.d0/2)*(uan1(2,i,j)+uap1(2,i,j)-1.2d0*(uap1(1,i,j)-uan1(1,i,j)))
     f1(2,i,j)=(1.d0/2)*(uan1(2,i,j)*vn1(i,j)+pn1(i,j)+uap1(2,i,j)*vp1(i,j)+pp1(i,j)&
&	          -1.2d0*(uap1(2,i,j)-uan1(2,i,j)))
     f1(3,i,j)=(1.d0/2)*((uan1(3,i,j)+pn1(i,j))*vn1(i,j)+(uap1(3,i,j)+pp1(i,j))*vp1(i,j)&
&	          -1.2d0*(uap1(3,i,j)-uan1(3,i,j)))
  enddo
  
!-------------------------------------------------------------------------
!step2  
  do i=0,2*m,2 
    do l=1,3
	  u2(l,i)=3*ua(l,i,j)/4+(u1(l,i)-(dt/(2*dx))*(f1(l,i+1,j)-f1(l,i-1,j)))/4
	enddo
  enddo
  
  do l=1,3
    u2(l,-2)=u2(l,2*m-2)
    u2(l,-4)=u2(l,2*m-4)
	u2(l,2*m+2)=u2(l,2)
	u2(l,2*m+4)=u2(l,4)
	u2(l,2*m+6)=u2(l,6)
  enddo
  
  do i=-2,2*m+2,2
    do l=1,3
    uan1(l,i+1,j)=-u2(l,i-2)/6+5*u2(l,i)/6+u2(l,i+2)/3
	uap1(l,i-1,j)=u2(l,i)/3+5*u2(l,i+2)/6-u2(l,i+4)/6
	enddo
  enddo
  
  do i=-2,2*m+2,2
    vn1(i+1,j)=uan1(2,i+1,j)/uan1(1,i+1,j)
	vp1(i-1,j)=uap1(2,i-1,j)/uap1(1,i-1,j)
	pn1(i+1,j)=(gamma-1)*(uan1(3,i+1,j)-uan1(2,i+1,j)*vn1(i+1,j)/2)
	pp1(i-1,j)=(gamma-1)*(uap1(3,i-1,j)-uap1(2,i-1,j)*vp1(i-1,j)/2)
  enddo
  
  do i=-1,2*m+1,2
     f1(1,i,j)=(1.d0/2)*(uan1(2,i,j)+uap1(2,i,j)-1.2d0*(uap1(1,i,j)-uan1(1,i,j)))
     f1(2,i,j)=(1.d0/2)*(uan1(2,i,j)*vn1(i,j)+pn1(i,j)+uap1(2,i,j)*vp1(i,j)+pp1(i,j)&
&	          -1.2d0*(uap1(2,i,j)-uan1(2,i,j)))
     f1(3,i,j)=(1.d0/2)*((uan1(3,i,j)+pn1(i,j))*vn1(i,j)+(uap1(3,i,j)+pp1(i,j))*vp1(i,j)&
&	          -1.2d0*(uap1(3,i,j)-uan1(3,i,j)))
  enddo
  
  do i=0,2*m,2 
    do l=1,3  
	  ua(l,i,j+1)=ua(l,i,j)/3+2*(u2(l,i)-(dt/(2*dx))*(f1(l,i+1,j)-f1(l,i-1,j)))/3
    enddo
  enddo
  
  do l=1,3
    ua(l,-2,j+1)=ua(l,2*m-2,j+1)
    ua(l,-4,j+1)=ua(l,2*m-4,j+1)
	ua(l,2*m+2,j+1)=ua(l,2,j+1)
	ua(l,2*m+4,j+1)=ua(l,4,j+1)
	ua(l,2*m+6,j+1)=ua(l,6,j+1)
  enddo
  
enddo
!-------------------------------------------------------------
!output needed to be trans to rou v and p form
open(111,file='1DEulerm160.txt')
do i=0,2*m,2
!  do l=1,3
  write(111,*) ua(1,i,n)
!  enddo
enddo
!----------------------------------------------------------------
!get error
do i=0,2*m,2
  do k=0,2
  uac(i,k)=1.d0+0.2*dsin(xx(i,k)-T)
  enddo
enddo
!-----------------------------------------------------------
!cell average of accurate value
do i=0,2*m,2
  do k=0,2
  uacc(i)=uacc(i)+w(k)*uac(i,k)/2
  enddo
enddo
  
do i=0,2*m,2
  error=error+dabs(uacc(i)-ua(1,i,n))*(2*dx)  
enddo

write(*,*) error

end 
