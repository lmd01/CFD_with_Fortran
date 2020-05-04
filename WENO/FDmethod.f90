program main
implicit none
integer:: i,j,k,l
double precision::T,X,pi
integer,parameter::m=80,n=120!10,20,40,80,160,320;
double precision,parameter::ep=1d-6
double precision::dx,dt
double precision::u(-7:2*m+7,0:n),f(-6:2*m+6,0:n),u1(-7:2*m+7),u2(-7:2*m+7)
double precision::f1(-6:2*m+6),f2(-6:2*m+6),h1(-6:2*m+6,0:1),h2(-6:2*m+6,0:1)
double precision::fp(-6:2*m+6,0:n),fn(-6:2*m+6,0:n)
double precision::fp1(-6:2*m+6),fn1(-6:2*m+6),fp2(-6:2*m+6),fn2(-6:2*m+6)
double precision::h(-6:2*m+6,0:n,0:1),hr(0:1,-3:2*m+3,0:n,0:2)
double precision::alpha(0:1,-2:2*m+2,0:n,0:2),beta(0:1,-2:2*m+2,0:n,0:2),b(0:1,0:2)
double precision::omega(0:1,-2:2*m+2,0:n,0:2),omega1(0:1,-2:2*m+2,0:2)
double precision::hp(-3:2*m+3,0:n),hn(-3:2*m+3,0:n)
double precision::hp1(-3:2*m+3,0:n),hn1(-3:2*m+3,0:n),hp2(-3:2*m+3,0:n),hn2(-3:2*m+3,0:n)
double precision::h1r(0:1,-3:2*m+3,0:2),beta1(0:1,-2:2*m+2,0:2),alpha1(0:1,-2:2*m+2,0:2)
pi=4.d0*datan(1.d0)
X=2.d0*pi
T=1.5d0!1.5d0,0.5d0
dx=X/(2*m)
dt=T/n
!write(*,*) ep
!-----------------------------------------------
!IC
do i=-7,2*m+7
  u(i,0)=0.5d0+dsin(i*dx)
enddo 
!FD method
do i=-6,2*m+6,2
  f(i,0)=u(i,0)**2/2
enddo

do i=-6,2*m+6,2
  fp(i,0)=(f(i,0)+1.5d0*u(i,0))/2
  fn(i,0)=(f(i,0)-1.5d0*u(i,0))/2
enddo

do i=-6,2*m+6,2
  h(i,0,0)=fp(i,0)!0 means minus
  h(i,0,1)=fn(i,0)!1 means plus
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
    hr(1,i+1,j,0)=h(i-4,j,0)/3d0-7d0*h(i-2,j,0)/6+11d0*h(i,j,0)/6
    hr(1,i+1,j,1)=-h(i-2,j,0)/6d0+5d0*h(i,j,0)/6+h(i+2,j,0)/3d0
    hr(1,i+1,j,2)=h(i,j,0)/3d0+5d0*h(i+2,j,0)/6-h(i+4,j,0)/6d0
    hr(0,i-1,j,0)=-h(i-4,j,1)/6d0+5d0*h(i-2,j,1)/6+h(i,j,1)/3d0
    hr(0,i-1,j,1)=h(i-2,j,1)/3d0+5d0*h(i,j,1)/6-h(i+2,j,1)/6d0
    hr(0,i-1,j,2)=11d0*h(i,j,1)/6-7d0*h(i+2,j,1)/6+h(i+4,j,1)/3d0
!--------------------------------------------------------
!the smooth indicators
    beta(1,i,j,0)=(13d0/12)*(h(i-4,j,0)-2*h(i-2,j,0)+h(i,j,0))**2+&
&	(1d0/4)*(h(i-4,j,0)-4*h(i-2,j,0)+3*h(i,j,0))**2
    beta(1,i,j,1)=(13d0/12)*(h(i-2,j,0)-2*h(i,j,0)+h(i+2,j,0))**2+&
&   (1d0/4)*(h(i-2,j,0)-h(i+2,j,0))**2
    beta(1,i,j,2)=(13d0/12)*(h(i,j,0)-2*h(i+2,j,0)+h(i+4,j,0))**2+&
&	(1d0/4)*(3*h(i,j,0)-4*h(i+2,j,0)+h(i+4,j,0))**2

    beta(0,i,j,0)=(13d0/12)*(h(i-4,j,1)-2*h(i-2,j,1)+h(i,j,1))**2+&
&	(1d0/4)*(h(i-4,j,1)-4*h(i-2,j,1)+3*h(i,j,1))**2
    beta(0,i,j,1)=(13d0/12)*(h(i-2,j,1)-2*h(i,j,1)+h(i+2,j,1))**2+&
&   (1d0/4)*(h(i-2,j,1)-h(i+2,j,1))**2
    beta(0,i,j,2)=(13d0/12)*(h(i,j,1)-2*h(i+2,j,1)+h(i+4,j,1))**2+&
&	(1d0/4)*(3*h(i,j,1)-4*h(i+2,j,1)+h(i+4,j,1))**2
!-------------------------------------------------------
!nonlinear weights
    do k=0,1
	  do l=0,2
        alpha(k,i,j,l)=b(k,l)/(beta(k,i,j,l)+ep)**2		
	  enddo	
	enddo
	
	do k=0,1
	  do l=0,2
	    omega(k,i,j,l)=alpha(k,i,j,l)/(alpha(k,i,j,0)+alpha(k,i,j,1)+&
&       alpha(k,i,j,2))
      enddo
	enddo	  
!---------------------------------------------------------------------
!get hp&hn
    do l=0,2
      hn(i+1,j)=omega(1,i,j,l)*hr(1,i+1,j,l)+hn(i+1,j)
	  hp(i-1,j)=omega(0,i,j,l)*hr(0,i-1,j,l)+hp(i-1,j)
    enddo
  enddo
!-----------------------------------------------------------------------
!get the numerical flux
  do i=-1,2*m+1,2
    f(i,j)=hn(i,j)+hp(i,j)
  enddo

!---------------------------------------------------------
!3rd-order R-K
    do i=0,2*m,2 !1st step
	  u1(i)=u(i,j)-(dt/(2*dx))*(f(i+1,j)-f(i-1,j)) 
	enddo
	
	do k=-6,-2,2 
    u1(k)=u1(2*m+k)
    enddo
    do k=2*m+2,2*m+6,2
    u1(k)=u1(k-2*m)
    enddo
	
	do i=-6,2*m+6,2
     f1(i)=u1(i)**2/2
    enddo

    do i=-6,2*m+6,2
      fp1(i)=(f1(i)+1.5d0*u1(i))/2
      fn1(i)=(f1(i)-1.5d0*u1(i))/2
    enddo

    do i=-6,2*m+6,2
      h1(i,0)=fp1(i)!0 means minus
      h1(i,1)=fn1(i)!1 means plus
    enddo
	
  do i=-2,2*m+2,2
    h1r(1,i+1,0)=h1(i-4,0)/3d0-7d0*h1(i-2,0)/6+11d0*h1(i,0)/6
    h1r(1,i+1,1)=-h1(i-2,0)/6d0+5d0*h1(i,0)/6+h1(i+2,0)/3d0
    h1r(1,i+1,2)=h1(i,0)/3d0+5d0*h1(i+2,0)/6-h1(i+4,0)/6d0
    h1r(0,i-1,0)=-h1(i-4,1)/6d0+5d0*h1(i-2,1)/6+h1(i,1)/3d0
    h1r(0,i-1,1)=h1(i-2,1)/3d0+5d0*h1(i,1)/6-h1(i+2,1)/6d0
    h1r(0,i-1,2)=11d0*h1(i,1)/6-7d0*h1(i+2,1)/6+h1(i+4,1)/3d0
	
    beta1(0,i,0)=(13d0/12)*(h1(i-4,1)-2*h1(i-2,1)+h1(i,1))**2+&
&	(1d0/4)*(h1(i-4,1)-4*h1(i-2,1)+3*h1(i,1))**2
    beta1(0,i,1)=(13d0/12)*(h1(i-2,1)-2*h1(i,1)+h1(i+2,1))**2+&
&   (1d0/4)*(h1(i-2,1)-h1(i+2,1))**2
    beta1(0,i,2)=(13d0/12)*(h1(i,1)-2*h1(i+2,1)+h1(i+4,1))**2+&
&	(1d0/4)*(3*h1(i,1)-4*h1(i+2,1)+h1(i+4,1))**2

    beta1(1,i,0)=(13d0/12)*(h1(i-4,0)-2*h1(i-2,0)+h1(i,0))**2+&
&	(1d0/4)*(h1(i-4,0)-4*h1(i-2,0)+3*h1(i,0))**2
    beta1(1,i,1)=(13d0/12)*(h1(i-2,0)-2*h1(i,0)+h1(i+2,0))**2+&
&   (1d0/4)*(h1(i-2,0)-h1(i+2,0))**2
    beta1(1,i,2)=(13d0/12)*(h1(i,0)-2*h1(i+2,0)+h1(i+4,0))**2+&
&	(1d0/4)*(3*h1(i,0)-4*h1(i+2,0)+h1(i+4,0))**2
 
    do k=0,1
	  do l=0,2
        alpha1(k,i,l)=b(k,l)/(beta1(k,i,l)+ep)**2	
	  enddo	
	enddo
	
	do k=0,1
	  do l=0,2
	    omega1(k,i,l)=alpha1(k,i,l)/(alpha1(k,i,0)+alpha1(k,i,1)+&
&       alpha1(k,i,2))
      enddo
	enddo	

    do l=0,2
      hn1(i+1,j)=omega1(1,i,l)*h1r(1,i+1,l)+hn1(i+1,j)
	  hp1(i-1,j)=omega1(0,i,l)*h1r(0,i-1,l)+hp1(i-1,j)
    enddo	
  enddo
!---------------------------------------------------------	
	do k=-1,2*m+1,2! 2nd step
        f1(k)=hn1(k,j)+hp1(k,j) 
    enddo 
	
	do i=0,2*m,2  
	  u2(i)=3*u(i,j)/4+(u1(i)-(dt/(2*dx))*(f1(i+1)-f1(i-1)))/4
	enddo

	do k=-6,-2,2 
    u2(k)=u2(2*m+k)
    enddo
    do k=2*m+2,2*m+6,2
    u2(k)=u2(k-2*m)
    enddo
	
	do i=-6,2*m+6,2
     f2(i)=u2(i)**2/2
    enddo

    do i=-6,2*m+6,2
      fp2(i)=(f2(i)+1.5d0*u2(i))/2
      fn2(i)=(f2(i)-1.5d0*u2(i))/2
    enddo

    do i=-6,2*m+6,2
      h2(i,0)=fp2(i)!0 means minus
      h2(i,1)=fn2(i)!1 means plus
    enddo
	
  do i=-2,2*m+2,2
    h1r(1,i+1,0)=h2(i-4,0)/3d0-7d0*h2(i-2,0)/6+11d0*h2(i,0)/6
    h1r(1,i+1,1)=-h2(i-2,0)/6d0+5d0*h2(i,0)/6+h2(i+2,0)/3d0
    h1r(1,i+1,2)=h2(i,0)/3d0+5d0*h2(i+2,0)/6-h2(i+4,0)/6d0
    h1r(0,i-1,0)=-h2(i-4,1)/6d0+5d0*h2(i-2,1)/6+h2(i,1)/3d0
    h1r(0,i-1,1)=h2(i-2,1)/3d0+5d0*h2(i,1)/6-h2(i+2,1)/6d0
    h1r(0,i-1,2)=11d0*h2(i,1)/6-7d0*h2(i+2,1)/6+h2(i+4,1)/3d0
	
    beta1(0,i,0)=(13d0/12)*(h2(i-4,1)-2*h2(i-2,1)+h2(i,1))**2+&
&	(1d0/4)*(h2(i-4,1)-4*h2(i-2,1)+3*h2(i,1))**2
    beta1(0,i,1)=(13d0/12)*(h2(i-2,1)-2*h2(i,1)+h2(i+2,1))**2+&
&   (1d0/4)*(h2(i-2,1)-h2(i+2,1))**2
    beta1(0,i,2)=(13d0/12)*(h2(i,1)-2*h2(i+2,1)+h2(i+4,1))**2+&
&	(1d0/4)*(3*h2(i,1)-4*h2(i+2,1)+h2(i+4,1))**2

    beta1(1,i,0)=(13d0/12)*(h2(i-4,0)-2*h2(i-2,0)+h2(i,0))**2+&
&	(1d0/4)*(h2(i-4,0)-4*h2(i-2,0)+3*h2(i,0))**2
    beta1(1,i,1)=(13d0/12)*(h2(i-2,0)-2*h2(i,0)+h2(i+2,0))**2+&
&   (1d0/4)*(h2(i-2,0)-h2(i+2,0))**2
    beta1(1,i,2)=(13d0/12)*(h2(i,0)-2*h2(i+2,0)+h2(i+4,0))**2+&
&	(1d0/4)*(3*h2(i,0)-4*h2(i+2,0)+h2(i+4,0))**2
 
    do k=0,1
	  do l=0,2
        alpha1(k,i,l)=b(k,l)/(beta1(k,i,l)+ep)**2	
	  enddo	
	enddo
	
	do k=0,1
	  do l=0,2
	    omega1(k,i,l)=alpha1(k,i,l)/(alpha1(k,i,0)+alpha1(k,i,1)+&
&       alpha1(k,i,2))
      enddo
	enddo	

    do l=0,2
      hn2(i+1,j)=omega1(1,i,l)*h1r(1,i+1,l)+hn2(i+1,j)
	  hp2(i-1,j)=omega1(0,i,l)*h1r(0,i-1,l)+hp2(i-1,j)
    enddo	
  enddo
!---------------------------------------------------------------	
	do k=-1,2*m+1,2! 3rd step
       f2(k)=hn2(k,j)+hp2(k,j)
    enddo
	
	do i=0,2*m,2  
	  u(i,j+1)=u(i,j)/3+2*(u2(i)-(dt/(2*dx))*(f2(i+1)-f2(i-1)))/3
	enddo
!-------------------------------------------------------------------
!reset the ghost point to fit the period BC
  do k=-6,-2,2 
    u(k,j+1)=u(2*m+k,j+1)
  enddo
  do k=2*m+2,2*m+6,2
    u(k,j+1)=u(k-2*m,j+1)
  enddo
  
  do i=-6,2*m+6,2
    f(i,j+1)=u(i,j+1)**2/2
  enddo

  do i=-6,2*m+6,2
    fp(i,j+1)=(f(i,j+1)+1.5d0*u(i,j+1))/2
    fn(i,j+1)=(f(i,j+1)-1.5d0*u(i,j+1))/2
  enddo

  do i=-6,2*m+6,2
    h(i,j+1,0)=fp(i,j+1)!0 means minus
    h(i,j+1,1)=fn(i,j+1)!1 means plus
  enddo  
enddo 
!-----------------------------------------------
!output
open(111,file='FDWENOshock_m80.txt')
do i=0,2*m,2
  write(111,*) u(i,n)
enddo
end
