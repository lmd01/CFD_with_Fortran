program main
implicit none
integer:: i,j,k
double precision::T,X,pi
integer,parameter::m=40,n=10!10,20,40,80,160,320
!Be careful when n is an odd number
double precision::dx,dt
double precision::u(-8*n-2:2*m+8*n+2,0:n),f(-8*n+3:2*m+8*n-3,0:n)
double precision::f1(-8*n+3:2*m+8*n-3),f2(-8*n+3:2*m+8*n-3)
double precision::uap(-8*n+3:2*m+8*n-3,0:n),uan(-8*n+3:2*m+8*n-3,0:n),ua(-8*n:2*m+8*n,0:n)
double precision::u1(-8*n+4:2*m+8*n-4),u2(-8*n+4:2*m+8*n-4)
pi=4.d0*datan(1.d0)
X=2.d0*pi
T=0.5d0!1.5d0,0.5d0
dx=X/(2*m)
dt=T/n
u=0.0
f=0.0
f1=0.0
f2=0.0
ua=0.0
uap=0.0
uan=0.0
u1=0.0
u2=0.0
!write(*,*) X,dx,m
!-----------------------------------------------
!IC&BC
do i=-8*n-2,2*m+8*n+2
  u(i,0)=0.5d0+dsin(i*dx)
enddo 
!-----------------------------------------------
!FV method
do i=-8*n,2*m+8*n,2
  ua(i,0)=(u(i+1,0)+2*u(i,0)+u(i-1,0))/4
enddo
!------------------------------------------------
do i=-8*n+2,2*m+8*n-4,2!2nd order approximation
  uan(i+1,0)=-ua(i-2,0)/6+5*ua(i,0)/6+ua(i+2,0)/3
  uap(i+1,0)=ua(i,0)/3+5*ua(i+2,0)/6-ua(i+4,0)/6
!  write(*,*) uan(i+1,0),uap(i+1,0)
enddo

!------------------------------------------------
!Godunov
do j=0,n-1
  do k=-8*n+3,2*m+8*n-3,2
    if(uan(k,j).lt.uap(k,j)) then
      f(k,j)=min((uan(k,j)**2)/2,(uap(k,j)**2)/2)
	else
	  f(k,j)=max((uan(k,j)**2)/2,(uap(k,j)**2)/2)
    endif	
  enddo
 !---------------------------------------------------------
    do i=-8*n+4,2*m+8*n-4,2 !3rd-order R-K
	  u1(i)=ua(i,j)-(dt/(2*dx))*(f(i+1,j)-f(i-1,j)) 
	  do k=-8*n+5,2*m+8*n-5,2
        if(u1(k-1).lt.u1(k+1)) then
          f1(k)=min((u1(k-1)**2)/2,(u1(k+1)**2)/2)
	    else
	      f1(k)=max((u1(k-1)**2)/2,(u1(k+1)**2)/2)
        endif		
      enddo 
	  u2(i)=3*ua(i,j)/4+(u1(i)-(dt/(2*dx))*(f1(i+1)-f1(i-1)))/4
	  do k=-8*n+7,2*m+8*n-7,2
        if(u2(k-1).lt.u2(k+1)) then
          f2(k)=min((u2(k-1)**2)/2,(u2(k+1)**2)/2)
	    else
	      f2(k)=max((u2(k-1)**2)/2,(u2(k+1)**2)/2)
        endif		
      enddo 
	  ua(i,j+1)=ua(i,j)/3+2*(u2(i)-(dt/(2*dx))*(f2(i+1)-f2(i-1)))/3
!----------------------------------------------------------------------
!update the uan & uap
     uan(i+1,j+1)=-ua(i-2,j+1)/6+5*ua(i,j+1)/6+ua(i+2,j+1)/3
     uap(i+1,j+1)=ua(i,j+1)/3+5*ua(i+2,j+1)/6-ua(i+4,j+1)/6 	  
!write(*,*) ua(i,j+1),f(i+1,j),i+1,f(i-1,j),i-1,j
    enddo	
enddo
!-----------------------------------------------
!error
open(111,file='3rdGodunov_m40.txt')
do i=0,2*m ,2
  write(111,*) ua(i,n)
enddo
end
