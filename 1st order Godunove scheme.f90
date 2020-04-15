program main
implicit none
integer:: i,j,k
double precision::T,X,pi
integer,parameter::m=80,n=20!10,20,40,80,160,320
double precision::dx,dt
double precision::u(-n:2*m+n,0:n),f(-n+1:2*m+n-1,0:n),ua(-n:2*m+n,0:n)
pi=4.d0*datan(1.d0)
X=2.d0*pi
T=1.5d0!1.5d0,0.5d0
dx=X/(2*m)
dt=T/n
!-----------------------------------------------
!IC&BC
do i=-n,2*m+n
  u(i,0)=0.5d0+dsin(i*dx)
enddo 
do i=-n+2,2*m+n-2,2
  ua(i,0)=(u(i+2,0)+2*u(i,0)+u(i-2,0))/4
enddo

!------------------------------------------------
!Godunov
do j=0,n-1
  do k=-n+1,2*m+n-1,2
    if(ua(k-1,j).lt.ua(k+1,j)) then
      f(k,j)=min((ua(k-1,j)**2)/2,(ua(k+1,j)**2)/2)
	else
	  f(k,j)=max((ua(k-1,j)**2)/2,(ua(k+1,j)**2)/2)
    endif
  enddo
    do i=-n+2,2*m+n-2,2
	  ua(i,j+1)=ua(i,j)-(dt/(2*dx))*(f(i+1,j)-f(i-1,j)) 
    enddo
enddo
!-----------------------------------------------
!error
open(111,file='Godunovshock_m80.txt')
do i=0,2*m,2
  write(111,*) ua(i,n)
enddo
end
