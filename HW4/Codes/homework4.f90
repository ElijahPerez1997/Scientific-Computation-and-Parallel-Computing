program hwk4
  use xycoord   !use the module xycoord to set the mapping 
  use trap2d
  implicit none !NOTE: use module-name must preceed the implicit none
  integer :: nr,ns,i,j,k,n(1:10)
  real(kind = 8) :: hr,hs,det,Eux,Euy,fout,int1,int2,int3,int4,error(1:10),error2(1:10), h(1:10)
  real(kind = 8), dimension(:), allocatable :: r,s
  real(kind = 8), dimension(:,:), allocatable :: u,ur,us,xr,xs,yr,ys,rx,ry,sx,sy,uxExact,uyExact,fx
  real(kind = 8), dimension(:,:), allocatable :: xc,yc,ux,uy,jac,uxx,uyy,uyr,uys,uxr,uxs,laplace

  j = 5
  do i = 1,10
    n(i) = j
    j = j*2
  end do

  do k = 1,6
	

  nr = n(k)
  ns = n(k)
  
  ! Allocate memory for various arrays
  allocate (r(0:nr),s(0:ns),u(0:nr,0:ns),ur(0:nr,0:ns),us(0:nr,0:ns),xr(0:nr,0:ns),uxExact(0:nr,0:ns),uyExact(0:nr,0:ns))
  allocate (xs(0:nr,0:ns),yr(0:nr,0:ns),ys(0:nr,0:ns),rx(0:nr,0:ns),ry(0:nr,0:ns),fx(0:nr,0:ns))
  allocate (sx(0:nr,0:ns),sy(0:nr,0:ns),ux(0:nr,0:ns),uy(0:nr,0:ns),uxr(0:nr,0:ns),uxs(0:nr,0:ns),uyr(0:nr,0:ns),uys(0:nr,0:ns))
  allocate (xc(0:nr,0:ns),yc(0:nr,0:ns),jac(0:nr,0:ns),uxx(0:nr,0:ns),uyy(0:nr,0:ns),laplace(0:nr,0:ns))

  hr = 2.d0/dble(nr)
  hs = 2.d0/dble(ns)
  h(k) = hr

  do i = 0,nr
     r(i) = -1.d0 + dble(i)*hr
  end do

  do i = 0,ns
     s(i) = -1.d0 + dble(i)*hs
  end do

  do j = 0,ns
     do i = 0,nr
        xc(i,j) = x_coord(r(i),s(j))
        yc(i,j) = y_coord(r(i),s(j))
     end do
  end do
  
  call  printdble2d(xc,nr,ns,'x.txt')
  call  printdble2d(yc,nr,ns,'y.txt')
  
  do j = 0,ns
     do i = 0,nr
        u(i,j) = sin(xc(i,j))*cos(yc(i,j))
	uxExact(i,j)=cos(xc(i,j))*cos(yc(i,j))
	uyExact(i,j)=-sin(xc(i,j))*sin(yc(i,j))
     end do
  end do
  
  ! Differentiate in the r-direction U_R
  do i = 0,ns
     call differentiate(u(0:nr,i),ur(0:nr,i),hr,nr)
  end do
  
  ! Differentiate in the s-direction U_S
  do i = 0,nr
     call differentiate(u(i,0:ns),us(i,0:ns),hs,ns)
  end do

  ! Differentiate in the r-direction
  do i = 0,ns
     call differentiate(xc(0:nr,i),xr(0:nr,i),hr,nr)
  end do

  ! Differentiate in the s-direction
  do i = 0,nr
     call differentiate(xc(i,0:ns),xs(i,0:ns),hs,ns)
  end do
  
  ! Differentiate in the s-direction
  do i = 0,nr
     call differentiate(yc(0:nr,i),yr(0:nr,i),hr,nr)
  end do
  
  ! Differentiate in the s-direction
  do i = 0,ns
     call differentiate(yc(i,0:ns),ys(i,0:ns),hs,ns)
  end do

  !Calculate rx, ry, sx, sy as the inverse
  do i = 0,ns
        do j = 0,nr
                det = xr(i,j)*ys(i,j) - xs(i,j)*yr(i,j)
                jac(i,j)=det
                rx(i,j) = ys(i,j)/det
                sx(i,j) = -yr(i,j)/det
                ry(i,j) = -xs(i,j)/det
                sy(i,j) = xr(i,j)/det
                !also calculate ux, uy
                ux(i,j) = ur(i,j)*rx(i,j) + us(i,j)*sx(i,j)
                uy(i,j) = ur(i,j)*ry(i,j) + us(i,j)*sy(i,j)
        end do
  end do
  
  
  !calculate second deriviatives

  

  !calculate error 
  do i=0,nr
	do j=0,nr
		fx(i,j)=(ux(i,j)+uy(i,j)-uxExact(i,j)-uyExact(i,j))**2
	end do
  end do
  int2=0
  int4=0
  do i=0,nr
	call trap(int1,r,nr)
	int2=int1+int2
  enddo
  int2=int2/dble(nr)
  do i=0,nr
	
	call trap(int3,fx(i,:),nr)
	int4=int3+int4
  enddo
  int4=int4/dble(nr)
  int4=int4+int2
  	write(*,*) int4,nr
	
  
  error(k) = int4**0.5d0  
  
  
  deallocate(r,s,u,ur,us,xr,xs,yr,ys,rx,ry,sx,sy,xc,yc,ux,uy,jac,uyr,uys,uyy,uxr,uxs,uxx,laplace,uxExact,uyExact,fx)

end do

open(1,file=trim('error.txt'),status='unknown')
do i = 1,10
  write(1,fmt='(E24.16)',advance='no') error(i)
end do

open(1,file=trim('h.txt'),status='unknown')
do i = 1,10
  write(1,fmt='(E24.16)',advance='no') h(i)
end do

!open(1,file=trim('errorLaplace.txt'),status='unknown')
!do i = 1,10
!  write(1,fmt='(E24.16)',advance='no') error2(i)
!end do

end program hwk4


