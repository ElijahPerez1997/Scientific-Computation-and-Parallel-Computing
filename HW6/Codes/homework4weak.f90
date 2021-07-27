program hwk4
!weak scale
!$ use omp_lib
  use xycoord   !use the module xycoord to set the mapping 
  use trap2d
  implicit none !NOTE: use module-name must preceed the implicit none
  integer :: nr,ns,i,j,k,chunk,tid,nthreads,count
  real(kind = 8) :: hr,hs,det,Eux,Euy,fout,int1,int2,int3,int4,sum
  real(kind = 8) :: error, start(1:10), finish(1:10), t_av
  real(kind = 8), dimension(:), allocatable :: r,s
  real(kind = 8), dimension(:,:), allocatable :: u,ur,us,xr,xs,yr,ys,rx,ry,sx,sy,uxExact,uyExact
  real(kind = 8), dimension(:,:), allocatable :: fx,xc,yc,ux,uy,jac,uxx,uyy,uyr,uys,uxr,uxs,laplace

  !OUTER LOOP
  do tid = 1,16
  write(*,*) 'Num of threads = ', tid

 k = nint((tid**0.5d0)*200.d0)
  nr = k
  ns = k
  write(*,*) 'Grid size = ', k

  ! Allocate memory for various arrays
  allocate (r(0:nr),s(0:ns),u(0:nr,0:ns),ur(0:nr,0:ns),us(0:nr,0:ns),xr(0:nr,0:ns),uxExact(0:nr,0:ns),uyExact(0:nr,0:ns))
  allocate (xs(0:nr,0:ns),yr(0:nr,0:ns),ys(0:nr,0:ns),rx(0:nr,0:ns),ry(0:nr,0:ns),fx(0:nr,0:ns))
  allocate (sx(0:nr,0:ns),sy(0:nr,0:ns),ux(0:nr,0:ns),uy(0:nr,0:ns),uxr(0:nr,0:ns),uxs(0:nr,0:ns),uyr(0:nr,0:ns),uys(0:nr,0:ns))
  allocate (xc(0:nr,0:ns),yc(0:nr,0:ns),jac(0:nr,0:ns),uxx(0:nr,0:ns),uyy(0:nr,0:ns),laplace(0:nr,0:ns))
  
  count = 1  

  do while (count < 11)

  !start timer
  call cpu_time(start(count))
  !$ start(count) = omp_get_wtime()

  hr = 2.d0/dble(nr)
  hs = 2.d0/dble(ns)
  !h = [2/20, 2/800]
  nthreads = 0
  !$ call omp_set_num_threads(tid)
  !$ nthreads = omp_get_num_threads()

  !$omp parallel

  !$omp do schedule(static) private(i)
  do i = 0,nr
     r(i) = -1.d0 + dble(i)*hr
  end do
  !$omp end do nowait

  !$omp do schedule(static) private(i)
  do i = 0,ns
     s(i) = -1.d0 + dble(i)*hs
  end do
  !$omp end do nowait

  !$omp do collapse(2) schedule(static) private(i)
  do j = 0,ns
     do i = 0,nr
        xc(i,j) = x_coord(r(i),s(j))
        yc(i,j) = y_coord(r(i),s(j))
     end do
  end do
  !$omp end do

  !$omp do collapse(2) schedule(static) private(i)
  do j = 0,ns
     do i = 0,nr
        u(i,j) = sin(xc(i,j))*cos(yc(i,j))
	uxExact(i,j)=cos(xc(i,j))*cos(yc(i,j))
	uyExact(i,j)=-sin(xc(i,j))*sin(yc(i,j))
     end do
  end do
  !$omp end do

  ! Differentiate in the r-direction U_R
  !$omp do schedule(static) private(i)
  do i = 0,ns
     call differentiate(u(0:nr,i),ur(0:nr,i),hr,nr,tid)
  end do
  !$omp end do nowait

  ! Differentiate in the s-direction U_S
  !$omp do schedule(static) private(i)
  do i = 0,nr
     call differentiate(u(i,0:ns),us(i,0:ns),hs,ns,tid)
  end do
  !$omp end do nowait

  ! Differentiate in the r-direction
  !$omp do schedule(static) private(i)
  do i = 0,ns
     call differentiate(xc(0:nr,i),xr(0:nr,i),hr,nr,tid)
  end do
  !$omp end do nowait

  ! Differentiate in the s-direction
  !$omp do schedule(static) private(i)
  do i = 0,nr
     call differentiate(xc(i,0:ns),xs(i,0:ns),hs,ns,tid)
  end do
  !$omp end do nowait
  
  ! Differentiate in the s-direction
  !$omp do schedule(static) private(i)
  do i = 0,nr
     call differentiate(yc(0:nr,i),yr(0:nr,i),hr,nr,tid)
  end do
  !$omp end do nowait
  
  ! Differentiate in the s-direction
  !$omp do schedule(static) private(i)
  do i = 0,ns
     call differentiate(yc(i,0:ns),ys(i,0:ns),hs,ns,tid)
  end do
  !$omp end do

  !Calculate rx, ry, sx, sy as the inverse

  !$omp single
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
  !$omp end single

  !calculate error 
  !$omp single
  do i=0,nr
	do j=0,nr
		fx(i,j)=(ux(i,j)+uy(i,j)-uxExact(i,j)-uyExact(i,j))**2
	end do
  end do
  !$omp end single

  int2=0
  int4=0

  !$omp do schedule(static) private(i,int1) REDUCTION(+:int2)
  do i=0,nr
	call trap(int1,fx(:,i),nr,tid)
	int2=int1+int2
  enddo
  !$omp end do

!$omp single
  int2=int2/dble(nr)

  do i=0,nr	
	call trap(int3,fx(i,:),nr,tid)
	int4=int3+int4
  enddo

  int4=int4/dble(nr)
  int4=int4+int2

  error = int4**0.5d0  

  !$omp end single
  
  !$omp end parallel
  call cpu_time(finish(count))
  !$ finish(count) = omp_get_wtime()
  !find time average
  sum = 0.d0
  if (count == 4) then
    do i=1,4
      sum = sum + (finish(i) - start(i))
    end do
    t_av = sum/dble(count)
  end if

  count = count + 1
  end do !while

  !record time
  write (*,*) 'Elapsed time = ', t_av  

  !record error
  write (*,*) 'Error = ', error  
  
  deallocate(r,s,u,ur,us,xr,xs,yr,ys,rx,ry,sx,sy,xc,yc,ux,uy,jac,uyr,uys,uyy,uxr,uxs,uxx,laplace,uxExact,uyExact,fx)

  end do !threads

end program hwk4





