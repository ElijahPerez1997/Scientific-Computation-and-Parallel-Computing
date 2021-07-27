program hw7
  use u_exact
  use f_mod
  use f2_fun
  use laplac
  implicit none
  integer :: n_x, ny, kk, j, k_, nt, i
  integer, dimension(1:5) :: NX
  real(kind = 8), dimension(1:5) :: E, H
  real(kind = 8) :: hx, hy, h_, t_final, dt, w, kx, ky, t, max
  real(kind = 8), dimension(:,:), allocatable :: x, y, u0, f2, L, f, u1, u2,u_ex_final, axy,x1,y1
!ALL INDEXES START AT 1

w = 10
kx = 6
ky = 4

NX = (/ 21, 41, 81, 161, 321 /)

!OUTER LOOP
do kk = 1,5
!kk=1
!write(*,*)'kk = ', kk
  ! Allocate memory for various arrays
  !x is 1-by-n_x, y is n_x-by-1
  allocate (x(1:1,1:NX(kk)),y(1:NX(kk),1:1),u0(1:NX(kk),1:NX(kk)),f2(1:NX(kk),1:NX(kk)))
  allocate (L(1:NX(kk),1:NX(kk)),f(1:NX(kk),1:NX(kk)),u1(1:NX(kk),1:NX(kk)),u2(1:NX(kk),1:NX(kk)),u_ex_final(1:NX(kk),1:NX(kk)))
  allocate (axy(1:NX(kk),1:NX(kk)),x1(1:NX(kk),1:NX(kk)),y1(1:NX(kk),1:NX(kk)))

  n_x = NX(kk)
  ny = n_x !equal number of grid points in x and y directions
  hx = 2/(dble(n_x)-1) !grid length in x direction
  hy = hx !grid length in y direction
  H(kk) = hx**2
  
  h_=(2)/dble(n_x-1) !step size, from -1 to 1
  do j = 0,(n_x-1)
    x(1,j+1) = -1+j*h_
  end do 

  do j = 0,(n_x-1)
    y(j+1,1) = -1+j*h_
  end do 



  !axy = MATMUL(transpose(sin(x)),transpose(cos(y))) +1
  !axy = transpose(axy) !OK
  x1=sin(x)
  y1=cos(y)
  do i=1,n_x
    do j=1,ny
      axy(j,i) = x1(1,i)*y1(j,1)
    end do
  end do
  axy = axy +1 !OK

  t_final = 2.d0
  dt = 0.5d0*hx
  nt = FLOOR(t_final/dt) + 1
  dt = t_final/dble(nt)

  u_ex_final = exact(t_final,x,y,w,kx,ky) !OK

  u0=exact(0.d0,x,y,w,kx,ky) !OK
  f2 = fun(x,y,w,kx,ky) !OK
  call lap(u0,hx,hy,n_x,ny,L, axy) !OK
  call forcing(0.d0,x,y,w,kx,ky,f) !OK
!write(*,*) 'f'
!do i=1,5
!  write(*,*) f(i,:)
!end do

  u1 = u0 + (dt*f2) + 0.5d0*(dt**2)*(L+f) !OK

  u1(1:1,:) = exact31(dt,x,-1.d0,w,kx,ky) !sets first row

  u1(ny:ny,:) = exact31(dt,x,1.d0,w,kx,ky) !sets last row

  u1(:,1:1) = exact21(dt,-1.d0,y,w,kx,ky) !sets first col !OK

  u1(:,n_x:n_x) = exact21(dt,1.d0,y,w,kx,ky) !sets last col

  do k_=0,nt-2
    t = (k_+1)*dt
    !obtain right hand side
    call forcing(t,x,y,w,kx,ky,f)
    call lap(u1,hx,hy,n_x,ny,L,axy)
    !March in time
    u2 = 2*u1-u0+(dt**2)*(L+f)
    !update BC
    u2(1:1,:) = exact31(t+dt,x,-1.d0,w,kx,ky)
    u2(ny:ny,:) = exact31(t+dt,x,1.d0,w,kx,ky)
    u2(:,1:1) = exact21(t+dt,-1.d0,y,w,kx,ky)
    u2(:,n_x:n_x) = exact21(t+dt,1.d0,y,w,kx,ky)
    !switch solution at different time levels
    u0 = u1
    u1 = u2
  end do

  E(kk)=maxval(maxval(abs(u2-u_ex_final),dim=2))
  
  !deallocate
  deallocate(x,y,x1,y1,u0,axy,f2,L,f,u1,u2,u_ex_final)

  !end outer loop
end do

!print error
  write(*,*) 'error'
do k_=1,5
  write(*,*) E(k_)
end do

end program hw7


