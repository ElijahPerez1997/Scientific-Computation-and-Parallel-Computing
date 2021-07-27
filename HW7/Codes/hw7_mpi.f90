program hw7
  use mpi
  use u_exact
  use f_mod
  use f2_fun
  use laplac
  implicit none

  integer :: n_x, ny, kk, j, k_, nt, i, ierr, nprocs, myid,z
  integer :: status(MPI_STATUS_SIZE)
  integer :: px,ix_off ! the x-proc index and the offset
  integer :: p_left,p_right,px_max
  integer :: nxl       !this is the local size
  integer :: remx, int_sum, request
  character :: str*2

  integer, dimension(1:1) :: NX
  real(kind = 8), dimension(1:5) :: E, H
  real(kind=8) :: mpi_t0, mpi_t1, mpi_dt
  real(kind = 8) :: hx, hy, t_final, dt, w, kx, ky, t, max
  real(kind = 8), dimension(:,:), allocatable :: x, y, u0, f2, L, f, u1, u2,u_ex_final, axy,x1,y1

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  mpi_t0=mpi_wtime() !start timer

  ! Label the processes from 1 to px_max
  px = myid + 1
  px_max = nprocs
 
  w = 10
  kx = 6
  ky = 4
  NX = (/ 400 /)
  kk=1
  n_x = NX(kk)
  ny = n_x !equal number of grid points in x and y directions
  hx = 2.d0/(dble(n_x-1)) !grid length in x direction
  hy = hx !grid length in y direction
  H(kk) = hx

  allocate(y(1:NX(kk),1:1))

  do j = 0,(ny-1)
    y(j+1,1) = -1+j*hx
  end do 

  ! Split up the grid in the x-direction
  nxl = n_x/px_max !number of grid points per processor
  remx = n_x-nxl*px_max !remainder grid points
  if (px <= remx) then
     nxl = nxl + 1 !add any remainders to proc 1
     ix_off = (px-1)*nxl
  else
     ix_off = (remx)*(nxl+1) + (px-(remx+1))*nxl
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_Reduce(nxl,int_sum,1,&
       MPI_INTEGER,MPI_SUM,&
       0,MPI_COMM_WORLD,ierr)
  if(myid == 0) then
     if(n_x .ne. int_sum) then
        write(*,*) 'Something is wrong with the number of points in x-direction: ',&
             n_x,int_sum
     end if
  end if

  ! Determine the neighbours of processor px
  p_left  = px-1 - 1
  p_right = px+1 - 1
  if (px .eq. px_max) p_right = MPI_PROC_NULL
  if (px .eq. 1) p_left = MPI_PROC_NULL  

  ! Allocate memory for various arrays
  !x is 1-by-n_x, y is n_x-by-1
  allocate (x(1:1,0:nxl+1),u0(1:ny,0:nxl+1),f2(1:ny,0:nxl+1))
  allocate (L(1:ny,0:nxl+1),f(1:ny,0:nxl+1),u1(1:ny,0:nxl+1),u2(1:ny,0:nxl+1),u_ex_final(1:ny,0:nxl+1))
  allocate (axy(1:ny,0:nxl+1),x1(1:ny,1:nxl),y1(1:ny,0:nxl+1))

  !initialize x grid
  do j = 1,nxl
    x(1,j) = -1.d0 + dble(j-1+ix_off)*hx
  end do 

  ! send to left recieve from right
  call MPI_Sendrecv(x(1,1),1,MPI_DOUBLE_PRECISION,p_left,123,&
       x(1,nxl+1),1,MPI_DOUBLE_PRECISION,p_right,123,MPI_COMM_WORLD,status,ierr)
  ! send to right recieve from left
  call MPI_Sendrecv(x(1,nxl),1,MPI_DOUBLE_PRECISION,p_right,125,&
       x(1,0),1,MPI_DOUBLE_PRECISION,p_left,125,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  t_final = 2.d0 !Chosen based on assignment reqs

  !Compute the timestep
  dt = 0.5d0*hx !dt < hx/max(a); CFL condition
  nt = FLOOR(t_final/dt) + 1 !same as grid size
  dt = t_final/dble(nt)

  x1=sin(x(1:1,1:nxl))
  y1=cos(y)
  do i=1,nxl
    do j=1,ny
      axy(j,i) = x1(1,i)*y1(j,1)
    end do
  end do
  axy = axy +1 

  u0(:,1:nxl)=exact(0.d0,x(:,1:nxl),y,w,kx,ky) !OK
  f2(:,1:nxl) = fun(x(:,1:nxl),y,w,kx,ky) !OK

  !send/rec axy
  ! Communicate between processors
  ! send to left recieve from right
  call MPI_Sendrecv(axy(:,1),ny,MPI_DOUBLE_PRECISION,p_left,321,&
          axy(:,nxl+1),ny,MPI_DOUBLE_PRECISION,p_right,321,MPI_COMM_WORLD,status,ierr)
  ! send to right recieve from left
  call MPI_Sendrecv(axy(:,nxl),ny,MPI_DOUBLE_PRECISION,p_right,125,&
          axy(:,0),ny,MPI_DOUBLE_PRECISION,p_left,125,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !send/rec u0
  ! send to left recieve from right
  call MPI_Sendrecv(u0(:,1),ny,MPI_DOUBLE_PRECISION,p_left,321,&
          u0(:,nxl+1),ny,MPI_DOUBLE_PRECISION,p_right,321,MPI_COMM_WORLD,status,ierr)
  ! send to right recieve from left
  call MPI_Sendrecv(u0(:,nxl),ny,MPI_DOUBLE_PRECISION,p_right,125,&
          u0(:,0),ny,MPI_DOUBLE_PRECISION,p_left,125,MPI_COMM_WORLD,status,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(nxl==n_x) then
    call lap(u0(:,1:nxl),hx,hy,nxl,ny,L(:,1:nxl),axy(:,1:nxl),ix_off,n_x)
  else if(ix_off == 0) then !first, send 1 extra col at end
    call lap(u0(:,1:nxl+1),hx,hy,nxl,ny,L(:,1:nxl),axy(:,1:nxl+1),ix_off,n_x)
  else if(ix_off == n_x-nxl) then !last, send 1 extra col at start
    call lap(u0(:,0:nxl),hx,hy,nxl,ny,L(:,1:nxl),axy(:,0:nxl),ix_off,n_x)
  else !middle, send 2 extra col
    call lap(u0(:,0:nxl+1),hx,hy,nxl,ny,L(:,1:nxl),axy(:,0:nxl+1),ix_off,n_x)
  end if

  call forcing(0.d0,x(:,1:nxl),y,w,kx,ky,f(:,1:nxl)) !OK
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i=1,ny
    do j=1,nxl
      u1(i,j) = u0(i,j) + (dt*f2(i,j)) + 0.5d0*(dt**2)*(L(i,j)+f(i,j))
    end do
  end do

  u1(1:1,1:nxl) = exact31(dt,x(:,1:nxl),-1.d0,w,kx,ky) !sets first row
  u1(ny:ny,1:nxl) = exact31(dt,x(:,1:nxl),1.d0,w,kx,ky) !sets last row
  if(myid==0) then
    u1(:,1:1) = exact21(dt,-1.d0,y,w,kx,ky) !sets first col !OK
  end if
  if(ix_off == n_x-nxl) then
    u1(:,nxl:nxl) = exact21(dt,1.d0,y,w,kx,ky) !sets last col
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr) 

  do k_=0,nt-2
    t = (k_+1)*dt
    !obtain right hand side
    call forcing(t,x(:,1:nxl),y,w,kx,ky,f(:,1:nxl))

    !send/rec axy, send to left recieve from right
    !non-blocking send
    call MPI_Isend(axy(:,1),ny,MPI_DOUBLE_PRECISION,p_left,321,MPI_COMM_WORLD,request,ierr)
    call MPI_RECV(axy(:,nxl+1),ny,MPI_DOUBLE_PRECISION,p_right,321,MPI_COMM_WORLD,status,ierr)
    ! send to right recieve from left
    call MPI_Isend(axy(:,nxl),1,MPI_DOUBLE_PRECISION,p_right,125,MPI_COMM_WORLD,request,ierr)
    call MPI_RECV(axy(:,0),1,MPI_DOUBLE_PRECISION,p_left,125,MPI_COMM_WORLD,status,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !send/rec u1, send to left recieve from right
    call MPI_Isend(u1(:,1),ny,MPI_DOUBLE_PRECISION,p_left,99,MPI_COMM_WORLD,request,ierr)
    call MPI_RECV(u1(:,nxl+1),ny,MPI_DOUBLE_PRECISION,p_right,99,MPI_COMM_WORLD,status,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    ! send to right recieve from left
    call MPI_Isend(u1(:,nxl),ny,MPI_DOUBLE_PRECISION,p_right,88,MPI_COMM_WORLD,request,ierr)
    call MPI_RECV(u1(:,0),ny,MPI_DOUBLE_PRECISION,p_left,88,MPI_COMM_WORLD,status,ierr)  
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(nxl==n_x) then
      call lap(u1(:,1:nxl),hx,hy,nxl,ny,L(:,1:nxl),axy(:,1:nxl),ix_off,n_x)
    else if(ix_off == 0) then !first, send 1 extra col at end
      call lap(u1(:,1:nxl+1),hx,hy,nxl,ny,L(:,1:nxl),axy(:,1:nxl+1),ix_off,n_x)
    else if(ix_off == n_x-nxl) then !last, send 1 extra col at start
      call lap(u1(:,0:nxl),hx,hy,nxl,ny,L(:,1:nxl),axy(:,0:nxl),ix_off,n_x)
    else !middle, send 2 extra col
      call lap(u1(:,0:nxl+1),hx,hy,nxl,ny,L(:,1:nxl),axy(:,0:nxl+1),ix_off,n_x)
    end if

    !March in time
    u2(:,1:nxl) = 2*u1(:,1:nxl)-u0(:,1:nxl)+(dt**2)*(L(:,1:nxl)+f(:,1:nxl))
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !update BC
    u2(1:1,1:nxl) = exact31(t+dt,x(:,1:nxl),-1.d0,w,kx,ky) 
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)    
    u2(ny:ny,1:nxl) = exact31(t+dt,x(:,1:nxl),1.d0,w,kx,ky) 
    if(myid==0) then
      u2(:,1:1) = exact21(t+dt,-1.d0,y,w,kx,ky)
    end if
    if(ix_off == n_x-nxl) then
      u2(:,nxl:nxl) = exact21(t+dt,1.d0,y,w,kx,ky)
    end if

    !switch solution at different time levels
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     do i=1,ny
       do j=1,nxl
         u0(i,j) = u1(i,j)
       end do
     end do
     do i=1,ny
       do j=1,nxl
         u1(i,j) = u2(i,j)
       end do
     end do

  end do

  !Calculate exact solution
  u_ex_final(:,1:nxl)= exact(t_final,x(:,1:nxl),y,w,kx,ky)

  ! Compute the error
  E(kk)=maxval(maxval(abs(u2(:,1:nxl)-u_ex_final(:,1:nxl)),dim=2))

  call MPI_Barrier(MPI_COMM_WORLD,ierr) 
  mpi_t1=mpi_wtime() !finish time
  mpi_dt=mpi_t1-mpi_t0 !total time elapsed

  if(myid == 0) write(*,*) "With", nprocs, "  nodes the wall clock time = ", mpi_dt, " sec."
 
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !print error
  write(*,*)'Processor ', myid, 'Error = ', E(kk)

  call MPI_Barrier(MPI_COMM_WORLD,ierr) 
  call mpi_finalize(ierr) !end MPI section

end program hw7


