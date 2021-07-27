module f_mod
save
contains

subroutine forcing(t,x,y,w,kx,ky,force) 
  implicit none
  real(kind = 8), dimension(:,:), INTENT(INOUT) :: force
  real(kind = 8), INTENT(IN), dimension(:,:) :: x, y
  real(kind = 8) :: x1(1:1,1:size(x)), y1(1:size(y),1:1)
(1:21,1:1),a4(1:1,1:3)
  real(kind = 8), dimension(:,:),allocatable :: auxx,c,auyy,auxx2,auyy2,utt,t1,t2,t3,t5,a1,a2,a3,a4
  real(kind = 8), INTENT(IN) :: t, w, kx, ky
  integer :: i,j,nx,ny,sx(1:2),sy(1:2),k

  sx = shape(x) !x is already localized
  nx = size(x) !ghost points ? sent without
  sy = shape(y)
  ny = size(y)

  allocate(auxx(1:ny,1:nx),c(1:ny,1:nx),auxx2(1:ny,1:nx),auyy(1:ny,1:nx),auyy2(1:ny,1:nx),t3(1:ny,1:nx),utt(1:ny,1:nx))
  allocate(t1(1:ny,1:1),t2(1:1,1:nx),t5(1:ny,1:1),a3(1:ny,1:1),a1(1:1,1:nx),a2(1:1,1:nx),a4(1:1,1:nx))

  t1 = sin(ky*y)
  t2 = w*t-kx*x
  x1 = sin(x)
  y1 = cos(y)

  do i=1,nx
    do j=1,ny
      t3(j,i) = x1(1,i)*y1(j,1)
    end do
  end do
  t3 = t3+1

  a1 = (-kx*cos(w*t-kx*x))
  a2 = cos(x)
  a3 = cos(y) 
  a4 = (kx**2)*sin(t2)

  do i=1,nx
      a1(1,i)= a1(1,i)*a2(1,i)
  end do

  do i=1,nx
    do j=1,ny
      auxx(j,i)= a1(1,i)*a3(j,1)
    end do
  end do

  do i=1,ny
    do j=1,nx
      auxx(i,j)= t1(i,1)*auxx(i,j)
    end do
  end do

  do i=1,ny
    do j=1,nx
      auxx2(i,j)= a4(1,j)*t1(i,1)
    end do
  end do

 do i=1,ny
    do j=1,nx
      auxx2(i,j)= auxx2(i,j)*t3(i,j)
    end do
  end do

  auxx=auxx-auxx2 !OK

  t5 = cos(ky*y)
  a2 = -sin(x)
  a3 = sin(y)
  a4 = ky*sin(t2)

  do i=1,ny
    do j=1,nx
      auyy(i,j)= a2(1,j)*a3(i,1)
    end do
  end do

  do i=1,ny
    do j=1,nx
      auyy(i,j)= t5(i,1)*auyy(i,j)
    end do
  end do

  a3 = ky*t1

  do i=1,ny
    do j=1,nx
      auyy2(i,j)= a3(i,1)*t3(i,j)
    end do
  end do

  auyy2 = auyy-auyy2

  do i=1,ny
    do j=1,nx
      auyy(i,j)= a4(1,j)*auyy2(i,j)
    end do
  end do

  a2 = -w**2*sin(t2)

 do i=1,ny
    do j=1,nx
      utt(i,j)= t1(i,1)*a2(1,j)
    end do
  end do

  force=utt-auxx-auyy

end subroutine forcing
end module f_mod