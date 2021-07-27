module f2_fun
save
contains

function fun(x,y,w,kx,ky) result(f)
  implicit none
  real(kind = 8), INTENT(IN), dimension(:,:) :: x, y
  real(kind = 8) :: x1(1:1,1:size(x)), y1(1:size(y),1:1)
  real(kind = 8), INTENT(IN) :: w, kx, ky
  real(kind=8), dimension(size(x),size(y)) :: f
  integer, dimension(1:2) :: sx,sy
  integer :: i,j,nx,ny

  sx = shape(x)
  sy = shape(y)
  nx = sx(2)
  ny = sy(1)

  x1 = w*cos(kx*x)
  y1 = sin(ky*y)
  
  do i=1,nx
    do j=1,ny
      f(j,i) = x1(1,i)*y1(j,1)
    end do
  end do

end function
end module f2_fun