module u_exact
save
contains

function exact(t,x,y,w,kx,ky) result(u) !neither arg is constant
  implicit none
  real(kind = 8), INTENT(IN), dimension(:,:) :: x, y
  real(kind = 8), INTENT(IN) :: t, w, kx, ky
  real(kind=8), dimension(size(x),size(y)) :: u
  real(kind = 8) :: x1(1:1,1:size(x)), y1(1:size(y),1:1)
  integer, dimension(1:2) :: sx,sy
  integer :: i,j,nx,ny

  sx = shape(x)
  sy = shape(y)
  nx = sx(2)
  ny = sy(1)

  x1 = sin(w*t - kx*x)
  y1 = sin(ky*y)

  do i=1,nx
    do j=1,ny
      u(j,i) = x1(1,i)*y1(j,1)
    end do
  end do

end function !checked this function. good.

function exact21(t,x,y,w,kx,ky) result(u1) !2nd arg is constant and u is 1d
  implicit none
  real(kind = 8), INTENT(IN), dimension(:,:) :: y
  real(kind = 8), INTENT(IN) :: t, w, kx, ky, x
  real(kind=8), dimension(size(y),1) :: u1

  u1 = (sin(w*t - kx*x))*((sin(ky*y)))

end function

function exact31(t,x,y,w,kx,ky) result(u2) !3rd arg is constant and u is 1d
  implicit none
  real(kind = 8), INTENT(IN), dimension(:,:) :: x
  real(kind = 8), INTENT(IN) :: t, w, kx, ky, y
  real(kind=8), dimension(1,size(x)) :: u2

 u2 = (sin(w*t - kx*x)) * (sin(ky*y))
end function !pretty sure this is ok

end module u_exact