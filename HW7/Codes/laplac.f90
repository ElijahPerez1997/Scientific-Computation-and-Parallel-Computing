module laplac
save
contains

subroutine lap(u,hx,hy,nxl,ny,L,a,ix_off,nx)
  implicit none
  real(kind = 8), dimension(:,:), INTENT(inOUT) :: L
  real(kind = 8), INTENT(IN), dimension(:,:) :: u, a
  real(kind=8), INTENT(IN) :: hx, hy
  real(kind=8) :: t1,t2,a1,a2,a3,a4,a5
  integer, INTENT(IN) :: nxl, nx, ny, ix_off
  integer :: j, i, z

  if(nxl==nx) then
    z = 1
  else if(ix_off == nx-nxl) then
    z = 0
  else
    z = -1
  end if

  do j=2,nxl-z !already localized !get rid of minus 1 to calculate extra col
    do i=2,ny-1 
      t1 = 2.d0*hx**2.d0 
      t2 = 2.d0*hy**2.d0 
      a1 = ((a(i+1,j)+a(i,j))/t1)
      a2 = ((a(i,j)+a(i-1,j))/t1)
      a3 = ((a(i,j+1)+a(i,j))/t2)
      a4 = ((a(i,j)+a(i,j-1))/t2)
      a5 = ((a(i+1,j)+2*a(i,j)+a(i-1,j))/t1)+((a(i,j+1)+2*a(i,j)+a(i,j-1))/t2)
      L(i,j)=a1*u(i+1,j)+a2*u(i-1,j)+a3*u(i,j+1)+a4*u(i,j-1)-(a5*u(i,j))    
    end do
  end do

  !return only the columns needed
  if(nxl==nx) then
    L = L
  else if(ix_off == 0) then !first proc, get rid of last col
    L = L(:,1:nxl) !sent (:,1:nxl+1)
  else if(ix_off == nx-nxl) then !last proc, get rid of 1st col
    L = L(:, 2:nxl+1) !sent (:,1:nxl+1)
  else !middle, get rid of two outer cols
    L = L(:,2:nxl+1) !sent (:,0:nxl+1)
  end if

end subroutine lap

end module laplac