module xycoord
  real(kind = 8), parameter :: pi = acos(-1.d0)
  save
  contains
  
real(kind=8) function x_coord(r,s)
    implicit none
    real(kind=8) r,s
  !function 1
    !x_coord = ((2.d0*s)**2+r**2)*sin(s+r)
  !function 2
    !x_coord = (r**2+2.d0*s)
  !function 3
    x_coord = r+0.1*s
end function x_coord

real(kind=8) function y_coord(r,s)
    implicit none
    real(kind=8) r,s
  !function 1
    !y_coord = (3.d0*s+r**2)*cos(r+s)
  !function 2
    !y_coord = (3.d0*r+s**2)
  !function 3
    y_coord = s
end function y_coord
    
end module xycoord
