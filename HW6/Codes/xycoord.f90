module xycoord
  real(kind = 8), parameter :: pi = acos(-1.d0)
  save
  contains
  
real(kind=8) function x_coord(r,s)
    implicit none
    real(kind=8) r,s
    x_coord = r+0.1*s
end function x_coord

real(kind=8) function y_coord(r,s)
    implicit none
    real(kind=8) r,s
    y_coord = s
end function y_coord
    
end module xycoord
