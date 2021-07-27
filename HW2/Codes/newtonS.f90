!
! A simple example: solving an equation f(x) = 0
! using Newton's method
!
! This is the template file used for the scripted version  
!
program newton
  
  implicit none
  double precision :: f,fp,x,dx, curr, e_n, e_np1, convergence, convergence_2, tol
  integer :: function_count, count

  ! Here we try to find the solution to f(x) = 0
  
  write(*,*) ' sin(x)+cos(x*x) '
  x = -0.5d0
  e_np1 = 0
  curr = 0
  count = 0
  tol = 10**(-15)
  
  do while (abs(x-curr) > tol)
    count = count + 1
    e_n = abs(x-curr)
    curr = x
    !curr is the current guess, x is the next guess ie. x_(n+1)
    f = ffun(x)
    fp = fpfun(x)
    dx = -f/fp
    x = x + dx
    IF (count > 1) THEN
    convergence = (abs(x-curr))/(e_n)
    e_np1 = abs(x-curr)
    convergence_2 = (abs(x-curr))/(e_n)**2
    END IF
    write(*,*) count,',',e_np1,',',convergence,',',convergence_2
  end do

write(*,*) 'function'

contains

  double precision function ffun(x)
    implicit none
    double precision :: x

    ffun = sin(x)+cos(x*x)

  end function ffun

  double precision function fpfun(x)
    implicit none
    double precision :: x

    fpfun = cos(x)-2.d0*x*sin(x*x)

  end function fpfun  

end program newton
