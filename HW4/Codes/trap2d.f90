module trap2d
save
contains

subroutine trap(t_int,fx,n)
    !declare variables
    implicit none
    integer :: i,j
    integer, INTENT(IN) :: n
    real(kind=8) :: a,b,h
    real(kind=8), INTENT(OUT) :: t_int
    real(kind=8), INTENT(IN), dimension(:) :: fx
    real(kind=8), dimension(:), allocatable :: x


i=0
j=0
b=1.d0
a=-1.d0
t_int=0.d0



        
        h=(b-a)/DBLE(n)
        allocate (x(0:n))
        do j=0,n
                x(j)=a+j*h
        enddo
 
        t_int=(fx(1)+fx(n))/2

        do i=1,n-1
                t_int= t_int + fx(i)
        enddo
        
        t_int=h*t_int
        

        deallocate(x)


end subroutine trap
end module trap2d