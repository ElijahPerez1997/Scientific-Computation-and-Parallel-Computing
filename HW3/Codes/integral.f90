program integral
        !define variables
        implicit none 

        !call trap here
        call trap
        !call gauss here
        call gauss

end program integral

!subroutine for f(x)
subroutine fx(k,x,y,n)
        !define variables
        implicit none
        integer :: i,n
        real (kind=8), INTENT(IN) :: k
        !x = inputted nodes, y = function vals
        !nodes = 10 while we test, change later
        real (kind=8), INTENT(INOUT) :: x(0:n), y(0:n)

        !this is the function we are integrating
        do i = 0,n
                y(i) = exp(cos(k*x(i)))
        enddo
end subroutine fx

!sub routine for trapezoidal
subroutine trap

integer :: n, i,j,k
real (kind=8), dimension(:), allocatable :: t_y,x,t_y2
real (kind=8), parameter :: pi=acos(-1.d0), pi2=pi**2, tol=(10d0)**(-10d0)
real (kind=8), parameter :: exact_pi = 2.532131755504017d0,  exact_pi2= 2.452283895d0
real (kind=8) :: h,t_int, t_int2, error, error2,a,b

n=2.d0
i=0
j=0
b=1.d0
a=-1.d0
error=1.d0
t_int=0.d0
!Trapezoidal rule for the function when k=pi
print *, "Trap Integral for k = pi "
do while (error > tol)

        
        h=(b-a)/DBLE(n)
        allocate (x(0:n))
        do j=0,n
                x(j)=a+j*h
        enddo
 
        t_int=(exp(cos((-1.d0)*pi))+exp(cos(pi)))/(2.d0)

        do i=1,n-1
                t_int= t_int + exp(cos(pi*x(i)))
        enddo
        
        t_int=h*t_int
        error = abs(t_int-exact_pi)
        write (*,*) n,',',error
        n=n+1
        deallocate(x)

enddo
!Trapezoidal rule for the function when k=pi^2
n=2.d0
i=0
j=0
b=1.d0
a=-1.d0
error=1.d0
t_int=0.d0
!Trapezoidal rule for the function when k=pi
print *, "Trap Integral for k = pi^2 "
do while (error > tol)

        
        h=(b-a)/DBLE(n)
        allocate (x(0:n))
        do j=0,n
                x(j)=a+j*h
        enddo

        t_int=(exp(cos(-pi2))+exp(cos(pi2)))/(2.d0)

        do i=1,n-1
        t_int= t_int + exp(cos(pi2*x(i)))
        enddo

        t_int=h*t_int
        error = abs(t_int-exact_pi2)
        write (*,*) n,',',error
        n=n+1
        deallocate(x)
        if (n==6000) then
               exit
        endif
enddo



end subroutine trap

!subroutine for Gauss
subroutine gauss
        !define variables
        implicit none
        integer :: n,i,j
        !dynamic arrays
        real (kind=8), dimension(:), allocatable :: g_y, g_node, w, g_y2
        real (kind=8), parameter :: pi=acos(-1.d0), pi2=pi**2, tol=(10d0)**(-10d0)
        real (kind=8) :: g_int, g_int2, error, error2
        real (kind=8), parameter :: exact_pi = 2.532131755504017d0, exact_pi2= 2.452283895d0        
        
        n = 1
        error = 1.d0
print *, "Integrate for K = pi"
        !increase n until error < 10^-10
do while (error > tol) 

        !arrays that depend on n
        !real (kind=8) :: g_y(0:n)
        !real (kind=8) :: g_node(0:n), w(0:n)
        allocate (g_y(0:n))
        allocate (g_node(0:n))
        allocate (w(0:n))

        !g_node is the points we will sub into fx
        call lglnodes(g_node,w,n)
        
        !evaluate fx at each node for pi
        call fx(pi,g_node,g_y,n)

        !calculate the integral as the sum of (y * w)
        g_int = 0.d0
        do i=0,n
                g_int = g_int + (g_y(i)*w(i))
        enddo

        !calculate error
        error = abs(g_int-exact_pi)

        !output n and error for csv file
        write (*,*) n,',',error

        !update n
        n = n+1
        deallocate (g_y,g_node,w)

        !troubleshoot infinite loop
        if (n==2000) then
                exit
        endif

enddo

!repeat process for k=pi^2
print *, "Integrate for K = pi^2"
n = 1
error2 = 1.d0
do while (error2 > tol)
        !arrays that depend on n
        !real (kind=8) :: g_y2(1:n)
        !real (kind=8) :: g_node(0:n), w(0:n)
        allocate (g_y2(0:n))
        allocate (g_node(0:n))
        allocate (w(0:n))

        !g_node is the points we will sub into fx
        call lglnodes(g_node,w,n)

        !evaluate fx at each node for pi^2
        call fx(pi2,g_node,g_y2,n)

        !calculate the integral as the sum of (y * w)
        g_int2 = 0.d0
        do i=0,n
                g_int2 = g_int2 + (g_y2(i)*w(i))
        enddo

        !calculate error
        error2 = abs(g_int2-exact_pi2)

        !output n and error for csv file
        write (*,*) n,',',error2

        !update n
        n = n+1
        
        !troubleshoot infinite loop
        if (n==2000) then
                exit
        endif

        deallocate(g_y2,g_node,w)
enddo

end subroutine gauss


