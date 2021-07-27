module trap2d
save
contains

subroutine trap(t_int,fx,n,tid)
!$ use omp_lib
    !declare variables
    implicit none
    integer :: i,j,chunk, tid
    integer, INTENT(IN) :: n
    real(kind=8) :: a,b,h
    real(kind=8), INTENT(OUT) :: t_int
    real(kind=8), INTENT(IN), dimension(:) :: fx

  i=0
  j=0
  b=1.d0
  a=-1.d0
  t_int=0.d0
        
  h=(b-a)/DBLE(n)

  !$ call omp_set_num_threads(tid)
  !$omp parallel
         
  t_int=(fx(1)+fx(n))/2

  !$omp do schedule(static) private(i) REDUCTION(+:t_int)
  do i=1,n-1
    t_int = t_int + fx(i)
  enddo
  !$omp end do
        
  t_int=h*t_int 

  !$omp end parallel     

end subroutine trap
end module trap2d