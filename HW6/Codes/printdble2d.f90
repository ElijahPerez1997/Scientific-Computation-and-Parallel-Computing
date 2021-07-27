subroutine printdble2d(u,nx,ny,str,tid)
!$ use omp_lib

  implicit none
  integer, intent(in) :: nx,ny,tid
  real(kind = 8), intent(in) :: u(0:nx,0:ny)
  character(len=*), intent(in) :: str
  integer :: i,j,chunk

  open(2,file=trim(str),status='unknown')

  !$ call omp_set_num_threads(tid)
  !$omp parallel
  chunk = 1
  !$omp do schedule(dynamic,chunk) private(i)
  do j=0,ny,1
     do i=0,nx,1
        write(2,fmt='(E24.16)',advance='no') u(i,j)
     end do
     write(2,'()')
  end do
  !$omp end do
  !$omp end parallel
  close(2)
end subroutine printdble2d
