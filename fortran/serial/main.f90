program main
  use mpi
  use gram_schmidt
  
  implicit none
  
  integer :: i, rank, ierr
  real(8),dimension(:,:),allocatable :: vector
  real(8) :: time_interval = 0.0_8, rms
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  
  if (0 == rank) then
    allocate(vector(0:N-1,0:M-1))
    call initialize(vector)
    time_interval = -MPI_WTIME()
  endif
  
  call gram_schmidt_process(vector)
  
  if (0 == rank) then
    time_interval = time_interval + MPI_WTIME()
    rms = rms_error(vector)
    write(*,'(A,e20.10)') "error: ", rms
    write(*,'(A,f10.3,A)') "elapsed time: ", time_interval, " sec"
  endif
  
  deallocate(vector)
  
  call MPI_FINALIZE(ierr)
  
contains
  
  subroutine initialize(vector)
    real(8), dimension(0:N-1,0:M-1) :: vector
    integer :: i, k
    integer, dimension(3) :: timeArray
    
    ! call itime(timeArray)
    ! call srand(timeArray(1)+timeArray(2)+timeArray(3))
    call srand(0)
    do i = 0, M - 1
      do k = 0, N - 1
        vector(k, i) = 2.0_8*rand() - 1.0_8
      enddo
    enddo

  end subroutine initialize
  
  real(8) function dot_product(A, B)
    real(8), dimension(0:N-1) :: A, B
    integer :: k
    real(8) :: sum
      
    sum = 0.0_8
    do k = 0, N - 1
      sum = sum + A(k)*B(k)
    enddo
    
    dot_product = sum
  end function dot_product
  
  
  real(8) function rms_error(vector)
    real(8), dimension(0:N-1,0:M-1) :: vector
    integer :: i, j
    real(8) :: sum = 0.0_8, err
    
    do i = 0, M - 1
      do j = 0, M - 1
        err = dot_product(vector(:,i), vector(:,j))
        if (i == j) err = err - 1.0_8
        err = err/N
        sum = sum + err*err
      enddo
    enddo    
    rms_error = SQRT(sum/M/M)
  end function rms_error
  
end program main
