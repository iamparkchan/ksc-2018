program main
  use mpi
  use param
  use gram_schmidt
  
  implicit none
  
  integer :: i, rank, ierr
  real(8), dimension(:,:), allocatable :: vector, basis
  real(8) :: time_interval
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  
  if (0 == rank) then
    allocate(vector(0:M-1,0:N-1))
    allocate(basis (0:M-1,0:N-1))
    call initialize(vector)
    basis = vector
    
    time_interval = -MPI_WTIME()
  endif
  
  call gram_schmidt_process(basis)
  
  if (0 == rank) then
    time_interval = time_interval + MPI_WTIME()
    write(*,'(A,f10.6,A)') "Elapsed Time: ", time_interval, " sec"
    write(*,'(A,e13.6)') "Orthogonality Test: ", orthogonality_test(basis)
    write(*,'(A,e13.6)') "Span Test: ", span_test(vector, basis)

    deallocate(vector)
    deallocate(basis)
  endif
  
  
  call mpi_finalize(ierr)
  
contains
  
  subroutine initialize(vector)
    real(8), dimension(0:M-1,0:N-1) :: vector
    integer :: i, k, seed_size, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    CALL system_clock(count=clock)
    seed = clock + 37*(/ (i - 1, i = 1, seed_size) /)
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(vector)
    vector = 2.0_8*vector - 1.0_8
  end subroutine initialize

  real(8) function orthogonality_test(basis)
    real(8), dimension(0:M-1,0:N-1) :: basis
    real(8) :: max, sum_val
    integer :: i, j, k
    
    max = 0.0_8
    do i = 0, N - 1
      do j = 0, N - 1
        sum_val = sum(basis(:,i)*basis(:,j))
        if (i == j) sum_val = sum_val - 1.0_8
        sum_val = abs(sum_val)
        if (sum_val > max) max = sum_val
      enddo
    enddo
    
    orthogonality_test = max
  end function orthogonality_test
  
  real(8) function span_test(vector, basis)
    real(8), dimension(0:M-1,0:N-1) :: vector, basis
    real(8) :: max, max_, sum_val
    real(8), dimension(:), allocatable :: v
    integer :: i, j
    
    max = 0.0_8
    allocate(v(0:M-1))
    do i = 0, N - 1
      v = vector(:, i)
      do j = 0, N - 1
        sum_val = sum(v*basis(:,j))
        v = v - sum_val*basis(:,j)
      enddo
      max_ = maxval(abs(v))
      if (max < max_) max = max_
    enddo
    
    deallocate(v)
    span_test = max
  end function span_test

end program main
