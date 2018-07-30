module gram_schmidt
  use mpi
  use param

  implicit none

  private :: dot_product, multiply_add, multiply, ceil_div
  
contains
  
  real(8) function dot_product(A, B, n, comm)
    real(8), dimension(0:n-1) :: A, B
    integer :: n, comm
    integer :: ierr
    integer :: k
    real(8) :: sum

    sum = 0.0_8
    !$omp parallel do reduction(+:sum)
    do k = 0, n - 1
      sum = sum + A(k)*B(k)
    enddo
    !$omp end parallel do
    
    call mpi_allreduce(mpi_in_place, sum, 1, mpi_real8, mpi_sum, comm, ierr)
    
    dot_product = sum
  end function dot_product
  
  subroutine multiply_add(A, c, B, n)
    real(8), dimension(0:n-1) :: A, B
    real(8) :: c
    integer :: n
    integer :: k
    
    !$omp parallel do
    do k = 0, n - 1
      A(k) = A(k) + c*B(k)
    enddo
    !$omp end parallel do
  end subroutine multiply_add
  
  subroutine multiply(A, c, n)
    real(8), dimension(0:n-1) :: A
    real(8) :: c
    integer :: n
    integer :: k
    
    !$omp parallel do
    do k = 0, n - 1
      A(k) = A(k)*c
    enddo
    !$omp end parallel do
  end subroutine multiply
  
  integer function ceil_div(x, y)
    integer :: x, y

    ceil_div = ((x) + (y) - 1) / (y)
  end function ceil_div
  
  subroutine gram_schmidt_process(vector)
    real(8), dimension(0:M-1,0:N-1) :: vector
    real(8), dimension(:,:), allocatable :: v
    integer :: i, j
    real(8) :: coef
    integer :: ierr, rank, size, unit, nprocs, world_group, work_group, range(3,1), comm, r, chunk
    integer, dimension(:), allocatable :: counts, displs
    
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, size, ierr)
    
    unit = ceil_div(M, size)
    nprocs = ceil_div(M, unit)
    
    call mpi_comm_group(mpi_comm_world, world_group, ierr)
    range(:,1) = (/ 0, nprocs - 1, 1 /)
    call mpi_group_range_incl(world_group, 1, range, work_group, ierr)
    call mpi_comm_create(mpi_comm_world, work_group, comm, ierr)
    call mpi_group_free(world_group, ierr)
    call mpi_group_free(work_group, ierr)
    
    if (rank >= nprocs) return
    call mpi_comm_rank(comm, rank, ierr)

    allocate(counts(0:nprocs-1))
    allocate(displs(0:nprocs-1))

    do r = 0, nprocs - 1
      if (r /= nprocs - 1) then
        counts(r) = unit
      else
        counts(r) = M - unit*r
      endif
      displs(r) = unit*r
    enddo
    
    chunk = counts(rank)
    
    allocate(v(0:chunk-1,0:N-1))

    do i = 0, N - 1
      call mpi_scatterv(vector(:,i), counts, displs, mpi_real8, v(:,i), chunk, mpi_real8, 0, comm, ierr)
    enddo
    
    do j = 0, N - 1
      do i = 0, j - 1
        coef = -dot_product(v(:,i), v(:,j), chunk, comm)
        call multiply_add(v(:,j), coef, v(:,i), chunk)
      enddo
      coef = 1.0_8 / sqrt(dot_product(v(:,j), v(:,j), chunk, comm))
      call multiply(v(:,j), coef, chunk)
    enddo

    do i = 0, N - 1
      call mpi_gatherv(v(:,i), chunk, mpi_real8, vector(:,i), counts, displs, mpi_real8, 0, comm, ierr)
    enddo
    
    deallocate(v)
    deallocate(counts)
    deallocate(displs)
    call mpi_comm_free(comm, ierr)
  end subroutine gram_schmidt_process
  
end module gram_schmidt
