module gram_schmidt
  use mpi

  implicit none
  
  integer, parameter :: N = 20000
  integer, parameter :: M = 400
  
  private :: dot_product, multiply_add, multiply
  
contains
  
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
  
  subroutine multiply_add(A, c, B)
    real(8), dimension(0:N-1) :: A, B
    real(8) :: c
    integer :: k
    
    do k = 0, N - 1
      A(k) = A(k) + c*B(k)
    enddo
  end subroutine multiply_add
  
  subroutine multiply(A, c)
    real(8), dimension(0:N-1) :: A
    real(8) :: c
    integer :: k
    
    do k = 0, N - 1
      A(k) = A(k)*c
    enddo
  end subroutine multiply
  
  subroutine gram_schmidt_process(vector)
    real(8), dimension(0:N-1,0:M-1) :: vector
    integer :: rank, ierr, i, j
    real(8) :: coef
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    if (0 == rank) then
      do j = 0, M - 1
        do i = 0, j - 1
          coef = -dot_product(vector(:,i), vector(:,j))
          call multiply_add(vector(:,j), coef, vector(:,i))
        enddo
        coef = 1.0_8 / SQRT(dot_product(vector(:,j), vector(:,j)))
        call multiply(vector(:,j), coef)
      enddo
    endif
  end subroutine gram_schmidt_process
  
end module gram_schmidt
