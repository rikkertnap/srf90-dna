! computes L2norm
 
real(dp) function l2norm (f,n)
      
    use precision_definition  
    implicit none

    integer, intent(in)  :: n 
    real(dp), intent(in) :: f(n)

    integer :: i
    real(dp) :: norm

    norm=0.0d0

    do i=1,n
        norm = norm + f(i)*f(i)
    enddo

    l2norm=dsqrt(norm)
    
    return

end function l2norm
      
