!---------------------------------------------------------
subroutine TDMA(N,Mat_a,Mat_b,Mat_c,RHS_s,x)
    implicit NONE
    real*8, dimension(N) :: Mat_a,Mat_b,Mat_c,RHS_s,x
    integer :: N
    
    integer :: i
    real*8, dimension(N) :: TDMA_r,TDMA_b,TDMA_s
       
    !------------------start TDMA--------------------
    TDMA_b(1) = Mat_b(1)
    TDMA_s(1) = RHS_s(1)
    
    do i=2, N, 1
        TDMA_r(i) = Mat_a(i)/TDMA_b(i-1)
        TDMA_b(i) = Mat_b(i) - (TDMA_r(i)*Mat_c(i-1))
        TDMA_s(i) = RHS_s(i) - (TDMA_r(i)*TDMA_s(i-1))
    
    end do
    
    
    
    x(N) =  TDMA_s(N)/TDMA_b(N)
    
    
    do i=N-1,1,-1
        x(i)=(TDMA_s(i)-Mat_c(i)*x(i+1))/TDMA_b(i)
        
    end do
end subroutine TDMA