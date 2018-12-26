!------------------------------------------------------------------------
subroutine user(AA,BB,Re,G,printout,beta)
    implicit NONE
    
    integer :: A_select,B_select,Re_select,method_select,G,printout
    real*8 :: AA,BB,Re,beta
    
    !------------------------------------------------
    !D/d
    write(*,*)'Choose the number you want for D/d.'
    write(*,*)'(1) 2'
    write(*,*)'(2) 6'
    write(*,*)'(3) else'
    read(*,*)A_select
    select case (A_select)
        case(1)
            AA = 2.d0
        case(2)
            AA = 6.d0
        case(3)
            read(*,*)AA
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select

    !---------------------------------------------------
    !w/d
    write(*,*)'Choose the number you want for w/d.'
    write(*,*)'(1) 0.1'
    write(*,*)'(2) 0.5'
    write(*,*)'(3) 1'
    write(*,*)'(4) 5'
    write(*,*)'(5) else'
    read(*,*)B_select
    select case (B_select)
        case(1)
            BB = 0.1d0
        case(2)
            BB = 0.5d0
        case(3)
            BB = 1.d0
        case(4)
            BB = 5.d0
        case(5)
            read(*,*)BB
        case default
            write(*,*)'Please follow the instructions~'
            stop
        end select
    
    !---------------------------------------------------
    !Re
    write(*,*)'Choose the number you want for Re.'
    write(*,*)'(1) 100'
    write(*,*)'(2) 500'
    write(*,*)'(3) 1000'
    write(*,*)'(4) 2000'
    write(*,*)'(5) else'
    read(*,*)Re_select
    select case (Re_select)
        case(1)
            Re = 100.d0
        case(2)
            Re = 500.d0
        case(3)
            Re = 1000.d0
        case(4)
            Re = 2000.d0
        case(5)
            read(*,*)Re
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
        
    !-------------------------------------------------
    write(*,*)'Choose a method to solve the problem.'
    write(*,*)'(1) FOU'
    write(*,*)'(2) QUICK'
    read(*,*)method_select
    select case (method_select)
        case(1)
            G = 0
            beta = 0
            printout = 1
        case(2)
            G = 1
            beta = 1.d0/8.d0
            printout = 2
        case default
            write(*,*)'Please follow the instructions~'
            stop
    end select
    
end subroutine user