! Written by Ting-Hsuan Hsu

! Start date 1 20 2018
! End date 1 23 2018

! 
!------------------------------------------------------------------------------------------------------------------------------
program Main
    implicit NONE
    
    integer :: G,printout,Nz1,Nr1,Nz2,Nr2
    real*8 :: AA,BB,Re,beta,del_r,del_z,w,conv,lu,ld,Time_Exe
    integer Time_Start,Time_End,Rate

    character :: again,method
    
    again = 'y'
    do while (again == 'y')
        
        call User(AA,BB,Re,G,printout,beta) !User input
        call Initial(AA,BB,Nz1,Nr1,Nz2,Nr2,w,conv,del_r,del_z,lu,ld)    !Initial condition
        call system_clock(Time_Start,Rate)  !Time calculation start
        call Calculate(Nz1,Nr1,Nz2,Nr2,AA,BB,Re,del_z,del_r,w,conv,lu,ld)  !Doing the main iteration
        call system_clock(Time_End,Rate)    !Time calculation end
        Time_Exe = (real(Time_End)-real(Time_Start))/real(Rate)
        write(*,*) '------------------------------------------------------------'
        Write(*,*)'The computational time(sec):', Time_Exe
        write(*,*) '------------------------------------------------------------'

        write(*,*)'Do you want to do it again?(y/n)'
        read(*,*)again
    end do
    
end program Main