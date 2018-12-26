!------------------------------------------------------------------------------------------------
subroutine Stream(Nz1,Nr1,Nz2,Nr2,AA,BB,lu,ld,w,del_r,del_z,str1,str2,str3,vor1,vor2,vor3,str_err)
    implicit NONE
    
    integer :: Nz1,Nr1,Nz2,Nr2,Pad
    real*8 :: del_r,del_z,lu,ld,AA,BB,w,str_err
    real*8, dimension(Nz1,Nr1) :: str1,vor1
    real*8, dimension(Nz2,Nr2) :: str2,str3,vor2,vor3
    
    integer :: j,k,N_iter
    real*8, dimension(Nz1) :: C_w_str1,C_p_str1,C_e_str1,RHS1,str1_j
    real*8, dimension(Nz2) :: C_w_str2,C_p_str2,C_e_str2,RHS2,str2_j
    real*8, dimension(Nz1,Nr1) :: str1_new
    real*8, dimension(Nz2,Nr2) :: str2_new
    character(len=20) :: position
    
    N_iter = 0
    str_err = 0.d0
    
    do while (N_iter < 1)
        
        !Center Block-------------------------------------------------------------------------
        do j = 1,Nr1,1
            if (j == 1) then
            !Axis:str = 0
                C_w_str1 = 0.d0
                C_p_str1 = 1.d0
                C_e_str1 = 0.d0
                RHS1 = 0.d0
            else if (j == Nr1) then
            !Top
                do k = 1,Nz1,1
                    if (k > lu/del_z + 1 .AND. k < (lu+1)/del_z + 1) then
                    !first block
                        position = 'B1_Top2'
                        call Coeff_str(Nz1,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str1,C_p_str1,C_e_str1,RHS1,position)
                    else if (k > (lu+2)/del_z + 1 .AND. k < (lu+3)/del_z + 1) then
                    !second block
                        position = 'B1_Top3'
                        call Coeff_str(Nz1,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str1,C_p_str1,C_e_str1,RHS1,position)
                    else
                    !On the wall:str = 1/8
                        C_w_str1(k) = 0.d0
                        C_p_str1(k) = 1.d0
                        C_e_str1(k) = 0.d0
                        RHS1(k) = 1.d0/8.d0
                    end if
                end do  !k
            else
                do k = 1,Nz1,1
                    if (k == 1) then
                    !Left(Inlet)
                        C_w_str1(k) = 0.d0
                        C_p_str1(k) = 1.d0
                        C_e_str1(k) = 0.d0
                        RHS1(k) = -2.d0*(del_r*(j-1))**4/AA**4 + (del_r*(j-1))**2/AA**2
                    else if (k == Nz1) then
                    !Right(Outlet)
                        C_w_str1(k) = 0.d0
                        C_p_str1(k) = 1.d0
                        C_e_str1(k) = 0.d0
                        RHS1(k) = str1(k-1,j)
                    else
                    !Center
                        position = 'B1_Center'
                        call Coeff_str(Nz1,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str1,C_p_str1,C_e_str1,RHS1,position)
                    end if
                end do  !k
            end if
            call TDMA(Nz1,C_w_str1,C_p_str1,C_e_str1,RHS1,str1_j)
            
            do k = 1,Nz1,1
                str1_new(k,j) = str1_j(k)
                if (abs(str1_new(k,j)-str1(k,j)) > str_err) then
                    str_err = abs(str1_new(k,j)-str1(k,j))
                end if
                str1(k,j) = w*str1_new(k,j) + (1-w)*str1(k,j)
            end do
            
        end do  !j
        
        !Left top Block-------------------------------------------------------------------------
        do j = 1,Nr2,1
            if (j == 1) then
            !Bottom
                do k = 1,Nz2,1
                    if (k == 1 .OR. k == Nz2) then
                    !On the wall:str = 1/8
                        C_w_str2(k) = 0.d0
                        C_p_str2(k) = 1.d0
                        C_e_str2(k) = 0.d0
                        RHS2(k) = 1.d0/8.d0
                    else
                    !Bottom
                        position = 'B2_Bottom'
                        call Coeff_str(Nz2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str2,C_p_str2,C_e_str2,RHS2,position)
                    end if
                end do  !k
            else if (j == Nr2) then
            !Top wall:str = 1/8
                C_w_str2 = 0.d0
                C_p_str2 = 1.d0
                C_e_str2 = 0.d0
                RHS2 = 1.d0/8.d0
            else
                do k = 1,Nz2,1
                    if (k == 1 .OR. k == Nz2) then
                    !Left & Right
                        C_w_str2(k) = 0.d0
                        C_p_str2(k) = 1.d0
                        C_e_str2(k) = 0.d0
                        RHS2(k) = 1.d0/8.d0
                    else
                    !Center
                        position = 'B2_Center'
                        call Coeff_str(Nz2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str2,C_p_str2,C_e_str2,RHS2,position)
                    end if
                end do  !k
            end if
            call TDMA(Nz2,C_w_str2,C_p_str2,C_e_str2,RHS2,str2_j)
            
            do k = 1,Nz2,1
                str2_new(k,j) = str2_j(k)
                if (abs(str2_new(k,j)-str2(k,j)) > str_err) then
                    str_err = abs(str2_new(k,j)-str2(k,j))
                end if
                str2(k,j) = w*str2_new(k,j) + (1-w)*str2(k,j)
            end do
            
        end do  !j
        
        !Right top Block------------------------------------------------------------------------
        do j = 1,Nr2,1
            if (j == 1) then
            !Bottom
                do k = 1,Nz2,1
                    if (k == 1) then
                    !On the wall:str = 1/8
                        C_w_str2(k) = 0.d0
                        C_p_str2(k) = 1.d0
                        C_e_str2(k) = 0.d0
                        RHS2(k) = 1.d0/8.d0
                    else if (k == Nz2) then
                    !On the wall:str = 1/8
                        C_w_str2(k) = 0.d0
                        C_p_str2(k) = 1.d0
                        C_e_str2(k) = 0.d0
                        RHS2(k) = 1.d0/8.d0
                    else
                    !Bottom
                        position = 'B3_Bottom'
                        call Coeff_str(Nz2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str2,C_p_str2,C_e_str2,RHS2,position)
                    end if
                end do  !k
            else if (j == Nr2) then
            !Top wall:str = 1/8
                C_w_str2 = 0.d0
                C_p_str2 = 1.d0
                C_e_str2 = 0.d0
                RHS2 = 1.d0/8.d0
            else
                do k = 1,Nz2,1
                    if (k == 1 .OR. k == Nz2) then
                    !Left & Right
                        C_w_str2(k) = 0.d0
                        C_p_str2(k) = 1.d0
                        C_e_str2(k) = 0.d0
                        RHS2(k) = 1.d0/8.d0
                    else
                    !Center
                        position = 'B3_Center'
                        call Coeff_str(Nz2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str2,C_p_str2,C_e_str2,RHS2,position)
                    end if
                end do  !k
            end if
            call TDMA(Nz2,C_w_str2,C_p_str2,C_e_str2,RHS2,str2_j)
            
            do k = 1,Nz2,1
                str2_new(k,j) = str2_j(k)
                if (abs(str2_new(k,j)-str3(k,j)) > str_err) then
                        str_err = abs(str2_new(k,j)-str3(k,j))
                end if
                str3(k,j) = w*str2_new(k,j) + (1-w)*str3(k,j)
            end do
            
        end do  !j
        
       N_iter = N_iter + 1 
        
    end do  !while
    
    !Pad = 100
    !
    !open (unit=Pad, file="stream.txt", status="UNKNOWN")
    !write(Pad,*) 'variables="z","r","stream"'
    !write(Pad,*) 'zone i=', Nz1,'j=',Nr1+Nr2,'DATAPACKING=POINT'
    !do j = 1,Nr1+Nr2,1
    !    do k = 1,Nz1,1
    !        if (j < Nr1+1) then
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str1(k,j)
    !        else if (j > Nr1 .AND. k > lu/del_z .AND. k < (lu+1)/del_z + 2) then
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str2(k - lu/del_z,j - Nr1)
    !        else if (j > Nr1 .AND. k > (lu+2)/del_z .AND. k < (lu+3)/del_z + 2) then
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str3(k - (lu+2)/del_z,j - Nr1)
    !        else
    !            write(Pad,'(2F30.12)')(k-1)*del_z,(j-1)*del_r,0
    !        end if
    !    end do
    !end do
    !close (Pad,status = 'Keep')
    

end subroutine Stream