!--------------------------------------------------------------------------------------------------------------------
subroutine Vorticity(Nz1,Nr1,Nz2,Nr2,AA,BB,Re,lu,ld,w,del_r,del_z,vor1,vor2,vor3,str1,str2,str3,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,vor_err)
    implicit NONE
    
    integer :: Nz1,Nr1,Nz2,Nr2,Pad
    real*8 :: del_r,del_z,lu,ld,AA,BB,w,Re,vor_err
    real*8, dimension(Nz1,Nr1) :: vor1,str1,Vr1,Vz1
    real*8, dimension(Nz2,Nr2) :: vor2,vor3,str2,str3,Vr2,Vz2,Vr3,Vz3
    
    integer :: j,k,N_iter
    real*8, dimension(Nz1) :: C_w_vor1,C_p_vor1,C_e_vor1,RHS1,vor1_j
    real*8, dimension(Nz2) :: C_w_vor2,C_p_vor2,C_e_vor2,RHS2,vor2_j
    real*8, dimension(Nz1,Nr1) :: vor1_new
    real*8, dimension(Nz2,Nr2) :: vor2_new
    character(len=20) :: position
    
    N_iter = 0
    vor_err = 0.d0
    
    do while (N_iter < 1)
        
        !Center Block-------------------------------------------------------------------------
        do j = 1,Nr1,1
            if (j == 1) then
            !Axis:vor = 0
                C_w_vor1 = 0.d0
                C_p_vor1 = 1.d0
                C_e_vor1 = 0.d0
                RHS1 = 0.d0
            else if (j == Nr1) then
            !Top
                do k = 1,Nz1,1
                    if (k == lu/del_z + 1 .OR. k == (lu+2)/del_z + 1 ) then
                    !Corner_C
                        C_w_vor1(k) = 0.d0
                        C_p_vor1(k) = 1.d0
                        C_e_vor1(k) = 0.d0
                        RHS1(k) = AA**3/(del_r*(j-1))*((str1(k,j) - str1(k,j-1))/del_r**2 + 1.d0/BB**2*(str1(k,j) - str1(k+1,j))/del_z**2)
                    else if (k == (lu+1)/del_z + 1 .OR. k == (lu+3)/del_z + 1 ) then
                    !Corner_F
                        C_w_vor1(k) = 0.d0
                        C_p_vor1(k) = 1.d0
                        C_e_vor1(k) = 0.d0
                        RHS1(k) = AA**3/(del_r*(j-1))*((str1(k,j) - str1(k,j-1))/del_r**2 + 1.d0/BB**2*(str1(k,j) - str1(k-1,j))/del_z**2)
                    else if (k > lu/del_z + 1 .AND. k < (lu+1)/del_z + 1) then
                    !first block
                        position = 'B1_Top2'
                        call Coeff_vor(Nz1,Nr1,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr1,Vz1,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor1,C_p_vor1,C_e_vor1,RHS1,position)
                    else if (k > (lu+2)/del_z + 1 .AND. k < (lu+3)/del_z + 1) then
                    !second block
                        position = 'B1_Top3'
                        call Coeff_vor(Nz1,Nr1,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr1,Vz1,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor1,C_p_vor1,C_e_vor1,RHS1,position)
                    else
                    !On the wall_BC.DE.FG
                        C_w_vor1(k) = 0.d0
                        C_p_vor1(k) = 1.d0
                        C_e_vor1(k) = 0.d0
                        RHS1(k) = 2.d0*AA**3/(del_r*(j-1))*(str1(k,j) - str1(k,j-1))/del_r**2
                    end if
                end do  !k
            else
                do k = 1,Nz1,1
                    if (k == 1) then
                    !Left(Inlet)
                        C_w_vor1(k) = 0.d0
                        C_p_vor1(k) = 1.d0
                        C_e_vor1(k) = 0.d0
                        !RHS1(k) = 2.d0*AA**3*(str1(k,j) - str1(k+1,j))/(del_r*(j-1))/BB**2/del_z**2 + 16*del_r*(j-1)/AA
                        RHS1(k) = 2.d0*AA**3*(str1(k,j) - str1(k+1,j))/(del_r*(j-1))/BB**2/del_z**2 + AA**3*(str1(k,j+1) - str1(k,j-1))/(del_r*(j-1))**2/del_r - AA**3*(str1(k,j+1) - 2*str1(k,j) + str1(k,j-1))/(del_r*(j-1))/del_r**2
                    else if (k == Nz1) then
                    !Right(Outlet)
                        C_w_vor1(k) = 0.d0
                        C_p_vor1(k) = 1.d0
                        C_e_vor1(k) = 0.d0
                        RHS1(k) = vor1(k-1,j)
                    else
                    !Center
                        position = 'B1_Center'
                        call Coeff_vor(Nz1,Nr1,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr1,Vz1,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor1,C_p_vor1,C_e_vor1,RHS1,position)
                    end if
                end do  !k
            end if
            !write(*,*)C_p_vor1,j
            !pause
            call TDMA(Nz1,C_w_vor1,C_p_vor1,C_e_vor1,RHS1,vor1_j)
            
            do k = 1,Nz1,1
                vor1_new(k,j) = vor1_j(k)
                if (abs(vor1_new(k,j)-vor1(k,j)) > vor_err) then
                    vor_err = abs(vor1_new(k,j)-vor1(k,j))
                end if
                vor1(k,j) = w*vor1_new(k,j) + (1-w)*vor1(k,j)
            end do
            
        end do  !j
        
        
        
        
        !Left top Block-------------------------------------------------------------------------
        do j = 1,Nr2,1
            if (j == 1) then
            !Bottom
                do k = 1,Nz2,1
                    if (k == 1) then
                    !Wall_CD
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str2(k,j) - str2(k+1,j))/del_z**2
                    else if (k == Nz2) then
                    !Wall_EF
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str2(k,j) - str2(k-1,j))/del_z**2
                    else
                    !Bottom
                        position = 'B2_Bottom'
                        call Coeff_vor(Nz2,Nr2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr2,Vz2,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor2,C_p_vor2,C_e_vor2,RHS2,position)
                    end if
                end do  !k
            else if (j == Nr2) then
            !Top wall
                do k = 1,Nz2,1
                    if (k == 1) then
                    !Corner_D
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = AA**3/(del_r*(j-1) + AA/2.d0 + del_r)*((str2(k,j) - str2(k,j-1))/del_r**2 + 1.d0/BB**2*(str2(k,j) - str2(k+1,j))/del_z**2)
                    else if (k == Nz2) then
                    !Corner_E
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = AA**3/(del_r*(j-1) + AA/2.d0 + del_r)*((str2(k,j) - str2(k,j-1))/del_r**2 + 1.d0/BB**2*(str2(k,j) - str2(k-1,j))/del_z**2)
                    else
                    !Top
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/(del_r*(j-1) + AA/2.d0 + del_r)*(str2(k,j) - str2(k,j-1))/del_r**2
                    end if
                end do
                
                
            else
                do k = 1,Nz2,1
                    if (k == 1) then
                    !Wall_CD
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str2(k,j) - str2(k+1,j))/del_z**2
                    else if (k == Nz2) then
                    !Wall_EF
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str2(k,j) - str2(k-1,j))/del_z**2
                    else
                    !Center
                        position = 'B2_Center'
                        call Coeff_vor(Nz2,Nr2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr2,Vz2,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor2,C_p_vor2,C_e_vor2,RHS2,position)
                    end if
                end do  !k
            end if
            call TDMA(Nz2,C_w_vor2,C_p_vor2,C_e_vor2,RHS2,vor2_j)
            
            do k = 1,Nz2,1
                vor2_new(k,j) = vor2_j(k)
                if (abs(vor2_new(k,j)-vor2(k,j)) > vor_err) then
                    vor_err = abs(vor2_new(k,j)-vor2(k,j))
                end if
                vor2(k,j) = w*vor2_new(k,j) + (1-w)*vor2(k,j)
            end do
            
        end do  !j
        
        !Right top Block------------------------------------------------------------------------
        do j = 1,Nr2,1
            if (j == 1) then
            !Bottom
                do k = 1,Nz2,1
                    if (k == 1) then
                    !Wall_CD
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str3(k,j) - str3(k+1,j))/del_z**2
                    else if (k == Nz2) then
                    !Wall_EF
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str3(k,j) - str3(k-1,j))/del_z**2
                    else
                    !Bottom
                        position = 'B3_Bottom'
                        call Coeff_vor(Nz2,Nr2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr3,Vz3,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor2,C_p_vor2,C_e_vor2,RHS2,position)
                    end if
                end do  !k
            else if (j == Nr2) then
            !Top wall
                do k = 1,Nz2,1
                    if (k == 1) then
                    !Corner_D
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = AA**3/(del_r*(j-1) + AA/2.d0 + del_r)*((str3(k,j) - str3(k,j-1))/del_r**2 + 1.d0/BB**2*(str3(k,j) - str3(k+1,j))/del_z**2)
                    else if (k == Nz2) then
                    !Corner_E
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = AA**3/(del_r*(j-1) + AA/2.d0 + del_r)*((str3(k,j) - str3(k,j-1))/del_r**2 + 1.d0/BB**2*(str3(k,j) - str3(k-1,j))/del_z**2)
                    else
                    !Top
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/(del_r*(j-1) + AA/2.d0 + del_r)*(str3(k,j) - str3(k,j-1))/del_r**2
                    end if
                end do
            else
                do k = 1,Nz2,1
                    if (k == 1) then
                    !Wall_CD
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str3(k,j) - str3(k+1,j))/del_z**2
                    else if (k == Nz2) then
                    !Wall_EF
                        C_w_vor2(k) = 0.d0
                        C_p_vor2(k) = 1.d0
                        C_e_vor2(k) = 0.d0
                        RHS2(k) = 2.d0*AA**3/BB**2/(del_r*(j-1) + AA/2.d0 + del_r)*(str3(k,j) - str3(k-1,j))/del_z**2
                    else
                    !Center
                        position = 'B3_Center'
                        call Coeff_vor(Nz2,Nr2,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr3,Vz3,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_vor2,C_p_vor2,C_e_vor2,RHS2,position)
                    end if
                end do  !k
            end if
            call TDMA(Nz2,C_w_vor2,C_p_vor2,C_e_vor2,RHS2,vor2_j)
            
            do k = 1,Nz2,1
                vor2_new(k,j) = vor2_j(k)
                if (abs(vor2_new(k,j)-vor3(k,j)) > vor_err) then
                    vor_err = abs(vor2_new(k,j)-vor3(k,j))
                end if
                vor3(k,j) = w*vor2_new(k,j) + (1-w)*vor3(k,j)
            end do
        
            N_iter = N_iter + 1 
        
        end do  !j
    end do  !while
    
    !Pad = 100
    !open (unit=Pad, file="vorticity.txt", status="UNKNOWN")
    !write(Pad,*) 'variables="z","r","vorticity"'
    !write(Pad,*) 'zone i=', Nz1,'j=',Nr1+Nr2,'DATAPACKING=POINT'
    !do j = 1,Nr1+Nr2,1
    !    do k = 1,Nz1,1
    !        if (j < Nr1+1) then
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,vor1(k,j)
    !        else if (j > Nr1 .AND. k > lu/del_z .AND. k < (lu+1)/del_z + 2) then
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,vor2(k - lu/del_z,j - Nr1)
    !        else if (j > Nr1 .AND. k > (lu+2)/del_z .AND. k < (lu+3)/del_z + 2) then
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,vor3(k - (lu+2)/del_z,j - Nr1)
    !        else
    !            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,0
    !        end if
    !    end do
    !end do
    !close (Pad,status = 'Keep')
    !pause
    
end subroutine Vorticity