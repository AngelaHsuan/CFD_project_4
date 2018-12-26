!-------------------------------------------------------------------------------
subroutine Coeff_str(Nz,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,str1,str2,str3,vor1,vor2,vor3,C_w_str,C_p_str,C_e_str,RHS,position)
    implicit NONE
    
    integer :: k,j,Nz,Nz1,Nz2,Nr1,Nr2
    real*8 :: AA,BB,del_r,del_z,lu
    real*8,dimension(Nz) :: C_w_str,C_p_str,C_e_str,RHS
    real*8, dimension(Nz1,Nr1) :: str1,vor1
    real*8, dimension(Nz2,Nr2) :: str2,str3,vor2,vor3
    character(len=20) :: position
    
    C_w_str(k) = AA**2/BB**2/del_z**2
    C_e_str(k) = AA**2/BB**2/del_z**2
    
    if (position == 'B1_Top2') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1))**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1))/del_r)*str2(k - lu/del_z,1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1))/del_r)*str1(k,j-1) - vor1(k,j)*del_r*(j-1)/AA
    else if (position == 'B1_Top3') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1))**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1))/del_r)*str3(k - (lu+2)/del_z,1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1))/del_r)*str1(k,j-1) - vor1(k,j)*del_r*(j-1)/AA
    else if (position == 'B1_Center') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1))**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1))/del_r)*str1(k,j+1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1))/del_r)*str1(k,j-1) - vor1(k,j)*del_r*(j-1)/AA
    else if (position == 'B2_Bottom') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str2(k,j+1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str1(k + lu/del_z,Nr1) - vor2(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    else if (position == 'B2_Center') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str2(k,j+1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str2(k,j-1) - vor2(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    else if (position == 'B3_Bottom') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str3(k,j+1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str1(k + (lu+2)/del_z,Nr1) - vor2(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    else if (position == 'B3_Center') then
        C_p_str(k) = -2.d0*AA**2/del_r**2 - 2.d0*AA**2/BB**2/del_z**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)**2
        RHS(k) = -( AA**2/del_r**2 - 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str3(k,j+1)&
                & - (AA**2/del_r**2 + 0.5d0*AA**2/(del_r*(j-1) + AA/2.d0 + del_r)/del_r)*str3(k,j-1) - vor3(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    else
        write(*,*)'Mistake in Coeff_str'
        stop
    end if
    
    !if (position == 'B1_Top2') then
    !    RHS(k) = -( AA**2/del_r**2)*str2(k - lu/del_z,1)&
    !            & - (AA**2/del_r**2)*str1(k,j-1) - vor1(k,j)*del_r*(j-1)/AA
    !else if (position == 'B1_Top3') then
    !    RHS(k) = -( AA**2/del_r**2)*str3(k - (lu+2)/del_z,1)&
    !            & - (AA**2/del_r**2)*str1(k,j-1) - vor1(k,j)*del_r*(j-1)/AA
    !else if (position == 'B1_Center') then
    !    RHS(k) = -( AA**2/del_r**2)*str1(k,j+1)&
    !            & - (AA**2/del_r**2)*str1(k,j-1) - vor1(k,j)*del_r*(j-1)/AA
    !else if (position == 'B2_Bottom') then
    !    RHS(k) = -( AA**2/del_r**2)*str2(k,j+1)&
    !            & - (AA**2/del_r**2)*str1(k + lu/del_z,Nr1) - vor2(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    !else if (position == 'B2_Center') then
    !    RHS(k) = -( AA**2/del_r**2)*str2(k,j+1)&
    !            & - (AA**2/del_r**2)*str2(k,j-1) - vor2(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    !else if (position == 'B3_Bottom') then
    !    RHS(k) = -( AA**2/del_r**2)*str3(k,j+1)&
    !            & - (AA**2/del_r**2)*str1(k + (lu+2)/del_z,Nr1) - vor3(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    !else if (position == 'B3_Center') then
    !    RHS(k) = -( AA**2/del_r**2)*str3(k,j+1)&
    !            & - (AA**2/del_r**2)*str3(k,j-1) - vor3(k,j)*(del_r*(j-1) + AA/2.d0 + del_r)/AA
    !else
    !    write(*,*)'Mistake in Coeff_str'
    !    stop
    !end if
    
end subroutine Coeff_str