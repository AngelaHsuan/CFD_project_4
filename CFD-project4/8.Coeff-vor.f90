!------------------------------------------------------------------------------------------------------
subroutine Coeff_vor(Nz,Nr,Nz1,Nz2,Nr1,Nr2,k,j,del_r,del_z,AA,BB,lu,Re,str1,str2,str3,vor1,vor2,vor3,Vr,Vz,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,C_w_str,C_p_str,C_e_str,RHS,position)
    implicit NONE
    
    integer :: k,j,Nz,Nr,Nz1,Nz2,Nr1,Nr2
    real*8 :: AA,BB,del_r,del_z,lu,Re
    real*8,dimension(Nz) :: C_w_str,C_p_str,C_e_str,RHS
    real*8, dimension(Nz1,Nr1) :: str1,vor1,Vr1,Vz1
    real*8, dimension(Nz2,Nr2) :: str2,str3,vor2,vor3,Vr2,Vz2,Vr3,Vz3
    character(len=20) :: position
    
    
    real*8 :: Vz_EI_p, Vz_EI_n, Vz_WI_p, Vz_WI_n, Vr_NI_p, Vr_NI_n, Vr_SI_p, Vr_SI_n,Vz_p,Vr_p
    real*8, dimension(Nz,Nr) :: Vz,Vr
    
    call vel_half(Nz,Nr,Nz1,Nr1,Nz2,Nr2,k,j,BB,del_z,lu,Vz_EI_p, Vz_EI_n, Vz_WI_p, Vz_WI_n, Vr_NI_p, Vr_NI_n, Vr_SI_p, Vr_SI_n,Vz_p,Vr_p,Vz1,Vr1,Vz2,Vr2,Vz3,Vr3,Vz,Vr,position)
    
    
    if (position == 'B1_Top2') then
        C_w_str(k) = -del_r*(j-1)*Vz_WI_p/AA/BB/del_z - del_r*(j-1)/BB**2/Re/del_z**2
        C_p_str(k) = (del_r*(j-0.5d0)*Vr_NI_p - del_r*(j-1.5d0)*Vr_SI_n)/AA/del_r + 2.d0*del_r*(j-1)/del_r**2/Re&
                        & + del_r*(j-1)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*del_r*(j-1)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1))
        C_e_str(k) = del_r*(j-1)*Vz_EI_n/AA/BB/del_z - del_r*(j-1)/BB**2/Re/del_z**2
        RHS(k) = -( del_r*(j-0.5d0)*Vr_NI_n/AA/del_r - del_r*(j-0.5d0)/Re/del_r**2 )*vor2(k - lu/del_z,1)&
                & - ( -del_r*(j-1.5d0)*Vr_SI_p/AA/del_r - del_r*(j-1.5d0)/Re/del_r**2 )*vor1(k,j-1)
    else if (position == 'B1_Top3') then
        C_w_str(k) = -del_r*(j-1)*Vz_WI_p/AA/BB/del_z - del_r*(j-1)/BB**2/Re/del_z**2
        C_p_str(k) = (del_r*(j-0.5d0)*Vr_NI_p - del_r*(j-1.5d0)*Vr_SI_n)/AA/del_r + 2.d0*del_r*(j-1)/del_r**2/Re&
                        & + del_r*(j-1)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*del_r*(j-1)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1))
        C_e_str(k) = del_r*(j-1)*Vz_EI_n/AA/BB/del_z - del_r*(j-1)/BB**2/Re/del_z**2
        RHS(k) = -( del_r*(j-0.5d0)*Vr_NI_n/AA/del_r - del_r*(j-0.5d0)/Re/del_r**2 )*vor3(k - (lu+2)/del_z,1)&
                & - ( -del_r*(j-1.5d0)*Vr_SI_p/AA/del_r - del_r*(j-1.5d0)/Re/del_r**2 )*vor1(k,j-1)
    else if (position == 'B1_Center') then
        C_w_str(k) = -del_r*(j-1)*Vz_WI_p/AA/BB/del_z - del_r*(j-1)/BB**2/Re/del_z**2
        C_p_str(k) = (del_r*(j-0.5d0)*Vr_NI_p - del_r*(j-1.5d0)*Vr_SI_n)/AA/del_r + 2.d0*del_r*(j-1)/del_r**2/Re&
                        & + del_r*(j-1)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*del_r*(j-1)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1))
        C_e_str(k) = del_r*(j-1)*Vz_EI_n/AA/BB/del_z - del_r*(j-1)/BB**2/Re/del_z**2
        RHS(k) = -( del_r*(j-0.5d0)*Vr_NI_n/AA/del_r - del_r*(j-0.5d0)/Re/del_r**2 )*vor1(k,j+1)&
                & - ( -del_r*(j-1.5d0)*Vr_SI_p/AA/del_r - del_r*(j-1.5d0)/Re/del_r**2 )*vor1(k,j-1)
    else if (position == 'B2_Bottom') then
        C_w_str(k) = -(del_r*(j-1) + AA/2.d0 + del_r)*Vz_WI_p/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        C_p_str(k) = ((del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_p - (del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_n)/AA/del_r + 2.d0*(del_r*(j-1) + AA/2.d0 + del_r)/del_r**2/Re&
                        & + (del_r*(j-1) + AA/2.d0 + del_r)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*(del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1) + AA/2.d0 + del_r)
        C_e_str(k) = (del_r*(j-1) + AA/2.d0 + del_r)*Vz_EI_n/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        RHS(k) = -( (del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_n/AA/del_r - (del_r*(j-0.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor2(k,j+1)&
                & - ( -(del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_p/AA/del_r - (del_r*(j-1.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor1(k + lu/del_z,Nr1)
    else if (position == 'B2_Center') then
        C_w_str(k) = -(del_r*(j-1) + AA/2.d0 + del_r)*Vz_WI_p/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        C_p_str(k) = ((del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_p - (del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_n)/AA/del_r + 2.d0*(del_r*(j-1) + AA/2.d0 + del_r)/del_r**2/Re&
                        & + (del_r*(j-1) + AA/2.d0 + del_r)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*(del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1) + AA/2.d0 + del_r)
        C_e_str(k) = (del_r*(j-1) + AA/2.d0 + del_r)*Vz_EI_n/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        RHS(k) = -( (del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_n/AA/del_r - (del_r*(j-0.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor2(k,j+1)&
                & - ( -(del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_p/AA/del_r - (del_r*(j-1.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor2(k,j-1)
    else if (position == 'B3_Bottom') then
        C_w_str(k) = -(del_r*(j-1) + AA/2.d0 + del_r)*Vz_WI_p/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        C_p_str(k) = ((del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_p - (del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_n)/AA/del_r + 2.d0*(del_r*(j-1) + AA/2.d0 + del_r)/del_r**2/Re&
                        & + (del_r*(j-1) + AA/2.d0 + del_r)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*(del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1) + AA/2.d0 + del_r)
        C_e_str(k) = (del_r*(j-1) + AA/2.d0 + del_r)*Vz_EI_n/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        RHS(k) = -( (del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_n/AA/del_r - (del_r*(j-0.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor3(k,j+1)&
                & - ( -(del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_p/AA/del_r - (del_r*(j-1.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor1(k + (lu+2)/del_z,Nr1)
    else if (position == 'B3_Center') then
        C_w_str(k) = -(del_r*(j-1) + AA/2.d0 + del_r)*Vz_WI_p/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        C_p_str(k) = ((del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_p - (del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_n)/AA/del_r + 2.d0*(del_r*(j-1) + AA/2.d0 + del_r)/del_r**2/Re&
                        & + (del_r*(j-1) + AA/2.d0 + del_r)*(Vz_EI_p - Vz_WI_n)/AA/BB/del_z + 2*(del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2 + 1.d0/Re/(del_r*(j-1) + AA/2.d0 + del_r)
        C_e_str(k) = (del_r*(j-1) + AA/2.d0 + del_r)*Vz_EI_n/AA/BB/del_z - (del_r*(j-1) + AA/2.d0 + del_r)/BB**2/Re/del_z**2
        RHS(k) = -( (del_r*(j-0.5d0) + AA/2.d0 + del_r)*Vr_NI_n/AA/del_r - (del_r*(j-0.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor3(k,j+1)&
                & - ( -(del_r*(j-1.5d0) + AA/2.d0 + del_r)*Vr_SI_p/AA/del_r - (del_r*(j-1.5d0) + AA/2.d0 + del_r)/Re/del_r**2 )*vor3(k,j-1)
    else
        write(*,*)'Mistake in Coeff_vor'
        stop
    end if
    
    !write(*,*)C_w_str(k),C_p_str(k),C_e_str(k),RHS(k)
    !pause
    
end subroutine Coeff_vor