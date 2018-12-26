!-------------------------------------------------------------------
subroutine vel_half(Nz,Nr,Nz1,Nr1,Nz2,Nr2,k,j,BB,del_z,lu,Vz_EI_p, Vz_EI_n, Vz_WI_p, Vz_WI_n, Vr_NI_p, Vr_NI_n, Vr_SI_p, Vr_SI_n,Vz_p,Vr_p,Vz1,Vr1,Vz2,Vr2,Vz3,Vr3,Vz,Vr,position)
    implicit NONE
    
    integer :: k,j,Nz,Nr,Nz1,Nr1,Nz2,Nr2
    real*8 :: BB,lu,del_z,Vz_EI_p,Vz_EI_n,Vz_WI_p,Vz_WI_n,Vr_NI_p,Vr_NI_n,Vr_SI_p,Vr_SI_n,Vz_p,Vr_p
    real*8, dimension(Nz,Nr) :: Vz,Vr
    real*8, dimension(Nz1,Nr1) :: Vr1,Vz1
    real*8, dimension(Nz2,Nr2) :: Vr2,Vz2,Vr3,Vz3
    character(len=20) :: position
    
    if (k.ne.Nz) then
        Vz_EI_p = (  0.5d0*( Vz(k+1,j)+Vz(k,j) ) + abs( 0.5d0*( Vz(k+1,j)+Vz(k,j) ) )  )/2
        Vz_EI_n = (  0.5d0*( Vz(k+1,j)+Vz(k,j) ) - abs( 0.5d0*( Vz(k+1,j)+Vz(k,j) ) )  )/2
    end if
    if (k.ne.1) then
        Vz_WI_p = (  0.5d0*( Vz(k-1,j)+Vz(k,j) ) + abs( 0.5d0*( Vz(k-1,j)+Vz(k,j) ) )  )/2
        Vz_WI_n = (  0.5d0*( Vz(k-1,j)+Vz(k,j) ) - abs( 0.5d0*( Vz(k-1,j)+Vz(k,j) ) )  )/2
    end if
    if (j.ne.Nr) then
        Vr_NI_p = (  0.5d0*( Vr(k,j+1)+Vr(k,j) ) + abs( 0.5d0*( Vr(k,j+1)+Vr(k,j) ) )  )/2
        Vr_NI_n = (  0.5d0*( Vr(k,j+1)+Vr(k,j) ) - abs( 0.5d0*( Vr(k,j+1)+Vr(k,j) ) )  )/2
    end if
    if (j.ne.1) then
        Vr_SI_p = (  0.5d0*( Vr(k,j-1)+Vr(k,j) ) + abs( 0.5d0*( Vr(k,j-1)+Vr(k,j) ) )  )/2
        Vr_SI_n = (  0.5d0*( Vr(k,j-1)+Vr(k,j) ) - abs( 0.5d0*( Vr(k,j-1)+Vr(k,j) ) )  )/2
    end if
    if (position == 'B2_Bottom') then
        Vr_SI_p = (  0.5d0*( Vr1(k + lu/del_z,Nr1)+Vr(k,j) ) + abs( 0.5d0*( Vr1(k + lu/del_z,Nr1)+Vr(k,j) ) )  )/2
        Vr_SI_n = (  0.5d0*( Vr1(k + lu/del_z,Nr1)+Vr(k,j) ) - abs( 0.5d0*( Vr1(k + lu/del_z,Nr1)+Vr(k,j) ) )  )/2
    end if
    if (position == 'B3_Bottom') then
        Vr_SI_p = (  0.5d0*( Vr1(k + (lu+2)/del_z,Nr1)+Vr(k,j) ) + abs( 0.5d0*( Vr1(k + (lu+2)/del_z,Nr1)+Vr(k,j) ) )  )/2
        Vr_SI_n = (  0.5d0*( Vr1(k + (lu+2)/del_z,Nr1)+Vr(k,j) ) - abs( 0.5d0*( Vr1(k + (lu+2)/del_z,Nr1)+Vr(k,j) ) )  )/2
    end if
    if (position == 'B1_Top2') then
        Vr_NI_p = (  0.5d0*( Vr2(k - lu/del_z,1)+Vr(k,j) ) + abs( 0.5d0*( Vr2(k - lu/del_z,1)+Vr(k,j) ) )  )/2
        Vr_NI_n = (  0.5d0*( Vr2(k - lu/del_z,1)+Vr(k,j) ) - abs( 0.5d0*( Vr2(k - lu/del_z,1)+Vr(k,j) ) )  )/2
    end if
    if (position == 'B1_Top3') then
        Vr_NI_p = (  0.5d0*( Vr3(k - (lu+2)/del_z,1)+Vr(k,j) ) + abs( 0.5d0*( Vr3(k - (lu+2)/del_z,1)+Vr(k,j) ) )  )/2
        Vr_NI_n = (  0.5d0*( Vr3(k - (lu+2)/del_z,1)+Vr(k,j) ) - abs( 0.5d0*( Vr3(k - (lu+2)/del_z,1)+Vr(k,j) ) )  )/2
    end if
    Vz_p = Vz(k,j)
    Vr_p = Vr(k,j)
    
end subroutine vel_half