!--------------------------------------------------------------------------
subroutine Initial(AA,BB,Nz1,Nr1,Nz2,Nr2,w,conv,del_r,del_z,lu,ld)
    implicit NONE
    
    
    integer :: Nz1,Nr1,Nz2,Nr2
    real*8 :: del_r,del_z,lu,ld,AA,BB,w,conv
    
    !Geometry
    lu = 2.d0
    ld = 0.5d0
    
    del_r = 0.05d0
    del_z = 0.05d0
    
    !Node number
    Nz1 = (lu + 3 + ld)/del_z + 1
    Nr1 = AA/2.d0/del_r + 1
    
    Nz2 = 1/del_z + 1
    Nr2 = 1/del_r
    
    !Iteration
    w = 0.4d0
    conv = 1D-6
    
    
end subroutine Initial