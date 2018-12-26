!--------------------------------------------------------------------------------
subroutine Calculate(Nz1,Nr1,Nz2,Nr2,AA,BB,Re,del_z,del_r,w,conv,lu,ld)
    implicit NONE
    
    integer :: Nz1,Nr1,Nz2,Nr2,Pad
    real*8 :: AA,BB,Re,del_z,del_r,w,conv,lu,ld
    
    integer :: j,k,N_iter
    real*8 :: str_err,vor_err
    real*8, dimension(Nz1,Nr1) :: str1,vor1,Vr1,Vz1
    real*8, dimension(Nz2,Nr2) :: str2,str3,vor2,vor3,Vr2,Vz2,Vr3,Vz3
    
    !Initial guess
    !do i = 1,500
    !    !Stream function
    !    call Stream(Nz1,Nr1,Nz2,Nr2,AA,BB,lu,ld,w,del_r,del_z,str1,str2,str3,vor1,vor2,vor3,str_err)
    !end do
    str1 = 0.d0
    str2 = 0.d0
    str3 = 0.d0
    vor1 = 0.d0
    vor2 = 0.d0
    vor3 = 0.d0
    
    str_err = 1
    vor_err = 1
    
    N_iter = 0
    
    !Iteration
    do while (str_err > conv .AND. vor_err > conv)
        !Stream function
        call Stream(Nz1,Nr1,Nz2,Nr2,AA,BB,lu,ld,w,del_r,del_z,str1,str2,str3,vor1,vor2,vor3,str_err)
        write(*,'(A22,I7,A11,1E12.5)') 'number of iteration :', N_iter,'error=', str_err
        
        call vel(Nz1,Nr1,Nz2,Nr2,AA,BB,del_z,del_r,lu,ld,str1,str2,str3,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3)
           
        !Vorticity
        call Vorticity(Nz1,Nr1,Nz2,Nr2,AA,BB,Re,lu,ld,w,del_r,del_z,vor1,vor2,vor3,str1,str2,str3,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3,vor_err)
        
        !Error
        write(*,'(A22,I7,A11,1E12.5)') 'number of iteration :', N_iter,'error=', vor_err
        N_iter = N_iter + 1
        
        if (vor_err > 1D100 .OR. str_err > 1D100) then
            write(*,*)'---------------------------------------------------------------'
            write(*,*)'The result divereged!!!!!!!!!!!'
            str_err = 0
            vor_err = 0
        end if
        
    end do  !while
    
    
    !-------------------------------------------------------------------------------
    
    Pad = 100
    open (unit=Pad, file="Vel.txt", status="UNKNOWN")
    write(Pad,*) 'variables="z","r","Vz","Vr"'
    write(Pad,*) 'zone i=', Nz1,'j=',Nr1+Nr2,'DATAPACKING=POINT'
    do j = 1,Nr1+Nr2,1
        do k = 1,Nz1,1
            if (j < Nr1+1) then
                write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,Vz1(k,j),Vr1(k,j)
            else if (j > Nr1 .AND. k > lu/del_z .AND. k < (lu+1)/del_z + 2) then
                write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,Vz2(k - lu/del_z,j - Nr1),Vr2(k - lu/del_z,j - Nr1)
            else if (j > Nr1 .AND. k > (lu+2)/del_z .AND. k < (lu+3)/del_z + 2) then
                write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,Vz3(k - (lu+2)/del_z,j - Nr1),Vr3(k - (lu+2)/del_z,j - Nr1)
            else
                write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,0,0
            end if
        end do
    end do
    close (Pad,status = 'Keep')
    
    
    
    Pad = 100
    
    open (unit=Pad, file="stream.txt", status="UNKNOWN")
    write(Pad,*) 'variables="z","r","stream"'
    write(Pad,*) 'zone i=', Nz1,'j=',Nr1+Nr2,'DATAPACKING=POINT'
    do j = 1,Nr1+Nr2,1
        do k = 1,Nz1,1
            if (j < Nr1+1) then
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str1(k,j)
            else if (j > Nr1 .AND. k > lu/del_z .AND. k < (lu+1)/del_z + 2) then
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str2(k - lu/del_z,j - Nr1)
            else if (j > Nr1 .AND. k > (lu+2)/del_z .AND. k < (lu+3)/del_z + 2) then
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str3(k - (lu+2)/del_z,j - Nr1)
            else
                write(Pad,'(2F30.12)')(k-1)*del_z,(j-1)*del_r,0
            end if
        end do
    end do
    close (Pad,status = 'Keep')
    !-------------
    open (unit=Pad, file="vorticity.txt", status="UNKNOWN")
    write(Pad,*) 'variables="z","r","vorticity"'
    write(Pad,*) 'zone i=', Nz1,'j=',Nr1+Nr2,'DATAPACKING=POINT'
    do j = 1,Nr1+Nr2,1
        do k = 1,Nz1,1
            if (j < Nr1+1) then
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,vor1(k,j)
            else if (j > Nr1 .AND. k > lu/del_z .AND. k < (lu+1)/del_z + 2) then
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,vor2(k - lu/del_z,j - Nr1)
            else if (j > Nr1 .AND. k > (lu+2)/del_z .AND. k < (lu+3)/del_z + 2) then
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,vor3(k - (lu+2)/del_z,j - Nr1)
            else
                write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,0
            end if
        end do
    end do
    close (Pad,status = 'Keep')
    
    Pad = 100
    
    open (unit=Pad, file="stream_c1.txt", status="UNKNOWN")
    
    write(Pad,*) 'variables="z","r","stream_c1"'
    write(Pad,*) 'zone i=', Nz2,'j=',Nr2,'DATAPACKING=POINT'
    do j = 1,Nr2,1
        do k = 1,Nz2,1
            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str2(k,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
    Pad = 100
    
    open (unit=Pad, file="stream_c2.txt", status="UNKNOWN")
    
    write(Pad,*) 'variables="z","r","stream_c2"'
    write(Pad,*) 'zone i=', Nz2,'j=',Nr2,'DATAPACKING=POINT'
    do j = 1,Nr2,1
        do k = 1,Nz2,1
            write(Pad,'(3F30.12)')(k-1)*del_z,(j-1)*del_r,str3(k,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine Calculate