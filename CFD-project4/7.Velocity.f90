!-----------------------------------------------------------------------------------
subroutine vel(Nz1,Nr1,Nz2,Nr2,AA,BB,del_z,del_r,lu,ld,str1,str2,str3,Vr1,Vz1,Vr2,Vz2,Vr3,Vz3)
    implicit NONE
	
	integer :: Nz1,Nr1,Nz2,Nr2,Pad
	real*8 :: AA,BB,del_z,del_r,lu,ld
	real*8, dimension(Nz1,Nr1) :: str1,vor1,Vr1,Vz1
    real*8, dimension(Nz2,Nr2) :: str2,str3,vor2,vor3,Vr2,Vz2,Vr3,Vz3
	
	integer :: j,k
    
    !Calculate velocity field
    do j = 2,Nr1-1
        do k = 2,Nz1-1
            Vr1(k,j) = -0.5d0*AA**2*(str1(k+1,j) - str1(k-1,j))/del_z/BB/(del_r*(j-1))
            Vz1(k,j) = 0.5d0*AA**2*(str1(k,j+1) - str1(k,j-1))/del_r/(del_r*(j-1))
		end do
	end do
	
    do j = 2,Nr2-1
        do k = 2,Nz2-1
            Vr2(k,j) = -0.5d0*AA**2*(str2(k+1,j) - str2(k-1,j))/del_z/BB/(del_r*(j-1) + AA/2.d0 + del_r)
            Vz2(k,j) = 0.5d0*AA**2*(str2(k,j+1) - str2(k,j-1))/del_r/(del_r*(j-1) + AA/2.d0 + del_r)
            Vr3(k,j) = -0.5d0*AA**2*(str3(k+1,j) - str3(k-1,j))/del_z/BB/(del_r*(j-1) + AA/2.d0 + del_r)
            Vz3(k,j) = 0.5d0*AA**2*(str3(k,j+1) - str3(k,j-1))/del_r/(del_r*(j-1) + AA/2.d0 + del_r)
        end do
    end do
    !B1_boundary
    Vr1(:,1) = Vr1(:,2)
    Vz1(:,1) = Vz1(:,2)
    do k = 1,Nz1
        if (k > lu/del_z + 1 .AND. k < (lu+1)/del_z + 1) then
            Vr1(k,Nr1) = Vr1(k,Nr1-1)
            Vz1(k,Nr1) = Vz1(k,Nr1-1)
        else if (k > (lu+2)/del_z + 1 .AND. k < (lu+3)/del_z + 1) then
            Vr1(k,Nr1) = Vr1(k,Nr1-1)
            Vz1(k,Nr1) = Vz1(k,Nr1-1)
        else
            Vr1(k,Nr1) = 0
            Vz1(k,Nr1) = 0
        end if
    end do
    !B2_boundary
    Vz2(1,:) = 0
    Vr2(1,:) = 0
    Vz2(Nz2,:) = 0
    Vr2(Nz2,:) = 0
    Vz2(:,Nr2) = 0
    Vr2(:,Nr2) = 0
    Vz2(:,1) = Vz2(:,2)
    Vr2(:,1) = Vr2(:,2)
    !B3_boundary
    Vz3(1,:) = 0
    Vr3(1,:) = 0
    Vz3(Nz2,:) = 0
    Vr3(Nz2,:) = 0
    Vz3(:,Nr2) = 0
    Vr3(:,Nr2) = 0
    Vz3(:,1) = Vz3(:,2)
    Vr3(:,1) = Vr3(:,2)
    !Inlet velocity
    do j = 1,Nr1
        Vz1(1,j) = -8.d0*(del_r*(j-1))**2/AA**2 + 2.d0
        Vr1(1,j) = 0
    end do
	!Outlet velocity
    Vz1(Nz1,:) = Vz1(Nz1-1,:)
    Vr1(Nz1,:) = Vr1(Nz1-1,:)
	
	
 !  Pad = 100
 !       open (unit=Pad, file="Vr.txt", status="UNKNOWN")
 !       write(Pad,*) 'variables="z","r","Vr"'
 !       write(Pad,*) 'zone i=', Nz1,'j=',Nr1+Nr2,'DATAPACKING=POINT'
 !       do j = 1,Nr1+Nr2,1
 !           do k = 1,Nz1,1
 !               if (j < Nr1+1) then
 !                   write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,Vr1(k,j)
 !               else if (j > Nr1 .AND. k > lu/del_z .AND. k < (lu+1)/del_z + 2) then
 !                   write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,Vr2(k - lu/del_z,j - Nr1)
 !               else if (j > Nr1 .AND. k > (lu+2)/del_z .AND. k < (lu+3)/del_z + 2) then
 !                   write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,Vr3(k - (lu+2)/del_z,j - Nr1)
 !               else
 !                   write(Pad,'(4F25.12)')(k-1)*del_z,(j-1)*del_r,0
 !               end if
 !           end do
 !       end do
 !       close (Pad,status = 'Keep')
 !       pause
    
end subroutine vel