program main_p
use module_Data     
implicit none
integer n_atoms
integer N !  N values; N-AGNR-I(n,m)
integer n1,n2,n3
integer t1 ,t2,re
real t_ini,t_end

do N = 7,7,2

n3=(N+4)
n2=(N-2)+4
n1=(N-2)+2
!-------------------------------------------------------
if (m==2) then
n_atoms = n_lowercase*2*n3+2*n2
call tight_binding1(n,n_atoms,n1,n2,n3) 
!-------------------------------------------------------
else if (m==3) then    
n_atoms = n_lowercase*2*n3+2*n1+2*n2
call tight_binding1(n,n_atoms,n1,n2,n3) 
else
n_atoms = n_lowercase*2*n3+2*n1+2*n2+(m-3)*2*n1
call tight_binding1(n,n_atoms,n1,n2,n3) 
end if
end do

end 
    
