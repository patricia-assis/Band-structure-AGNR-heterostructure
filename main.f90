program main
use module_Data   
implicit none
integer n_atoms
integer N !  N values; N-AGNR-I(n,m)
integer n1,n2,n3
integer n_lowercase, m

do N = 3,5,2 !->N-AGNR-I(n,m)
do n_lowercase = 1, 1 !N-AGNR-I(->n,m)
do m =3,3!!N-AGNR-I(n,->m)
n3=(N+4)
n2=(N-2)+4
n1=(N-2)+2
!-------------------------------------------------------
if (m==2) then
n_atoms = n_lowercase*2*n3+2*n2
call tight_binding1(n,n_atoms,n1,n2,n3,n_lowercase,m)
!-------------------------------------------------------
else if (m==3) then    
n_atoms = n_lowercase*2*n3+2*n1+2*n2
call tight_binding1(n,n_atoms,n1,n2,n3,n_lowercase,m)
else
n_atoms = n_lowercase*2*n3+2*n1+2*n2+(m-3)*2*n1
call tight_binding1(n,n_atoms,n1,n2,n3,n_lowercase,m)
end if
end do
end do
end do
end 
    
    
   
