subroutine tight_binding1(n,n_atoms,n1,n2,n3)
use module_Data
implicit none
integer n
integer n1,n2,n3
integer n_atoms
complex*16 v1,v2,v3
double precision k(-ngridy:ngridy)
complex*16 harm(n_atoms,n_atoms), hzh(n_atoms,n_atoms)
integer lda, lwork, info
complex*16 work(2*n_atoms-1)
double precision rwork(3*n_atoms-2), w(n_atoms)
character*256 file_name,file_name2,file_name3

lda=n_atoms
lwork=2*n_atoms-1
t = -3.0d0
         
do jj=1, n_atoms
w(jj) =0.0d0
end do
                     
do ii =-ngridy, ngridy
k(ii) = (1.0d0)*((pi*dfloat(ii))/(ngridy))
end do

do ii =-ngridy,ngridy !ngridy, ngridy

if(m==2) then
v1 = t*exp(-zi*k(ii)*acc2)
v2 = t*exp(-zi*k(ii)*acc2/2)
v3 = t*exp(zi*k(ii)*acc2/2)

write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,i2.2,  '.dat')") n, n_lowercase, m
open(unit = 10, file = file_name, status = "unknown", action = "write")
call Hamiltonian_m2(n_atoms,n1,n2,n3,harm,hzh,v1,v2,v3)

else if (m==3) then 
v1 = t*exp(-zi*k(ii)*acc3)
v2 = t*exp(-zi*k(ii)*acc3/2)
v3 = t*exp(zi*k(ii)*acc3/2)

write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,i2.2,  '.dat')") n, n_lowercase, m
open(unit = 10, file = file_name, status = "unknown", action = "write")
call Hamiltonian_m3(n_atoms,n1,n2,n3,harm,hzh,v1,v2,v3)

else
v1 = t*exp(-zi*k(ii)*acc4)
v2 = t*exp(-zi*k(ii)*acc4/2)
v3 = t*exp(zi*k(ii)*acc4/2)

write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,i2.2,  '.dat')") n, n_lowercase, m
open(unit = 10, file = file_name, status = "unknown", action = "write")
call Hamiltonian_m4(n_atoms,n1,n2,n3,harm,hzh,v1,v2,v3,n)
end if
       
call zheev('v','u',n_atoms,hzh,lda,w,work,lwork,rwork,info) 
write(unit = 10, fmt = "(F18.12,1000(F18.12,x))") k(ii),(w(jj),jj=1,n_atoms)

if(info.gt.0)then
write(*,*)'The algorithm failed to compute eigenvalues.'
stop
end if
end do
!
end subroutine       
