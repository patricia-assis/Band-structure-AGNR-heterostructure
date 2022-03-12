subroutine tight_binding1(n,n_atoms,n1,n2,n3,n_lowercase,m)
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
double precision g2,g3,g4
double precision at

double precision  acc2,acc3,acc4
integer n_lowercase, m ! N-AGNR-I(n_lowercase,m)

at = 1.0d0

g2=dfloat(m+1+n_lowercase*3)
g3=dfloat((n_lowercase+m-1)*3)
g4=dfloat((n_lowercase+m-1)*3)
acc2=at/g2!complex hopping amplitude phases for m = 2
acc3=at/g3!complex hopping amplitude phases for m = 3
acc4=at/g4!complex hopping amplitude phases for m >3

do i = 1, n_atoms
do j = 1, n_atoms
harm(i,j) = 0.0d0
hzh(i,j) = 0.0d0
end do
end do


lda=n_atoms
lwork=2*n_atoms-1
t = -3.0d0
         
do jj=1, n_atoms
w(jj) =0.0d0
end do
                     
do ii =-ngridy, ngridy
!do ii =-0.0d0, 0.0d0
k(ii) = (1.0d0)*((pi*dfloat(ii))/(ngridy))
end do

do ii =-ngridy, ngridy

if(m==2) then
v1 = t*exp(-zi*k(ii)*acc2)
v2 = t*exp(-zi*k(ii)*acc2/2)
v3 = t*exp(zi*k(ii)*acc2/2)

write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,  '.dat')") n, m
!write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,i2.2,  '.dat')") n, n_lowercase, m
open(unit = 10, file = file_name, status = "unknown", action = "write")
call Hamiltonian_m2(n_atoms,n1,n2,n3,harm,hzh,v1,v2,v3,n_lowercase,m)

else if (m==3) then 
v1 = t*exp(-zi*k(ii)*acc3)
v2 = t*exp(-zi*k(ii)*acc3/2)
v3 = t*exp(zi*k(ii)*acc3/2)

write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,i2.2,  '.dat')") n, n_lowercase, m
open(unit = 10, file = file_name, status = "unknown", action = "write")
call Hamiltonian_m3(n_atoms,n1,n2,n3,harm,hzh,v1,v2,v3,n_lowercase,m)

else
v1 = t*exp(-zi*k(ii)*acc4)
v2 = t*exp(-zi*k(ii)*acc4/2)
v3 = t*exp(zi*k(ii)*acc4/2)

write(unit = file_name, fmt="('.\Eigenvalues_',    '_',  i4.4,i2.2,i2.2,  '.dat')") n, n_lowercase, m
open(unit = 10, file = file_name, status = "unknown", action = "write")
call Hamiltonian_m4(n_atoms,n1,n2,n3,harm,hzh,v1,v2,v3,n,n_lowercase,m)
end if
       
call zheev('v','u',n_atoms,hzh,lda,w,work,lwork,rwork,info) 
!write(unit = 10, fmt = "(F18.12,1000(F18.12,x))") k(ii),(w(jj),jj=1,n_atoms)
write(unit = 10, fmt = "(F18.12,1000(F18.12,x))") dfloat(n_lowercase),(w(jj),jj=1,n_atoms)

if(info.gt.0)then
write(*,*)'The algorithm failed to compute eigenvalues.'
stop
end if
end do
!
end subroutine       
