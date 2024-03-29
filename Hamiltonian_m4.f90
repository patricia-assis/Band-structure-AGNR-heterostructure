subroutine Hamiltonian_m4(n_atoms,b1,b2,b3,z,h,d11,d12,d13,n,n_lowercase,m)
use module_Data    
implicit none      
integer n_atoms     
integer b1,b2,b3     
complex*16 d11,d12,d13,d14     
complex*16 z(n_atoms,n_atoms), h(n_atoms,n_atoms)    
integer m1, m2, m3, m4, m5, m6, m7, m8
integer f1, f2, f3, f4, f5     
integer x     
integer n, v
character*256 file_name,file_name2,file_name3   
integer n_lowercase, m ! N-AGNR-I(n_lowercase,m)
!******************************************************************!
!******************************************************************!
do  i = 1,b1-1,2
j=i+1
z(i,j)  = d13
z(j,i)  = dconjg(d13)
end do        
!******************************************************************!
!******************************************************************!
do  i = 2,b1,2
j=i+1
z(i,j)    = d12
z(j,i)= dconjg(d12)
end do
!******************************************************************!
!******************************************************************! 
m1 = b1
do  i = 1,b2-1,2
j=i+1
z(i+m1,j+m1)   = d13
z(j+m1,i+m1)= dconjg(d13)
end do         
!******************************************************************!
!******************************************************************!
do  i = 2,b2,2
j=i+1
z(i+m1,j+m1) = d12
z(j+m1,i+m1) = dconjg(d12)
end do
!******************************************************************!
!******************************************************************!
m2 = b2 + b1
do x = 0,n_lowercase-1,1
do  i = 1,2*b3,4
j=i+1
z(i+x*(2*b3)+m2,j+x*(2*b3)+m2)  = d11
z(j+x*(2*b3)+m2,i+x*(2*b3)+m2)  = dconjg(d11)
end do         
!******************************************************************!
!******************************************************************!
do  i = 2,2*b3-1,2
j=i+1
z(i+x*(2*b3)+m2,j+x*(2*b3)+m2)  = d12
z(j+x*(2*b3)+m2,i+x*(2*b3)+m2)= dconjg(d12)
end do
!******************************************************************!
!******************************************************************!
do  i = 1,2*b3-2,2
j=i+3
z(i+x*(2*b3)+m2,j+x*(2*b3)+m2) = d13
z(j+x*(2*b3)+m2,i+x*(2*b3)+m2)= dconjg(d13)
end do
!******************************************************************!
!******************************************************************!
if (x.ne.0) then
do  i= 3,2*b3,4
j=i+1+2*b3
z(i+(x-1)*(2*b3)+m2,j+(x-1)*(2*b3)+m2) = d11
z(j+(x-1)*(2*b3)+m2,i+(x-1)*(2*b3)+m2) = dconjg( d11)
end do
end if
end do
!******************************************************************!
!******************************************************************!
m3 = b1+b2+b3*2*n_lowercase
do  i = 1,b2-1,2
j=i+1
z(i+m3,j+m3)  = d12
z(j+m3,i+m3)  = dconjg(d12)
end do         
!
do  i = 2,b2,2
j=i+1
z(i+m3,j+m3)    = d13
z(j+m3,i+m3)= dconjg(d13)
end do
!******************************************************************!
!******************************************************************!
m4=b1+b2+b3*2*n_lowercase
do  i = 1,b1-1,2
j=i+1
z(i+b2+m4,j+b2+m4)   = d12
z(j+b2+m4,i+b2+m4)= dconjg(d12)
end do         
!
do  i = 2,b1,2
j=i+1
z(i+b2+m4,j+b2+m4)   = d13
z(j+b2+m4,i+b2+m4)= dconjg(d13)
end do
!******************************************************************!
!******************************************************************!
m5 = 2*b2 + 2*b1+2*b3*n_lowercase
do x = 0,(m-3)-1,1
do  i = 1,2*b1,4
j=i+1
z(i+x*(2*b1)+m5,j+x*(2*b1)+m5)  = d11
z(j+x*(2*b1)+m5,i+x*(2*b1)+m5)  = dconjg(d11)
end do         
!
do  i = 2,2*b1-1,2
j=i+1
z(i+x*(2*b1)+m5,j+x*(2*b1)+m5)  = d12
z(j+x*(2*b1)+m5,i+x*(2*b1)+m5)= dconjg(d12)
end do
!
do  i = 1,2*b1-2,2
j=i+3
z(i+x*(2*b1)+m5,j+x*(2*b1)+m5) = d13
z(j+x*(2*b1)+m5,i+x*(2*b1)+m5)= dconjg(d13)
end do
!
if (x.ne.0) then
do  i= 3,2*b1,4
j=i+1+2*b1
z(i+(x-1)*(2*b1)+m5,j+(x-1)*(2*b1)+m5) = d11
z(j+(x-1)*(2*b1)+m5,i+(x-1)*(2*b1)+m5) = dconjg( d11)
end do
end if
end do
!******************************************************************!
!******************************************************************!
do  i = 1,b1,2
j=i+1
z(i,j+b1)   = d11
z(j+b1,i)= dconjg(d11)
end do
!
f1 = b1
do  i= 1,b2,2
j=i
z(i+f1,j+b2+(i+2)+f1) = d11
z(j+b2+(i+2)+f1,i+f1) = dconjg(d11)
end do
!******************************************************************!
!******************************************************************!
f2 =b1+b2+(2*n_lowercase-2)*b3-1
do  i= 3,b3,2
j=i+2*b3-1
z(2*(i-1)+f2,j+f2) = d11
z(j+f2,2*(i-1)+f2) = dconjg(d11)
end do
!******************************************************************!
!******************************************************************!
f3=b1+b2+b3*2*n_lowercase
do  i = 2,b2,2
j=i+b1+1
z(i+f3,j+f3)   = d11
z(j+f3,i+f3)= dconjg(d11)
end do
!
f4 = b1+2*b2+2*b3*n_lowercase
do  i= 2,b1-1,2
j=i+b2
z(i+f4,j+(i+2)+f4-4) = d11
z(j+(i+2)+f4-4,i+f4) = dconjg(d11)
end do
!
v = b1+2*b2+2*n_lowercase*b3+((m-3)*2-1)*b1+3
do  i= 1,b1-1,2
j=i+1
z(2*(i-1)+v,j) = d11
z(j,2*(i-1)+v) = dconjg(d11)
end do
!******************************************************************!
!******************************************************************!
do i=1,n_atoms
do j=1,n_atoms
h(i,j)=z(i,j)
end do
end do
!******************************************************************!
!******************************************************************!
end subroutine
