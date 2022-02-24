!******************************************************************!
!                              MODULE                              !
!******************************************************************!
!------------------------------------------------------------------!
!========================Data input================================!
!------------------------------------------------------------------!
module module_Data     
implicit none       
integer i,j,ii,jj,iii,kk,ngridy !counters 
integer n_lowercase ! N-AGNR-I(n_lowercase,m)
integer m ! N-AGNR-I(n_lowercase,m)
parameter(m=3,n_lowercase=1)    
parameter(ngridy=1000)
double precision  pi
parameter(pi=dacos(-1.0d0))
double precision g2,g3,g4
double precision at
parameter(at=1.0d0)
double precision  acc2,acc3,acc4
parameter(g2=dfloat(m+1+n_lowercase*3))
parameter(g3=dfloat((n_lowercase+m-1)*3))
parameter(g4=dfloat((n_lowercase+m-1)*3))
parameter(acc2=at/g2)!complex hopping amplitude phases for m = 2
parameter(acc3=at/g3)!complex hopping amplitude phases for m = 3
parameter(acc4=at/g4)!complex hopping amplitude phases for m >3
complex*16 zi
double precision t
parameter(zi=(0.0d0,1.0d0))
end module
!*************************************************************************
