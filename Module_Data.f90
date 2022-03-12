!******************************************************************!
!                              MODULE                              !
!******************************************************************!
!------------------------------------------------------------------!
!========================Data input================================!
!------------------------------------------------------------------!
module module_Data     
implicit none       
integer i,j,ii,jj,iii,kk,ngridy !counters 
!integer m ! N-AGNR-I(n_lowercase,m)
!parameter(m=3,n_lowercase=1)    
parameter(ngridy=1000)
double precision  pi
parameter(pi=dacos(-1.0d0))
complex*16 zi
double precision t
parameter(zi=(0.0d0,1.0d0))
end module
!*************************************************************************
