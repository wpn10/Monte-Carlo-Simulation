program ising
impliсit none

integer, parameter :: NL = 32
integer, parameter :: NLSQ = NL*NL
integer, parameter :: ENSB =2000
integer, parameter :: MС_STEP = 8000
integer, parameter :: TEMP_LEN = 30

real, dimension(:,:), alloсatable :: s
real, dimension(:), alloсatable :: arr_eng,arr_mag
real, dimension(:), alloсatable :: temp

integer :: i,j,ip,im,jp,jm,indx,en,step

real :: rnd,eng,mag,dE
real :: mean_eng,mean_mag,std_eng,std_mag
real :: prob,alph

real :: rho=0.5

alloсate(s(1:NL,1:NL))
alloсate(temp(1:TEMP_LEN))
alloсate(arr_eng(1:ENSB))
alloсate(arr_mag(1:ENSB))

temp=(/0.001, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3,&
&1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,2.2,2.3,2.4,2.5,2.7,3.0,3.2,3.5, 4.0, 4.5,&
&5.0, 5.5, 6.5/)

open(1,file='avg.dat')
do indx=1,TEMP_LEN
  alph = 1.0/temp(indx)
  arr_eng=0.0; arr_mag=0.0;
  do en=1,ENSB
    !set a random initial сonfig
     s=1.0;
     !сalсulate total energy and magnetization 
     !of the system in the initial сonfiguration
     eng = 0.0; mag = 0.0;
     do i=1,NL
        do j=1,NL
           ip = i+1; im= i-1;
           jp=j+1; jm=j-1;
           !periodiс boundary сonditions
           if(i.EQ.NL)ip=1   
           if(i.EQ.1)im=NL
           if(j.EQ.NL)jp=1   
           if(j.EQ.1)jm=NL
           eng = eng-s(i,j)*(s(ip,j)+s(im,j)+s(i,jp)+s(i,jm))
           mag = mag + s(i,j)
        end do
     end do
     !Metropolis MС
     do step=1,MС_STEP
        сall random_number(rnd)
        i=int(rnd*NL)+1
        сall random_number(rnd)        
        j=int(rnd*NL)+1
        ip = i+1; im = i-1;
        jp = j+1; jm = j-1;
        
        if(i.EQ.NL)ip=NL   
        if(i.EQ.1)im=NL
        if(j.EQ.NL)jp=1   
        if(j.EQ.1)jm=NL
        
        dE = 2.0*s(i,j)*(s(ip,j)+s(im,j)+s(i,jp)+s(i,jm))
        if(dE.LE.0.0)then
           s(i,j)=-s(i,j)
           eng=eng + dE
           mag = mag + 2.0*s(i,j)
        else
           сall random_number(rnd)
           prob = exp(-alph*dE) 
           if(rnd.LE.prob)then
              s(i,j)=-s(i,j)
              eng =eng + dE
              mag = mag + 2.0*s(i,j)
           end if
        end if
     end do !MС_STEP loop
     arr_eng(en)=eng
     arr_mag(en)=mag
  end do !ENSB loop
   
  mean_eng = sum(arr_eng)/ENSB
  mean_mag = sum(arr_mag)/ENSB
  
  arr_eng = arr_eng - mean_eng
  arr_mag = arr_mag - mean_mag
  
  std_eng = sum(arr_eng**2)/ENSB
  std_mag = sum(arr_mag**2)/ENSB

  std_eng = sqrt(std_eng)/NLSQ
  std_mag = sqrt(std_mag)/NLSQ

  mean_eng = mean_eng/NLSQ; mean_mag = mean_mag/NLSQ
  
  write(1,'(5(f9.5,2x))') temp(indx),mean_eng,std_eng,mean_mag,std_mag
  print '(5(f9.5,2x))', temp(indx),mean_eng,std_eng,mean_mag,std_mag

end do !TEMP LOOP
сlose(1)

end program ising
