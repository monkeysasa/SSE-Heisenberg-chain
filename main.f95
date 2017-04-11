!--------------------!
 module configuration
!--------------------!
 save

 integer :: lx       ! system length
 integer :: nh       ! number of H-operators in string
 integer :: mm       ! maximum string length

 real(8) :: beta     ! inverse temperature
 real(8) :: aprob    ! part of the acceptance probability for adding operator
 real(8) :: dprob    ! part of the acceptance probability for removing operator

 integer, allocatable :: spin(:)      ! spin state
 integer, allocatable :: bsites(:,:)  ! list of sites bsites(1,b),bsites(2,b) at bond b
 integer, allocatable :: opstring(:)  ! operator string

 integer, allocatable :: frstspinop(:) ! first operation on each site in linked vertex list
 integer, allocatable :: lastspinop(:) ! last operation on each site in linked vertex list
 integer, allocatable :: vertexlist(:) ! list of vertex links

 end module configuration
!------------------------!

!----------------------!
 module measurementdata
!----------------------!
 save

 real(8) :: enrg1=0.d0
 real(8) :: ususc=0.d0
 real(8) :: data1(2)=0.d0
 real(8) :: data2(2)=0.d0

!--------------------------!
 end module measurementdata
!--------------------------!

!============================!
 program basic_heisenberg_sse
!============================!
 use configuration; implicit none

 integer :: i,j,nbins,msteps,isteps

 open(10,file='read.in',status='old')
 read(10,*)lx,beta
 read(10,*)nbins,msteps,isteps
 close(10)

 call initran(1)
 call makelattice()
 call initconfig()                      !Initialize a random configuration

 aprob=0.5d0*beta*lx
 dprob=1.d0/(0.5d0*beta*lx)

 do i=1,isteps
    call diagonalupdate()
    call linkvertices()
    call loopupdate()
    call adjustcutoff(i)
 enddo

 do j=1,nbins
    do i=1,msteps
       call diagonalupdate()
       call linkvertices()
       call loopupdate()
       call measure()
    enddo
    call writeresults(msteps,j)
 enddo

 call deallocateall()

 end program basic_heisenberg_sse
!================================!

!---------------------------!
 subroutine diagonalupdate()
!---------------------------!
 use configuration; implicit none

 integer :: i,b,op
 real(8) :: ran

 do i=0,mm-1
    op=opstring(i)
    if (op==0) then
       b=int(ran()*lx)+1
       if (spin(bsites(1,b))/=spin(bsites(2,b))) then
          if (ran()*(mm-nh)<=aprob) then
             opstring(i)=2*b
             nh=nh+1
          endif
       endif
    elseif (mod(op,2)==0) then
       if (ran()<=dprob*(mm-nh+1)) then
          opstring(i)=0
          nh=nh-1
       endif
    else
       b=op/2
       spin(bsites(1,b))=-spin(bsites(1,b))
       spin(bsites(2,b))=-spin(bsites(2,b))
    endif
 enddo

 end subroutine diagonalupdate
!-----------------------------!

!-------------------------!
 subroutine linkvertices()
!-------------------------!
 use configuration; implicit none

 integer :: b,op,s1,s2,v0,v1,v2

 frstspinop(:)=-1
 lastspinop(:)=-1

 do v0=0,4*mm-1,4
    op=opstring(v0/4)
    if (op/=0) then
       b=op/2
       s1=bsites(1,b)
       s2=bsites(2,b)
       v1=lastspinop(s1)
       v2=lastspinop(s2)
       if (v1/=-1) then
          vertexlist(v1)=v0
          vertexlist(v0)=v1
       else
          frstspinop(s1)=v0
       endif
       if (v2/=-1) then
          vertexlist(v2)=v0+1
          vertexlist(v0+1)=v2
       else
          frstspinop(s2)=v0+1
       endif
       lastspinop(s1)=v0+2
       lastspinop(s2)=v0+3
    else
       vertexlist(v0:v0+3)=-1
    endif
 enddo
 do s1=1,lx
    v1=frstspinop(s1)
    if (v1/=-1) then
        v2=lastspinop(s1)
        vertexlist(v2)=v1
        vertexlist(v1)=v2
    endif
 enddo

 end subroutine linkvertices
!---------------------------!

!-----------------------!
 subroutine loopupdate()
!-----------------------!
 use configuration; implicit none

 integer :: i,v0,v1,v2
 real(8) :: ran

 do v0=0,4*mm-1,2
    if (vertexlist(v0)<0) cycle
    v1=v0
    if (ran()<0.5d0) then
       do
          opstring(v1/4)=ieor(opstring(v1/4),1)
          vertexlist(v1)=-2
          v2=ieor(v1,1)
          v1=vertexlist(v2)
          vertexlist(v2)=-2
          if (v1==v0) exit
       enddo
    else
       do
          vertexlist(v1)=-1
          v2=ieor(v1,1)
          v1=vertexlist(v2)
          vertexlist(v2)=-1
          if (v1==v0) exit
       enddo
    endif
 enddo

 do i=1,lx
    if (frstspinop(i)/=-1) then
       if (vertexlist(frstspinop(i))==-2) spin(i)=-spin(i)
    else
       if (ran()<0.5) spin(i)=-spin(i)
    endif
 enddo

 end subroutine loopupdate
!-------------------------!

!--------------------!
 subroutine measure()
!--------------------!
 use configuration; use measurementdata; implicit none

 enrg1=enrg1+dfloat(nh)
 ususc=ususc+dfloat(sum(spin)/2)**2

 end subroutine measure
!----------------------!

!------------------------------------!
 subroutine writeresults(msteps,bins)
!------------------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,msteps,bins
 real(8) :: wdata1(2),wdata2(2)

 enrg1=enrg1/msteps
 ususc=ususc/msteps

 enrg1=-enrg1/(beta*lx)+0.25d0
 ususc=beta*ususc/lx

 data1(1)=data1(1)+enrg1
 data1(2)=data1(2)+ususc

 data2(1)=data2(1)+enrg1**2
 data2(2)=data2(2)+ususc**2

 wdata1(:)=data1(:)/bins
 wdata2(:)=data2(:)/bins
 wdata2(:)=sqrt(abs(wdata2(:)-wdata1(:)**2)/bins)

 open(10,file='results.txt',status='replace')
 write(10,*)' Cut-off L : ',mm
 write(10,*)' Number of bins completed : ',bins
 write(10,*)' ========================================='
 write(10,10)'  E/N       : ',wdata1(1),wdata2(1)
 write(10,10)'  magnet    : ',wdata1(2),wdata2(2)
 write(10,*)' ========================================='
 10 format(1x,a,2f14.8)
 close(10)

 enrg1=0.d0
 ususc=0.d0

 end subroutine writeresults
!---------------------------!

!-----------------------------!
 subroutine adjustcutoff(step)
!-----------------------------!
 use configuration; implicit none

 integer, allocatable :: stringcopy(:)
 integer :: mmnew,step

 mmnew=nh+nh/3
 if (mmnew<=mm) return

 allocate(stringcopy(0:mm-1))
 stringcopy(:)=opstring(:)
 deallocate(opstring)
 allocate(opstring(0:mmnew-1))
 opstring(0:mm-1)=stringcopy(:)
 opstring(mm:mmnew-1)=0
 deallocate(stringcopy)

 mm=mmnew
 deallocate (vertexlist)
 allocate(vertexlist(0:4*mm-1))

 open(unit=10,file='results.txt',status='replace')
 write(10,*)' Step: ',step,'  Cut-off L: ',mm
 close(10)

 end subroutine adjustcutoff
!---------------------------!

!-----------------------!
 subroutine initconfig()
!-----------------------!
 use configuration; implicit none

 integer :: i
 real(8) :: ran

 allocate(spin(lx))
 do i=1,lx
    spin(i)=2*int(2.*ran())-1
 enddo

 mm=20
 allocate(opstring(0:mm-1))
 opstring(:)=0
 nh=0

 allocate(frstspinop(lx))
 allocate(lastspinop(lx))
 allocate(vertexlist(0:4*mm-1))

 end subroutine initconfig
!-------------------------!

!------------------------!
 subroutine makelattice()
!------------------------!
 use configuration; implicit none

 integer :: s,s1

 allocate(bsites(2,lx))

 do s=1,lx
    bsites(1,s)=s
    if (s<lx) then
        bsites(2,s)=s+1
    else
        bsites(2,s)=1
    end if
 enddo

 end subroutine makelattice
!--------------------------!

!--------------------------!
 subroutine deallocateall()
!--------------------------!
 use configuration; implicit none

 deallocate (spin)
 deallocate (bsites)
 deallocate (opstring)
 deallocate (frstspinop)
 deallocate (lastspinop)
 deallocate (vertexlist)

 end subroutine deallocateall
!----------------------------!

!----------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!
