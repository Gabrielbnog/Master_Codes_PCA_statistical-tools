!####################################################################################
subroutine probe(i,j,k,m,nsnapshots,dt,path)

  use mod_field
  implicit none
  integer(kind=4), intent(in) :: i
  integer(kind=4), intent(in) :: j
  integer(kind=4), intent(in) :: k
  integer(kind=4), intent(in) :: m
  integer(kind=4), intent(in) :: nsnapshots
  real(kind=8), intent(in)    :: dt
  character(len=*), intent(in) :: path

  integer(kind=4) :: t

  character(len=4) char1, char2, char3, char4, char5

  if (allocated(zone(m)%q) .eqv. .false.) stop ' Tem cagada na alocacao de "zone(m)%q" ...'

  write(char1,'(i4.4)') m
  write(char2,'(i4.4)') i
  write(char3,'(i4.4)') j
  write(char4,'(i4.4)') k
  write(char5,'(i4.4)') nsnapshots

  write(*,*) 'Probe at :: zone = '//char1//', i = '//char2//', j = '//char3//', k = '//char4

  open(50,file = trim(path)//'/signal_m'//char1//'_i'//char2//'_j'//char3//'_k'//char4//'_n'//char5//'.dat')
  do t = 1,nsnapshots
    write(50,*) float(t)*dt, zone(m)%q(i,j,k,t)
  enddo
  close(50)

  return

end subroutine
!####################################################################################

!####################################################################################
subroutine circular_probe(j,k,m,nsnapshots,dt,path)

  use mod_field
  implicit none
  integer(kind=4), intent(in) :: j
  integer(kind=4), intent(in) :: k
  integer(kind=4), intent(in) :: m
  integer(kind=4), intent(in) :: nsnapshots
  real(kind=8), intent(in)    :: dt
  character(len=*), intent(in) :: path

  integer(kind=4) :: i, t

  character(len=4) char1, char2, char3, char4, char5
  
  real(kind=8) :: aux, pi
  
  pi = dacos(-1.0d0)

  write(char1,'(i4.4)') m
  write(char2,'(i4.4)') zone(m)%nx-1
  write(char3,'(i4.4)') j
  write(char4,'(i4.4)') k
  write(char5,'(i4.4)') nsnapshots

  write(*,*) 'Circular probe at :: zone = '//char1//', i = '//char2//', j = '//char3//', k = '//char4
  
  open(40,file = trim(path)//'/time.dat')
  do t = 1,nsnap
    write(40,'(f18.12)') float(t)*dt
  enddo
  close(40)
  
  open(45,file = trim(path)//'/theta.dat')
  do i = 1,zone(m)%nx-1
    aux = datan2(zone(m)%y(i,j,k),zone(m)%x(i,j,k))*180.0d0/pi
    if (aux .lt. 0.0d0) then
      aux = aux + 360.0d0
    endif
    write(45,'(f18.12)') aux
  enddo
  close(45)
  
  open(50,file = trim(path)//'/circular_m'//char1//'_i'//char2//'_j'//char3//'_k'//char4//'_n'//char5//'.dat')
  do i = 1,zone(m)%nx-1
    write(50,'(9999f18.12)') zone(m)%q(i,j,k,:)
  enddo
  close(50)

  return

end subroutine
!####################################################################################

!####################################################################################
subroutine probe_along_Z(nx,ny,nz,nt,u,v,w,detadx,detady,dqsidx,dqsidy,position,yplus)

  implicit none
  real(kind=8),    intent(in)  :: u(nx,ny,nz,nt),v(nx,ny,nz,nt),w(nx,ny,nz,nt),yplus(ny)
  integer(kind=4), intent(in)  :: nx,ny,nz,nt,position
  real(kind=8),    intent(in)  :: detadx(nx,ny,nz),detady(nx,ny,nz),dqsidx(nx,ny,nz),dqsidy(nx,ny,nz)
  integer(kind=4)    ::  i,t,temp
  real(kind=8)       :: vtan(nx,ny,nz,nt),vnorm(nx,ny,nz,nt)
  character(len=2) char1
  
  print*, 'VALOR DE POSITION'  
  print*, position
  i = (nx/2) +1  
  print*, 'Probe em yplus ::', yplus(position)
  temp = yplus(position)
  print*, 'VALOR DE temp'  
  print*, temp
  
  write(char1,'(i2.2)') temp
  call variants(nx,ny,nz,nt,u,v,detadx,detady,dqsidx,dqsidy,vtan,vnorm)

  open(unit=100, file='./PDF/Probe_u_tangencial090_'//char1// '.dat',form='formatted')
  open(unit=200, file='./PDF/Probe_v_normal090_'//char1//'.dat',form='formatted')
  open(unit=300, file='./PDF/Probe_w090_'//char1//'.dat',form='formatted')
   
  do t = 1,nt
    write(100,'(9999f18.12)') , vtan(i,position,:,t)
    write(200,'(9999f18.12)') , vnorm(i,position,:,t)
    write(300,'(9999f18.12)') , w(i,position,:,t)
  
  
  enddo
  
  close(100)
  close(200)
  close(300)

  return

end subroutine
!####################################################################################
