!####################################################################################
!####################################################################################
subroutine velocity_profile(nx,ny,nz,nt,x,y,Re,Ma,u,v,r,utau,mu,dens,yplus)
!Funcao que calcula o perfil de pelocidade uplusyplus, dando como saida yplus.
!Para o calculo da derivada du/dn a função se utiliza de apenas velocidade V, sem a componente tangencial
  
  implicit none
  real(kind=8),               intent(out) :: yplus(ny),utau,mu,dens
  integer(kind=4),            intent(in)  :: nx,ny,nz,nt
  real(kind=8),               intent(in)  :: x(nx,ny,nz),y(nx,ny,nz)
  real(kind=8),               intent(in)  :: Re,Ma,u(nx,ny,nz,nt),v(nx,ny,nz,nt),r(nx,ny,nz,nt)
  integer(kind=4) :: j,i
  real(kind=8) :: tau_wall, dudn,Reinv,rho
  real(kind=8) :: u1, u2, u3 
  real(kind=8) :: dx, dy,  ds1, ds2, dn,cont
  real(kind=8), allocatable :: uplus(:)
  real(kind=8), allocatable :: u_aux(:,:), v_aux(:,:)   
  real(kind=8) :: dummy(nx,ny,nz,nt)
  
  allocate(uplus(1:ny),u_aux(1:nx,1:ny),v_aux(1:nx,1:ny)) 
  
  i = (nx/2) +1

  print*, ' '
  print*, ' Velocity profile at x =',  x(i,1,1)

  !calcular Cf usando du/dy
  !usar formula de segunda ordem (página 196 do Hirsch)
  dx = x(i,2,1) - x(i,1,1)
  dy = y(i,2,1) - y(i,1,1)
  ds1 = dsqrt(dx**2 + dy**2)
  
  dx = x(i,3,1) - x(i,2,1)
  dy = y(i,3,1) - y(i,2,1)
  ds2 = dsqrt(dx**2 + dy**2)
  
  ! mean utau
   Reinv = 1.0d0/(Re/Ma)
   mu = Reinv

  u1 = sum(u(i,1,:,:))/dble(nt*nz)
  u2 = sum(u(i,2,:,:))/dble(nt*nz)
  u3 = sum(u(i,3,:,:))/dble(nt*nz)
  
  dudn = ((ds1+ds2)/(ds1*ds2))*(u2-u1) - ((ds1/ds2)/(ds1+ds2))*(u3-u1)
  tau_wall = mu*dudn  
  rho =  sum(sum(r(i,1,:,:),dim=2),dim=1)/(dble(nz)*dble(nt))
  dens = rho
  utau = dsqrt(tau_wall/rho)
  
  ! !utau para Reconstrução Re400k
  ! print*, ' '
  ! print*, 'FORÇANDO Utau!!!!!'
  ! print*, 'ROTINA PROFILE!'
  ! print*, ' '
  ! ! utau =  0.00649620504d0   !x/c = 0.40
  ! ! utau =  0.00579252276d0 !x/c = 0.50
  ! ! utau =  0.00541901091d0 !x/c = 0.60
  !  utau =  0.00520872615d0 !x/c = 0.70
  ! ! utau =  0.00495132175d0 !x/c = 0.80
  ! ! utau =  0.00455792437d0 !x/c = 0.90

  ! Computing the fluctuation data
  !Na variavel dummy é guardado a flutuação, mas no caso nao se utiliza
  !nas variaveis locais auxiliares sao guardadas os valores das velocidades medias
  
  call Double_Decomposition(u,nx,ny,nz,nt,dummy,u_aux)  
  call Double_Decomposition(v,nx,ny,nz,nt,dummy,v_aux)
  
  open(unit=203, file='./Mean_profile/uplusyplus.dat', form='formatted')
  open(unit=405, file='./yyplus.dat', form='formatted')
  cont=0
  do j = 1,ny
    
    u1 = u_aux(i,j) 
    uplus(j) = u1/utau
   
    dx = x(i,j,1) - x(i,1,1)
    dy = y(i,j,1) - y(i,1,1)
    dn = dsqrt(dx**2 + dy**2)
    
    
    rho =  sum(sum(r(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt))      
    yplus(j) = dn*utau/(mu/rho)
    cont = cont+1 

    write(203,*) yplus(j), uplus(j)  
    write(405,*) cont,y(i,j,1), yplus(j)
  enddo
  close(203)
  close(405)

  return

end subroutine
!####################################################################################
!####################################################################################
subroutine ReynoldsStress_profile(nx,ny,nz,nt,x,u,v,w,yplus,utau)
  
  implicit none
  integer(kind=4),            intent(in)  ::  nx,ny,nz,nt
  real(kind=8),               intent(in)  ::  x(nx,ny,nz),u(nx,ny,nz,nt),v(nx,ny,nz,nt)
  real(kind=8),               intent(in)  ::  w(nx,ny,nz,nt),yplus(ny),utau

  integer(kind=4) :: j, k,t,i
  real(kind=8) :: uu, uv, uw, vv, vw, ww  
  real(kind=8) :: dummy1(nx,ny),dummy2(nx,ny),dummy3(nx,ny)
  real(kind=8) :: ufluc(nx,ny,nz,nt),vfluc(nx,ny,nz,nt),wfluc(nx,ny,nz,nt)
  
  i=(nx/2)+1

  print*, ' '
  print*, ' Reynolds stress at x =',  x(i,1,1)
  print*, ' '
    
  ! Computing the fluctuation data 
  call Double_Decomposition(u,nx,ny,nz,nt,ufluc,dummy1) 
  call Double_Decomposition(v,nx,ny,nz,nt,vfluc,dummy2)
  call Double_Decomposition(w,nx,ny,nz,nt,wfluc,dummy3)
 
  open(unit=300, file='./Reynolds/ReynoldsStress.dat', form='formatted')
  ! Computing the Reynolds stress tensor
  do j = 1,ny  
      uu = 0.0d0
      uv = 0.0d0
      uw = 0.0d0
      vv = 0.0d0
      vw = 0.0d0
      ww = 0.0d0
        do k = 1,nz  
          do t = 1,nt
            uu = uu + ufluc(i,j,k,t)*ufluc(i,j,k,t)/utau**2
            uv = uv + ufluc(i,j,k,t)*vfluc(i,j,k,t)/utau**2
            uw = uw + ufluc(i,j,k,t)*wfluc(i,j,k,t)/utau**2
            vv = vv + vfluc(i,j,k,t)*vfluc(i,j,k,t)/utau**2
            vw = vw + vfluc(i,j,k,t)*wfluc(i,j,k,t)/utau**2
            ww = ww + wfluc(i,j,k,t)*wfluc(i,j,k,t)/utau**2
          enddo
        enddo
        uu = uu/dble(nz*nt)
        uv = uv/dble(nz*nt)
        uw = uw/dble(nz*nt)
        vv = vv/dble(nz*nt)
        vw = vw/dble(nz*nt)
        ww = ww/dble(nz*nt)
    
    write(300,'(7f24.12)') yplus(j),uu, uv, uw, vv, vw, ww
    
  enddo  
  
  close(300)
  return
end subroutine

!##########################################################################################################
