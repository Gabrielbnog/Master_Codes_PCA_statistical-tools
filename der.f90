!###########################################################################################################
!!!Função que retorna as derivadas em x,y e z
subroutine ddxyz(f,dfdx,dfdy,dfdz,nx,ny,nz)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), dimension(nx,ny,nz), intent(in) :: f
    real(8), dimension(nx,ny,nz), intent(out) :: dfdx,dfdy,dfdz
    integer ::i,j,k,der
    
    !Ordem de precisao das derivadas, podem ser escolhidos os esquemas explicitos de segunda ordem, compacto de 6th ordem ou 10th ordem 
    
    der = 10

    if (der == 2) then
        call ddx(f,dfdx,nx,ny,nz)
        call ddy(f,dfdy,nx,ny,nz)
        call ddz(f,dfdz,nx,ny,nz)
    endif
    
    if (der == 6 .or. der == 10) then
        
        do k = 1,nz
            do j = 1,ny
                if (der == 6) then
                    call compact_scheme_6th(nx,1.0d0,f(:,j,k),dfdx(:,j,k))
                endif        
                if (der == 10) then
                    call compact_scheme_10th(nx,1.0d0,f(:,j,k),dfdx(:,j,k))
                endif
            enddo
        enddo
    do k = 1,nz
        do i = 1,nx
            if (der == 6) then
                call compact_scheme_6th(ny,1.0d0,f(i,:,k),dfdy(i,:,k))
            endif    
            if (der == 10) then
                call compact_scheme_10th(ny,1.0d0,f(i,:,k),dfdy(i,:,k))
            endif  
        enddo
    enddo
    do j = 1,ny
        do i = 1,nx
            if (der == 6) then
                call compact_scheme_6th(nz,1.0d0,f(i,j,:),dfdz(i,j,:))
            endif
            if (der == 10) then
                call compact_scheme_10th(nz,1.0d0,f(i,j,:),dfdz(i,j,:))
            endif    
        enddo
    enddo
   endif
end subroutine
!###########################################################################################################
!###########################################################################################################
!Função que calcula os termos de metrica
subroutine metric(x,y,z,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dxdqsi,dydqsi,&
                  jacobian,nx,ny,nz)

implicit none
 integer, intent(in) :: nx,ny,nz
 real(8), dimension(nx,ny,nz), intent(in) :: x,y,z
 real(8), dimension(nx,ny,nz), intent(out) :: dqsidx,dqsidy,dqsidz
 real(8), dimension(nx,ny,nz), intent(out) :: detadx,detady,detadz
 real(8), dimension(nx,ny,nz), intent(out) :: dphidx,dphidy,dphidz,jacobian
 real(8), dimension(nx,ny,nz), intent(out) :: dxdqsi,dydqsi
 real(8), dimension(nx,ny,nz) :: dxdeta,dxdphi
 real(8), dimension(nx,ny,nz) :: dydeta,dydphi
 real(8), dimension(nx,ny,nz) :: dzdqsi,dzdeta,dzdphi

  call ddxyz(x,dxdqsi,dxdeta,dxdphi,nx,ny,nz)
  call ddxyz(y,dydqsi,dydeta,dydphi,nx,ny,nz)
  call ddxyz(z,dzdqsi,dzdeta,dzdphi,nx,ny,nz)
  
  !  Get the determinant of J
  jacobian = -dxdphi*dydeta*dzdqsi+dxdeta*dydphi*dzdqsi+dxdphi*dydqsi*dzdeta &
             -dxdqsi*dydphi*dzdeta-dxdeta*dydqsi*dzdphi+dxdqsi*dydeta*dzdphi

  !  Get inv(J)
  dqsidx = (-dydphi*dzdeta + dydeta*dzdphi)/jacobian
  dqsidy = ( dxdphi*dzdeta - dxdeta*dzdphi)/jacobian
  dqsidz = (-dxdphi*dydeta + dxdeta*dydphi)/jacobian

  detadx = ( dydphi*dzdqsi - dydqsi*dzdphi)/jacobian
  detady = (-dxdphi*dzdqsi + dxdqsi*dzdphi)/jacobian
  detadz = ( dxdphi*dydqsi - dxdqsi*dydphi)/jacobian
     
  dphidx = (-dydeta*dzdqsi + dydqsi*dzdeta)/jacobian
  dphidy = ( dxdeta*dzdqsi - dxdqsi*dzdeta)/jacobian
  dphidz = (-dxdeta*dydqsi + dxdqsi*dydeta)/jacobian
end subroutine
!###########################################################################################################
!###########################################################################################################
!Função que calcula gradiente juntamente com os termos de metrica
subroutine gradS(f,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), dimension(nx,ny,nz), intent(in) :: f,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz
    real(8), dimension(nx,ny,nz), intent(out) :: dfdx,dfdy,dfdz
    real(8), dimension(nx,ny,nz) :: dfdqsi,dfdeta,dfdphi
     call ddxyz(f,dfdqsi,dfdeta,dfdphi,nx,ny,nz)
     dfdx = dfdqsi*dqsidx + dfdeta*detadx + dfdphi*dphidx
     dfdy = dfdqsi*dqsidy + dfdeta*detady + dfdphi*dphidy
     dfdz = dfdqsi*dqsidz + dfdeta*detadz + dfdphi*dphidz
end subroutine
!###########################################################################################################
!###########################################################################################################
!Função que calcula derivada em x segunda ordem
subroutine ddx(f,dfdx,nx,ny,nz)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), dimension(nx,ny,nz), intent(in) :: f
    real(8), dimension(nx,ny,nz), intent(out) :: dfdx
    integer :: i
     if(nx==1)then
      dfdx = 0.0d0
     else
      do i=2,nx-1
        dfdx(i,:,:) = (f(i+1,:,:) - f(i-1,:,:))/2.0d0
      enddo
      dfdx(1,:,:)  = f(2,:,:) -f(1,:,:)
      dfdx(nx,:,:) = f(nx,:,:) -f(nx-1,:,:)
    endif
     return
end subroutine
!###########################################################################################################
!###########################################################################################################
!Função que calcula derivada em y segunda ordem
subroutine ddy(f,dfdy,nx,ny,nz)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), dimension(nx,ny,nz), intent(in)  :: f
    real(8), dimension(nx,ny,nz), intent(out) :: dfdy
    integer :: i
     if(ny==1)then
      dfdy = 0.0d0
     else
      do i=2,ny-1
        dfdy(:,i,:) = (f(:,i+1,:) - f(:,i-1,:))/2.0d0
      enddo
      dfdy(:,1,:)  = f(:,2,:) - f(:,1,:)
      dfdy(:,ny,:) = f(:,ny,:) -f(:,ny-1,:)
    endif
end subroutine
!###########################################################################################################
!###########################################################################################################
!Função que calcula derivada em z segunda ordem
subroutine ddz(f,dfdz,nx,ny,nz)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(8), dimension(nx,ny,nz), intent(in)  :: f
    real(8), dimension(nx,ny,nz), intent(out) :: dfdz
    integer :: i
     if(nz==1)then
      dfdz = 0.0d0
     else
      do i=2,nz-1
        dfdz(:,:,i) = (f(:,:,i+1) - f(:,:,i-1))/2.0d0
      enddo
      dfdz(:,:,1) = f(:,:,2) -f(:,:,1) 
      dfdz(:,:,nz)= f(:,:,nz) -f(:,:,nz-1)
    endif
end subroutine
!###########################################################################################################
!###########################################################################################################
!Função que calcula derivada sexta ordem compacta
subroutine compact_scheme_6th(n,h,vec,diff)

	! This subroutine finds the derivative of a function (vec) using a 6th order compact scheme.
	! It is being considered a non-periodic boundary condition.

implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: vec(n), diff(n)
real(kind=8), intent(in) :: h 
integer :: j
real(kind=8), dimension(n) :: diagA, diagB, diagC
real(kind=8), dimension(n) :: rhs
real(kind=8) :: a, b, c
real(kind=8) :: alpha, beta, alpha1, alpha2
! Define the matrix coefficients
! j = 1 and n
alpha1 = 5.0d0
! j=2 and n-1
alpha2 = 2.0d0/11.0d0
! else
alpha = 1.0d0/3.0d0
beta  = 0.0d0
a = 14.0d0/18.0d0
b = 1.0d0/36.0d0
c = 0.0d0
! Builting tridiagonal vectors (diagonal vectors)
     call diagvec_6th(diagA,diagB,diagC,n,alpha,alpha1,alpha2)
 ! Compute RHS
! j = 1
    rhs(1) = ((-197.0d0/60.0d0)*vec(1) + (-5.0d0/12.0d0)*vec(2) + 5.0d0*vec(3) + (-5.0d0/3.0d0)*vec(4) + (5.0d0/12.0d0)*vec(5) + (-1.0d0/20.0d0)*vec(6))/h
    j = 2
    rhs(j) = ((-20.0d0/33.0d0)*vec(j-1) + (-35.0d0/132.0d0)*vec(j) + (34.0d0/33.0d0)*vec(j+1) + (-7.0d0/33.0d0)*vec(j+2) + (2.0d0/33.0d0)*vec(j+3) + (-1.0d0/132.0d0)*vec(j+4))/h

! j = n
    rhs(n) = ((197.0d0/60.0d0)*vec(n) + (5.0d0/12.0d0)*vec(n-1) + (-5.0d0)*vec(n-2) + (5.0d0/3.0d0)*vec(n-3) + (-5.0d0/12.0d0)*vec(n-4) + (1.0d0/20.0d0)*vec(n-5))/h
    j=n-1
    rhs(j) = ((20.0d0/33.0d0)*vec(j+1) + (35.0d0/132.0d0)*vec(j) + (-34.0d0/33.0d0)*vec(j-1) + (7.0d0/33.0d0)*vec(j-2) + (-2.0d0/33.0d0)*vec(j-3) + (1.0d0/132.0d0)*vec(j-4))/h  
     do j=3,n-2
    	rhs(j) = (vec(j+1) - vec(j-1))*(a/h) + (vec(j+2) - vec(j-2))*(b/h) !+ (vec(j+3) - vec(j-3))*(c/h)
     enddo
! Solve the linear system to find the derivative
     call tridiagonal(n,diagA,diagB,diagC,rhs,diff) 
return

end subroutine compact_scheme_6th
!###########################################################################################################
!###########################################################################################################
subroutine diagvec_6th(A,B,C,n,alpha,alpha1,alpha2)

implicit none
integer i ,n
real(8) A(n), B(n), C(n)
real(8) alpha, alpha1, alpha2

  B = 1.0d0

    do i=1,n
      A(i) = alpha
      C(i) = alpha
    enddo

    C(1) = alpha1
    C(2) = alpha2

    A(1)  = 0.0d0
    A(2)  = alpha2
    
    A(n) = alpha1
    A(n-1) = alpha2

    C(n)   = 0.0d0
    C(n-1) = alpha2
    return
end subroutine diagvec_6th
!###########################################################################################################
!###########################################################################################################
!########################## Solve the tridiagonal system ##########################
subroutine tridiagonal(n,A1,B1,C1,D1,x)

	! This subroutine solve tridiagonal matrix    B[A,D,C] {x} = {b}

	! Ward Cheney and David Kincaid, "Numerical Mathematics and Computing",
	! Thomson Brooks/Cole, 6th ed., chap. 7 - pg 284 (2008).

	! |B(1) C(1)                                        |
	! |A(1) B(2) C(2)                                   |
	! |     A(2) B(3) C(3)                              |
	! |          A(3) B(4) C(4)                         |
	! |                                                 |
	! |                      A(n-3) B(n-2) C(n-2)       |
	! |                             A(n-2) B(n-1) C(n-1)| 
	! |                                    A(n-1) B(N  )|

Implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: A1(n),B1(n),C1(n),D1(n)
integer :: i
real(kind=8) :: A(n),B(n),C(n),D(n),X(n)
		   
  ! Coeficiente c 
  C(1) = C1(1)/B1(1) 
  do i = 2,n-1
    C(i) = C1(i)/(B1(i) - A1(i)*C(i-1))
  end do
      
  ! Coeficiente d 
  D(1) = D1(1)/B1(1)
  do i = 2,n
    D(i) = (D1(i) - A1(i)*D(i-1))/(B1(i) - A1(i)*C(i-1))
  end do

  ! Substituicao regressiva
  x(n) = D(n)
  do i = n-1,1,-1
    x(i) = D(i) - C(i)*x(i+1)
  end do

return

end subroutine tridiagonal
!###########################################################################################################
!###########################################################################################################
subroutine compact_scheme_10th(n,h,vec,diff)

	! This subroutine finds the derivative of a function (vec) using a 10th order compact scheme.
	! It is being considered a non-periodic boundary condition.

implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: vec(n), diff(n)
real(kind=8), intent(in) :: h
integer :: j
real(kind=8), dimension(n) :: diagE, diagA, diagD, diagC, diagF
real(kind=8), dimension(n) :: rhs
real(kind=8) :: a, b, c, a2, a3, b3, a4, b4, c4
real(kind=8) :: alpha, beta, alpha1, alpha2, alpha3, beta3, alpha4, beta4


	! Define the matrix coefficients

! j = 1 and n
alpha1 = 2.0d0
! j=2 and n-1
alpha2 = 0.25d0
a2 = 3.0d0/4.0d0 
! j=3 and n-2
alpha3 = 4.7435d0/10.67175d0
beta3  = 0.2964375d0/10.67175d0
a3 = 7.905d0/10.67175d0
b3 = 1.23515625d0/10.67175d0
! j=4 and n-3
alpha4 = 4.63271875d0/9.38146875d0
beta4  = 0.451390625d0/9.38146875d0
a4 = 6.66984375d0/9.38146875d0
b4 = 1.53/9.38146875d0
c4 = 0.015/9.38146875d0
! else
alpha = 0.5d0
beta  = 0.05d0
a = 17.0d0/24.0d0
b = 101.0d0/600.0d0
c = 0.01d0/6.0d0
! Builting pentadiagonal vectors (diagonal vectors)
    call diagvec(diagE,diagA,diagD,diagC,diagF,n,beta,beta3,beta4,alpha,alpha1,alpha2,alpha3,alpha4)
! Compute RHS
! j = 1
rhs(1) = (-2.5d0*vec(1) + 2.0d0*vec(2) + 0.5d0*vec(3))/h
j = 2
rhs(j) = (vec(j+1) - vec(j-1))*a2/h
j = 3
rhs(j) = (vec(j+1) - vec(j-1))*a3/h + (vec(j+2) - vec(j-2))*b3/h
j = 4
rhs(j) = (vec(j+1) - vec(j-1))*a4/h + (vec(j+2) - vec(j-2))*b4/h + (vec(j+3) - vec(j-3))*c4/h   
! j = n
rhs(n) = (2.5d0*vec(n) - 2.0d0*vec(n-1) - 0.5d0*vec(n-2))/h
j=n-1
rhs(j) = (vec(j+1) - vec(j-1))*a2/h  
j=n-2
rhs(j) = (vec(j+1) - vec(j-1))*a3/h + (vec(j+2) - vec(j-2))*b3/h
j=n-3
rhs(j) = (vec(j+1) - vec(j-1))*a4/h + (vec(j+2) - vec(j-2))*b4/h + (vec(j+3) - vec(j-3))*c4/h
    do j=5,n-4
    	rhs(j) = (vec(j+1) - vec(j-1))*a/h + (vec(j+2) - vec(j-2))*b/h + (vec(j+3) - vec(j-3))*c/h
    enddo
    ! Solve the linear system to find the derivative
    call penta(n,diagE,diagA,diagD,diagC,diagF,rhs,diff) 
    return

end subroutine compact_scheme_10th
!###########################################################################################################
!###########################################################################################################
subroutine diagvec(E1,A2,D3,C4,F5,N,beta,beta3,beta4,alpha,alpha1,alpha2,alpha3,alpha4)

implicit none
integer i ,n
real(8) E1(n), A2(n), D3(n), C4(n), F5(n)
real(8) alpha, beta, alpha1, alpha2, alpha3, beta3, alpha4, beta4


d3 = 1.0d0
do i=1,n
	e1(i) = beta
	a2(i) = alpha
	c4(i) = alpha
	f5(i) = beta
enddo
c4(1) = alpha1
c4(2) = alpha2
c4(3) = alpha3
c4(4) = alpha4
f5(1) = 0.0d0
f5(2) = 0.0d0
f5(3) = beta3
f5(4) = beta4
a2(1) = alpha2
a2(2) = alpha3
a2(3) = alpha4
e1(1) = beta3
e1(2) = beta4
e1(n)   = 0.0d0
e1(n-1) = 0.0d0
e1(n-2) = 0.0d0
e1(n-3) = 0.0d0
e1(n-4) = beta3
e1(n-5) = beta4
a2(n)   = 0.0d0
a2(n-1) = alpha1
a2(n-2) = alpha2
a2(n-3) = alpha3
a2(n-4) = alpha4
c4(n)   = 0.0d0
c4(n-1) = alpha2
c4(n-2) = alpha3
c4(n-3) = alpha4
f5(n)   = 0.0d0
f5(n-1) = 0.0d0
f5(n-2) = beta3
f5(n-3) = beta4


return

end subroutine diagvec
!###########################################################################################################
!###########################################################################################################
!########################## Solve the pentadiagonal system ##########################
subroutine penta(n,E1,A1,D1,C1,F1,b,x)

	! This subroutine solve pentadiagonal matrix    B[E,A,D,C,F] {x} = {b}

	! Ward Cheney and David Kincaid, "Numerical Mathematics and Computing",
	! Thomson Brooks/Cole, 6th ed., chap. 7 - pg 284 (2008).

	! |D(1) C(1) F(1)                                   |
	! |A(1) D(2) C(2) F(2)                              |
	! |E(1) A(2) D(3) C(3) F(3)                         |
	! |     E(2) A(3) D(4) C(4) F(4)                    |
	! |                                                 |
	! |               E(n-4) A(n-3) D(n-2) C(n-2) F(n-2)|
	! |                      E(n-3) A(n-2) D(n-1) C(n-1)| 
	! |                             E(n-1) A(n-1) D(N  )|

Implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: E1(n),A1(n),D1(n),C1(n),F1(n)
integer :: i
real(kind=8) :: E(n),A(n),D(n),C(n),F(n),B(n),X(n), xmult
		    
E=E1; A=A1; D=D1; C=C1; F=F1
do i = 2,n-1
	xmult = A(i-1)/D(i-1)
	D(i) = D(i) - xmult*C(i-1)
	C(i) = C(i) - xmult*F(i-1)
	b(i) = b(i) - xmult*b(i-1)
	xmult = E(i-1)/D(i-1)
	A(i) = A(i) - xmult*C(i-1)
	D(i+1) = D(i+1) - xmult*F(i-1)
	B(i+1) = b(i+1) - xmult*b(i-1)
enddo
xmult = A(n-1)/D(n-1)
D(n) = D(n) - xmult*C(n-1)
x(n) = (b(n) - xmult*b(n-1))/D(n)
x(n-1) = (b(n-1) - C(n-1)*X(n))/D(n-1)
do i = n-2,1,-1     
	x(i) = (b(i) - F(i)*x(i+2) - C(i)*x(i+1))/D(i)
enddo
return

end subroutine penta
!###########################################################################################################
!###########################################################################################################

! !#########################################################################################
! Rotina que calcula as velocidade contravariantes e covariantes para a velocidade
! Esta rotina irá transformar a velocidade em suas componentes tangenciais e normais
! See NASA report Hung - Definition of contravariant velocity componenty 
! See too book's Pulliam for tangencial and normal velocity with metrics terms  

subroutine variants(nx,ny,nz,nt,u,v,detadx,detady,dqsidx,dqsidy,vtan,vnorm)

  implicit none
  integer(kind=4), intent(in)   :: nx,ny,nz,nt                                                                          
  real(kind=8),    intent(in)   :: u(nx,ny,nz,nt),v(nx,ny,nz,nt)                                                                                                           
  real(kind=8),    intent(in)   :: detadx(nx,ny,nz),detady(nx,ny,nz),dqsidx(nx,ny,nz),dqsidy(nx,ny,nz)
  real(kind=8),    intent(out)  :: vtan(nx,ny,nz,nt),vnorm(nx,ny,nz,nt)
  real(kind=8)                  :: normx(nx),normy(nx),normtx(nx),normty(nx)
  integer(kind=4)               :: i,t,cont

  i = (nx/2) +1

!Componentes do vetor normal e tangencial(Utilizar somente para j=1 --> acompanha a normal do aerofolio)

!Componentes do vetor normal
! normx(:) = detadx(:,1,1)/((dqsidx(:,1,1)**2+dqsidy(:,1,1)**2))**(1.0d0/2.0d0)
  normx(:) = detadx(:,1,1)/((detadx(:,1,1)**2+detady(:,1,1)**2))**(1.0d0/2.0d0)
  normy(:) = detady(:,1,1)/((detadx(:,1,1)**2+detady(:,1,1)**2))**(1.0d0/2.0d0)

!Componentes do vetor tangencial(Pode-se utilizar uma matriz de rotação)
normtx(:) = normy(:)
normty(:) = -normx(:)

do t=1,nt
  do cont= 1,nx
    ! vnorm(cont,:,:,t) = (normx(cont)*u(cont,:,:,t)+normy(cont)*v(cont,:,:,t))/((normx(cont)**2+normy(cont)**2))**(1.0d0/2.0d0)
    vnorm(cont,:,:,t) = (normx(cont)*u(cont,:,:,t)+normy(cont)*v(cont,:,:,t))
    vtan(cont,:,:,t) =  (normtx(cont)*u(cont,:,:,t)+normty(cont)*v(cont,:,:,t))
  enddo
enddo

endsubroutine variants