!###########################################################################################################
!!!Função que calcula o triangulo de Lumley(Estado de anisotropia da turbulencia)
subroutine Lumley(nx,ny,nz,nt,u,v,w,yplus,detadx,detady,dqsidx,dqsidy)

  implicit none
  integer(kind=4), intent(in)  :: nx,ny,nz,nt
  real(kind=8),    intent(in)  :: u(nx,ny,nz,nt),v(nx,ny,nz,nt),w(nx,ny,nz,nt),yplus(ny)
  real(kind=8),    intent(in)  :: detadx(nx,ny,nz),detady(nx,ny,nz),dqsidx(nx,ny,nz),dqsidy(nx,ny,nz)

  integer :: i,j,p,m,k,cont,lim
  real(kind=8)       :: qsi(ny),eta(ny),qsi1(101),eta1(101),eta2(101),eta3(101)
  real(kind=8)       :: ufluc(nx,ny,nz,nt),vfluc(nx,ny,nz,nt),wfluc(nx,ny,nz,nt)
  real(kind=8)       :: var1(nz,nt),var2(nz,nt),b(3,3)
  real(kind=8)       :: umean(nx,ny),vmean(nx,ny),wmean(nx,ny),traco,aux,invdois,invtres
  real(kind=8)       :: b11(ny),b12(ny),b13(ny),b21(ny),b22(ny),b23(ny),b31(ny),b32(ny),b33(ny)
  real(kind=8)       :: vtan(nx,ny,nz,nt),vnorm(nx,ny,nz,nt)
  
  i = (nx/2) +1                            

!################################################################################      
!########################Calculando as flutuacoes################################
!################################################################################
!As variaveis u, v e w possuem 4 dimensoes (x,y,z,tempo)

  open(unit=900, file='./Lumley/Lumley.dat', form='formatted')
  open(unit=700, file='./Lumley/Lumleyfit.dat', form='formatted')
  open(unit=400, file='./Lumley/ReynoldsNormalizadobij.dat', form='formatted')

!Primeiro ponto nao possui velocidade para o calculo do tensor de anisotropia
!Ultimo ponto no laço eh onde há o termino da camada limite, na variavel "lim"

!   do contador =  1, ny
!     if (yplus(contador) <= 500) then
!         lim = contador
!     endif
!   enddo

!Criando curvas que limitam o triangulo de Lumley
do cont=1,101
    qsi1(cont) = -1.0d0/6.0d0 + (cont-1)*0.005
    eta1(cont) = qsi1(cont)
    eta2(cont) = -qsi1(cont)
    eta3(cont) = ((1.0d0/27.0d0) + (2.0d0*qsi1(cont)**3))**(1.0d0/2.0d0)
    write(700,*) qsi1(cont),eta1(cont),eta2(cont),eta3(cont)
enddo

!Computa o trangulo de Lumley e o tensor de Reynolds com u_ij
if (.FALSE.) then
    
    print*, '  '
    print*, ' !!!!Computing Lumley triangle with u_ij!!!!' 
    print*, '  '
    
    call Double_Decomposition(u,nx,ny,nz,nt,ufluc,umean)  
    call Double_Decomposition(v,nx,ny,nz,nt,vfluc,vmean)
    call Double_Decomposition(w,nx,ny,nz,nt,wfluc,wmean)
  
    do j=2,lim
   
    traco = sum(sum(ufluc(i,j,:,:)*ufluc(i,j,:,:)+vfluc(i,j,:,:)*vfluc(i,j,:,:)+ & 
            wfluc(i,j,:,:)*wfluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt))
            
    do p=1,3  
        
        if (p == 1) then
            var1 = ufluc(i,j,:,:)
        endif
        if (p == 2) then
            var1 = vfluc(i,j,:,:)
        endif
        if (p == 3) then
            var1 = wfluc(i,j,:,:)
        endif
        do m=1,3    
        
            if (m == 1) then
                var2 = ufluc(i,j,:,:)
            endif
            if (m == 2) then
                var2 = vfluc(i,j,:,:)
            endif
            if (m == 3) then
                var2 = wfluc(i,j,:,:)
            endif
         
        aux = sum(sum(var1(:,:)*var2(:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) 

        !Calculando o tensor de Reynolds normalizado
            if (p /= m) then
                b(p,m) =   aux/traco
            endif
            if (p == m) then
                b(p,m) =   aux/traco - 1.0d0/3.0d0
            endif
        enddo
    enddo

    !Calculando os invariantes do tensor de Reynolds normalizado
 !Segundo invariante
    invdois=0.0d0
    do p=1,3
        do m=1,3
            invdois = invdois - (1.0d0/2.0d0)*(b(p,m)*b(m,p))
        enddo
    enddo

!Terceiro invariante
 invtres=0.0d0
    do p=1,3
        do m=1,3
            do k=1,3
                invtres = invtres + (b(p,m)*b(m,k)*b(k,p))/3.0d0
            enddo
        enddo
    enddo
   
    qsi(j) = ((invtres/2.0d0)/dabs(invtres/2.0d0))* dabs(invtres/2.0d0)**(1.0/3.0)
    eta(j) = (-invdois/3.0d0)**(1.0d0/2.0d0)
    b11(j) = b(1,1)
    b12(j) = b(1,2)
    b13(j) = b(1,3)
    b21(j) = b(2,1)
    b22(j) = b(2,2)
    b23(j) = b(2,3)
    b31(j) = b(3,1)
    b32(j) = b(3,2)
    b33(j) = b(3,3)

    write(400,*) yplus(j),b11(j),b12(j),b13(j),b21(j),b22(j),b23(j),b31(j),b32(j),b33(j)
    write(900,*) qsi(j),eta(j)

enddo
  close(400)
  close(900)

endif
!!!############################################################################################################
!!Computa o trangulo de Lumley e o tensor de Reynolds com velocidades normais e tangenciais

if (.TRUE.) then
    
    print*, '  '
    print*, ' !!!!Computing Lumley triangle with normal and tangencial components!!!!' 
    print*, '  '

    !Funcao no der.f90 que calcula as velocidades tangenciais e normais
    call variants(nx,ny,nz,nt,u,v,detadx,detady,dqsidx,dqsidy,vtan,vnorm)

    call Double_Decomposition(vtan,nx,ny,nz,nt,ufluc,umean)  
    call Double_Decomposition(vnorm,nx,ny,nz,nt,vfluc,vmean)
    call Double_Decomposition(w,nx,ny,nz,nt,wfluc,wmean)
 
    do j=2,ny
   
        traco = sum(sum(ufluc(i,j,:,:)*ufluc(i,j,:,:)+vfluc(i,j,:,:)*vfluc(i,j,:,:)+ & 
                wfluc(i,j,:,:)*wfluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt))
                
        do p=1,3  
            
            if (p == 1) then
                var1 = ufluc(i,j,:,:)
            endif
            if (p == 2) then
                var1 = vfluc(i,j,:,:)
            endif
            if (p == 3) then
                var1 = wfluc(i,j,:,:)
            endif
            do m=1,3    
            
                if (m == 1) then
                    var2 = ufluc(i,j,:,:)
                endif
                if (m == 2) then
                    var2 = vfluc(i,j,:,:)
                endif
                if (m == 3) then
                    var2 = wfluc(i,j,:,:)
                endif
             
            aux = sum(sum(var1(:,:)*var2(:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) 
    
            !Calculando o tensor de Reynolds normalizado
                if (p /= m) then
                    b(p,m) =   aux/traco
                endif
                if (p == m) then
                    b(p,m) =   aux/traco - 1.0d0/3.0d0
                endif
            enddo
        enddo

        !Calculando os invariantes do tensor de Reynolds normalizado
     !Segundo invariante
        invdois=0.0d0
        do p=1,3
            do m=1,3
                invdois = invdois - (1.0d0/2.0d0)*(b(p,m)*b(m,p))
            enddo
        enddo

    !Terceiro invariante
     invtres=0.0d0
        do p=1,3
            do m=1,3
                do k=1,3
                    invtres = invtres + (b(p,m)*b(m,k)*b(k,p))/3.0d0
                enddo
            enddo
        enddo

        qsi(j) = ((invtres/2.0d0)/dabs(invtres/2.0d0))* dabs(invtres/2.0d0)**(1.0/3.0)
        eta(j) = (-invdois/3.0d0)**(1.0d0/2.0d0)
        b11(j) = b(1,1)
        b12(j) = b(1,2)
        b13(j) = b(1,3)
        b21(j) = b(2,1)
        b22(j) = b(2,2)
        b23(j) = b(2,3)
        b31(j) = b(3,1)
        b32(j) = b(3,2)
        b33(j) = b(3,3)

        write(400,*) yplus(j),b11(j),b12(j),b13(j),b21(j),b22(j),b23(j),b31(j),b32(j),b33(j)
        write(900,*) qsi(j),eta(j)

    enddo
      close(400)
      close(900)
endif

end subroutine
