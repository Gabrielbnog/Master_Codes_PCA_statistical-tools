! !####################################################################################
! !Rotina que calcula todas as derivadas necessárias para o balanço de energia cinética turbulenta
! !para o caso do aerofólio que possui somente a direção z como homogenea de turbulencia
subroutine BalanceTKE(yplus,mu,utau,dens,nx,ny,nz,nt,u,v,w,p,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz)
                     
  implicit none
   integer(kind=4), intent(in)  :: nx,ny,nz,nt                                                                             
   real(kind=8),    intent(in)  :: mu,utau,dens
   real(kind=8),    intent(in)  :: u(nx,ny,nz,nt),v(nx,ny,nz,nt),w(nx,ny,nz,nt),p(nx,ny,nz,nt),yplus(ny)                                                                                                           
   real(kind=8),    intent(in)  :: dqsidx(nx,ny,nz),dqsidy(nx,ny,nz),dqsidz(nx,ny,nz)
   real(kind=8),    intent(in)  :: detadx(nx,ny,nz),detady(nx,ny,nz),detadz(nx,ny,nz)
   real(kind=8),    intent(in)  :: dphidx(nx,ny,nz),dphidy(nx,ny,nz),dphidz(nx,ny,nz)
   
   real(kind=8)       :: ufluc(nx,ny,nz,nt),vfluc(nx,ny,nz,nt),wfluc(nx,ny,nz,nt),pfluc(nx,ny,nz,nt)                        
   real(kind=8)       :: duuvdy(ny),dvvvdy(ny),dwwvdy(ny)
   real(kind=8)       :: duuudx(ny),dvvudx(ny),dwwudx(ny)
   real(kind=8)       :: dvpdy(ny),dupdx(ny)
   real(kind=8)       :: dkdy(ny),dkdx(ny)
   real(kind=8)       :: dUmdx(ny),dVmdx(ny),dUmdy(ny),dVmdy(ny)
   real(kind=8)       :: d2kdy2(ny),d2kdx2(ny)
   
   real(kind=8)       :: duuvdyt(nx,ny,nz,nt),dvvvdyt(nx,ny,nz,nt),dwwvdyt(nx,ny,nz,nt)
   real(kind=8)       :: duuudxt(nx,ny,nz,nt),dvvudxt(nx,ny,nz,nt),dwwudxt(nx,ny,nz,nt) 
   real(kind=8)       :: dvpdyt(nx,ny,nz,nt),dupdxt(nx,ny,nz,nt)
   real(kind=8)       :: dkdyt(nx,ny,nz,nt),dkdxt(nx,ny,nz,nt)
   real(kind=8)       :: dudxt(nx,ny,nz,nt),dvdxt(nx,ny,nz,nt),dudyt(nx,ny,nz,nt),dvdyt(nx,ny,nz,nt)
   real(kind=8)       :: duflucdyt(nx,ny,nz,nt),dwflucdyt(nx,ny,nz,nt),dvflucdxt(nx,ny,nz,nt),dwflucdxt(nx,ny,nz,nt)
   real(kind=8)       :: dUmdxt(nx,ny,nz,nt),dVmdxt(nx,ny,nz,nt),dUmdyt(nx,ny,nz,nt),dVmdyt(nx,ny,nz,nt)
   real(kind=8)       :: d2kdy2t(nx,ny,nz,nt),d2kdx2t(nx,ny,nz,nt),duflucdxt(nx,ny,nz,nt)
   real(kind=8)       :: dvflucdyt(nx,ny,nz,nt),duflucdzt(nx,ny,nz,nt),dvflucdzt(nx,ny,nz,nt),dwflucdzt(nx,ny,nz,nt)
   
   real(kind=8)       :: aux(nx,ny,nz),umean(nx,ny),vmean(nx,ny),wmean(nx,ny),pmean(nx,ny) 
   real(kind=8)       :: dfdx(nx,ny,nz),dfdy(nx,ny,nz),dfdz(nx,ny,nz)
   real(kind=8)       :: Transp(ny),Dp(ny),Pr(ny),C(ny),Dv(ny),D(ny),B(ny)
   real(kind=8)       :: Pr11(ny),Pr12(ny),Pr21(ny), Pr22(ny)
   integer(kind=4)    :: i,t,j

   i = (nx/2) +1                          
   print*, '  '
   print*, ' !!!!Computing the derivatives for the TKE!!!!' 
   print*, '  '
!################################################################################      
!########################Calculando as flutuacoes################################
!################################################################################
!As variaveis u, v e w e p possuem 4 dimensoes (x,y,z,tempo)
!Flutuacao de velocidade u e velocidade media umean 

  call Double_Decomposition(u,nx,ny,nz,nt,ufluc,umean)  
  call Double_Decomposition(v,nx,ny,nz,nt,vfluc,vmean)
  call Double_Decomposition(w,nx,ny,nz,nt,wfluc,wmean)
  call Double_Decomposition(p,nx,ny,nz,nt,pfluc,pmean)
      
  open(unit=800, file='./TKE/Balance/BalanceTKE.dat', form='formatted')
  open(unit=700, file='./TKE/Balance/Producao.dat', form='formatted')
 
  do t = 1,nt
!################################################################################  
!#####################      Transporte Turbulento            ####################
!################################################################################  
    !derivada - d(uuv)/dy 
    aux(:,:,:) = ufluc(:,:,:,t)*ufluc(:,:,:,t)*vfluc(:,:,:,t)
    call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    duuvdyt(:,:,:,t) = dfdy(:,:,:)
    
    !derivada - d(vvv)/dy 
    aux(:,:,:) = vfluc(:,:,:,t)*vfluc(:,:,:,t)*vfluc(:,:,:,t)
    call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    dvvvdyt(:,:,:,t) = dfdy(:,:,:)

    !derivada - d(wwv)/dy 
    aux(:,:,:) = wfluc(:,:,:,t)*wfluc(:,:,:,t)*vfluc(:,:,:,t)
    call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    dwwvdyt(:,:,:,t) = dfdy(:,:,:)
   
    !derivada - d(uuu)/dx 
    aux(:,:,:) = ufluc(:,:,:,t)*ufluc(:,:,:,t)*ufluc(:,:,:,t)
    call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    duuudxt(:,:,:,t) = dfdx(:,:,:)

    !derivada - d(vvu)/dx 
    aux(:,:,:) = vfluc(:,:,:,t)*vfluc(:,:,:,t)*ufluc(:,:,:,t)
    call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    dvvudxt(:,:,:,t) = dfdx(:,:,:)

    !derivada - d(wwu)/dx 
    aux(:,:,:) = wfluc(:,:,:,t)*wfluc(:,:,:,t)*ufluc(:,:,:,t)
    call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
    dwwudxt(:,:,:,t) = dfdx(:,:,:)
  
!################################################################################  
!#####################      Difusão de Pressão          #########################
!################################################################################

  !derivada - d(vp)/dy 
  aux(:,:,:) = vfluc(:,:,:,t)*pfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvpdyt(:,:,:,t) = dfdy(:,:,:)

  !derivada - d(up)/dx 
  aux(:,:,:) = ufluc(:,:,:,t)*pfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dupdxt(:,:,:,t) = dfdx(:,:,:)

!################################################################################  
!################      Energia Cinética Turbulenta K         ####################
!################  (1/2) do traço do tensor de Reynolds      ####################
!################              Convecção                     ####################
!################################################################################
!################################################################################
!################################################################################

!derivada - d(k)/dy 
  aux(:,:,:) = (ufluc(:,:,:,t)*ufluc(:,:,:,t)+vfluc(:,:,:,t)*vfluc(:,:,:,t)+wfluc(:,:,:,t)*wfluc(:,:,:,t))/2.0d0
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dkdyt(:,:,:,t) = dfdy(:,:,:)

!derivada - d(k)/dx 
  aux(:,:,:) = (ufluc(:,:,:,t)*ufluc(:,:,:,t)+vfluc(:,:,:,t)*vfluc(:,:,:,t)+wfluc(:,:,:,t)*wfluc(:,:,:,t))/2.0d0
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dkdxt(:,:,:,t) = dfdx(:,:,:)

!################################################################################  
!#####################         Produção       ###################################
!################################################################################  
  
  !derivada - dUm/dx  
  aux(:,:,:) = u(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dudxt(:,:,:,t) = dfdx(:,:,:)

  aux(:,:,:) = ufluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  duflucdxt(:,:,:,t) = dfdx(:,:,:)
  dUmdxt(:,:,:,t) = dudxt(:,:,:,t) - duflucdxt(:,:,:,t)
  
  !derivada - dVm/dx  
  aux(:,:,:) = v(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvdxt(:,:,:,t) = dfdx(:,:,:)

  aux(:,:,:) = vfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvflucdxt(:,:,:,t) = dfdx(:,:,:)
  dVmdxt(:,:,:,t) = dvdxt(:,:,:,t) - dvflucdxt(:,:,:,t)

  !derivada - dUm/dy 
  aux(:,:,:) = u(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dudyt(:,:,:,t) = dfdy(:,:,:)

  aux(:,:,:) = ufluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  duflucdyt(:,:,:,t) = dfdy(:,:,:)
  dUmdyt(:,:,:,t) = dudyt(:,:,:,t) - duflucdyt(:,:,:,t)
 
  !derivada - dVm/dy 
  aux(:,:,:) = v(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvdyt(:,:,:,t) = dfdy(:,:,:)

  aux(:,:,:) = vfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvflucdyt(:,:,:,t) = dfdy(:,:,:)
  dVmdyt(:,:,:,t) = dvdyt(:,:,:,t) - dvflucdyt(:,:,:,t)
 
!################################################################################  
!#####################    Difusão Viscosa            ############################
!################################################################################
  !derivada - d2k/dy2 
  aux(:,:,:) = dkdyt(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  d2kdy2t(:,:,:,t) = dfdy(:,:,:)
   
  !derivada - d2k/dx2 
  aux(:,:,:) = dkdxt(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  d2kdx2t(:,:,:,t) = dfdx(:,:,:)

! #################################################################################  
!#####################  Dissipação                #################################
!##################################################################################
  !derivada - du/dx
  aux(:,:,:) = ufluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  duflucdxt(:,:,:,t) = dfdx(:,:,:)
 
  !derivada - dv/dx
  aux(:,:,:) = vfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvflucdxt(:,:,:,t) = dfdx(:,:,:)

  !derivada - dw/dx
  aux(:,:,:) = wfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dwflucdxt(:,:,:,t) = dfdx(:,:,:)

  !derivada - du/dy
  aux(:,:,:) = ufluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  duflucdyt(:,:,:,t) = dfdy(:,:,:)
 
  !derivada - dv/dy
  aux(:,:,:) = vfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvflucdyt(:,:,:,t) = dfdy(:,:,:)

  !derivada - dw/dy
  aux(:,:,:) = wfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dwflucdyt(:,:,:,t) = dfdy(:,:,:)
  
  !derivada - du/dz
  aux(:,:,:) = ufluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  duflucdzt(:,:,:,t) = dfdz(:,:,:)
 
  !derivada - dv/dz
  aux(:,:,:) = vfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dvflucdzt(:,:,:,t) = dfdz(:,:,:)

  !derivada - dw/dz
  aux(:,:,:) = wfluc(:,:,:,t)
  call gradS(aux,dqsidx,dqsidy,dqsidz,detadx,detady,detadz,dphidx,dphidy,dphidz,dfdx,dfdy,dfdz,nx,ny,nz)
  dwflucdzt(:,:,:,t) = dfdz(:,:,:)
  
enddo   
!################################################################
!!Realizando a média das derivadas na envergadura e no tempo
!################################################################

  duuvdy(:) = sum(sum(duuvdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dvvvdy(:) = sum(sum(dvvvdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dwwvdy(:) = sum(sum(dwwvdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  duuudx(:) = sum(sum(duuudxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dvvudx(:) = sum(sum(dvvudxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dwwudx(:) = sum(sum(dwwudxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  
  dvpdy(:) = sum(sum(dvpdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dupdx(:) = sum(sum(dupdxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))

  dkdy(:) = sum(sum(dkdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dkdx(:) = sum(sum(dkdxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
 
  dUmdx(:) = sum(sum(dUmdxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dVmdx(:) = sum(sum(dVmdxt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dUmdy(:) = sum(sum(dUmdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  dVmdy(:) = sum(sum(dVmdyt(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  
  d2kdy2(:) = sum(sum(d2kdy2t(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  d2kdx2(:) = sum(sum(d2kdx2t(i,:,:,:),dim=3),dim=2)/(dble(nz)*dble(nt))
  
  do j = 1,ny
    !Termo de transporte turbulento
      Transp(j)  = - (1.0d0/2.0d0)*((duuvdy(j)+dvvvdy(j)+dwwvdy(j)+duuudx(j)+dvvudx(j)+dwwudx(j))*mu)/(utau**4)        
 
    !Termo de difusão de pressão
      Dp(j)      = - (1.0d0/dens)*((dvpdy(j) + dupdx(j))*mu)/(utau**4)  
 
    !Termo de produção  
      Pr(j)      = - (((sum(sum(ufluc(i,j,:,:)*ufluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dUmdx(j) + & 
                       (sum(sum(vfluc(i,j,:,:)*ufluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dVmdx(j) + &           
                       (sum(sum(ufluc(i,j,:,:)*vfluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dUmdy(j) + &
                       (sum(sum(vfluc(i,j,:,:)*vfluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dVmdy(j))*mu)/(utau**4)  
              
       Pr11(j) = - (((sum(sum(ufluc(i,j,:,:)*ufluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dUmdx(j))*mu)/(utau**4)  
       Pr12(j) = - (((sum(sum(vfluc(i,j,:,:)*ufluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dVmdx(j))*mu)/(utau**4)  
       Pr21(j) = - (((sum(sum(ufluc(i,j,:,:)*vfluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dUmdy(j))*mu)/(utau**4)  
       Pr22(j) = - (((sum(sum(vfluc(i,j,:,:)*vfluc(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))*dVmdy(j))*mu)/(utau**4)  
       
       Pr(j) = Pr11(j)+Pr12(j)+Pr21(j)+Pr22(j)

    !Termo de convecção
      C(j)       = - ((dkdy(j)*vmean(i,j)+dkdx(j)*umean(i,j))*mu)/(utau**4)  

    !Termo de difusão viscosa
      Dv(j)      =   ((mu**2)*(d2kdy2(j)+d2kdx2(j)))/(utau**4)
   
    !Dissipação
      D(j)   = - mu*(sum(sum(duflucdxt(i,j,:,:)*duflucdxt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(dvflucdxt(i,j,:,:)*dvflucdxt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(dwflucdxt(i,j,:,:)*dwflucdxt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(duflucdyt(i,j,:,:)*duflucdyt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(dvflucdyt(i,j,:,:)*dvflucdyt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(dwflucdyt(i,j,:,:)*dwflucdyt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(duflucdzt(i,j,:,:)*duflucdzt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(dvflucdzt(i,j,:,:)*dvflucdzt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)) + &
                     sum(sum(dwflucdzt(i,j,:,:)*dwflucdzt(i,j,:,:),dim=2),dim=1)/(dble(nz)*dble(nt)))/(utau**4/mu)
      
    !Balanco        
    B(j)= Transp(j)+Dp(j)+Pr(j)+C(j)+Dv(j)+D(j)

      write(800,*) yplus(j),Transp(j),Dp(j),Pr(j),C(j),Dv(j),D(j),B(j)
      write(700,*) yplus(j), Pr11(j),Pr12(j),Pr21(j),Pr22(j),Pr(j)
  enddo
  
     print*, 'Output in BalanceTKE curve of TKE'
     close(800)
     close(700)
     return

end subroutine
         
