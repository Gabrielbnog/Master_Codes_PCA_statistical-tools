!####################################################################################
! ####################################################################################
subroutine Double_Decomposition(q,nx,ny,nz,nt,qfluc,qmean2D)
  !Funcao que dado uma variavel q conservada "total" qualquer, retorna a flutuacao 
  !e somente a media , nas variaveis qfluc e qmean2D
    
  implicit none
  integer(kind=4), intent(in) :: nx,ny,nz,nt
  real(kind=8),    intent(in) :: q(nx,ny,nz,nt)
  real(kind=8),    intent(out) :: qmean2D(nx,ny),qfluc(nx,ny,nz,nt)
  real(kind=8) :: qmean3D(nx,ny,nz),qdummy(nx,ny)
  integer(kind=4) :: i,j,t

  qmean3D = 0.0d0
  qmean2D = 0.0d0
  
  i = (nx/2) +1     

    qmean3D = sum(q(1:nx,1:ny,1:nz,:),dim=4)/dble(nt)
    qmean2D  = sum(qmean3D,dim=3)/dble(nz)                         
   
  !Calculando a flutuação
  do j = 1,ny
      do i = 1,nx 
        qfluc(i,j,1:nz,:) = q(i,j,1:nz,:) - qmean2D(i,j)        
      enddo
  enddo
  
end subroutine Double_Decomposition
!####################################################################################
subroutine Reynolds_Decomposition(export_flag)

  !se a média for temporal apenas

  use mod_field
  use mod_pod_modes, only : working_zone
  use mod_CGNS
  implicit none
  logical, intent(in) :: export_flag
  
  integer(kind=4) :: j, k, m, t

  integer(kind=4) nx1,nx2,nx3

  write(*,'(A)') ' Computing fluctuation data ...'

  m = working_zone

    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    nx3 = kmax(m)
  
    allocate(zone(m)%qmean3D(imin(m):imax(m),1:jmax(m),1:kmax(m)))
             zone(m)%qmean3D(imin(m):imax(m),1:jmax(m),1:kmax(m)) = 0.0d0

    zone(m)%qmean3D(imin(m):imax(m),1:jmax(m),1:kmax(m)) = &
         sum(zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap),dim=4)/dble(nsnap)

    do t = 1,nsnap
      
      do k = 1,nx3
        do j = 1,nx2
          zone(m)%q(:,j,k,t) = zone(m)%q(:,j,k,t) - zone(m)%qmean3D(:,j,k)
        enddo
      enddo
      
    enddo
    
  if (export_flag .eqv. .false.) return
    
    !statistics file
    call create_file_CGNS('./mean_3D.cgns','3D')
    call write_mean_soln_3D_CGNS(output_zone,[nx1,nx2,nx3],zone(m)%qmean3D(imin(m):imax(m),1:jmax(m),1:kmax(m)),trim(working_var)//'_ave','mean_3D.cgns')
   
    !call write_iblank_2D_CGNS(m,[nx1,nx2],zone(m)%iblank(imin(m):imax(m),1:jmax(m)),'./mean_3D.cgns')
   
    call write_link_CGNS(output_zone,'mean_3D.cgns','grid_mod.cgns')
     
  return  

end subroutine Reynolds_Decomposition
!####################################################################################
