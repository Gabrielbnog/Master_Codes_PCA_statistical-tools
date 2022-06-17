subroutine POD_harmonic

  use mod_pod_modes
  use mod_field, only : nsnap
  use mod_signal_process
  implicit none

!~~~~~~~~~~~~~~
!~ One must look on a different code to see what is the index of the Fourier mode associated to the desired frequency
!~ Current possibilities:

!~~~    binsize = 64
!~~~    mode_number = 5
      
!~~~    binsize = 128
!~~~    mode_number = 9
      
!~~~    binsize = 256
!~~~    mode_number = 19
      
!~~~    binsize = 512
!~~~    mode_number = 37
    
      call compute_cell_volume(working_zone)
      
      ! Setando os limites de cada bin da transformada de Fourier
      call compute_bins(1,1,nsnap,.false.)
      
      ! << true >> for Harmonic "Towne's" Spectral POD
      !    -> exports the Fourier-transformed solution matrix for SVD computation
      ! if false, the result can be used for something different than SpectralPOD
      call compute_temporal_fourier_modes(.true.)
      
      ! And then skip the correlation matrix because it is exported ...
      
      ! The spatial modes U are directly obtained by the SVD
      ! The "temporal" modes V give the phase shift (I believe so)
      ! The singular values S give the magnitude weighting

  return

end subroutine
!####################################################################################  

!####################################################################################
subroutine compute_temporal_fourier_modes(flag)

  use mod_field
  use mod_signal_process
  use mod_pod_modes, only : working_zone
  use mod_CGNS  
  implicit none
  logical, intent(in) :: flag

  complex(kind=8) :: coeff
  
  integer(kind=4) :: i, j, k, k1, k2, l, m
  !integer(kind=4) :: idummy
  integer(kind=4) :: nx1, nx2, nx3

  character(len=10) :: dft_file
  character(len=64) :: pod_file

  allocate(dsignal(1:binsize))
  
  m = working_zone
  
    write(*,'(A,i0,A)') ' Working with Harmonic "Spectral" POD for mode ', mode_number, ' !!!'
    print*, 2.0d0*3.141592d0*dble(mode_number-1)*(1.0d0/(dble(binsize)*dt))

    nx1 = imax(m)
    nx2 = jmax(m)
    nx3 = kmax(m)
    allocate(zone(m)%q_hat(1:nx1,1:nx2,1:nx3,1:nbins))
             zone(m)%q_hat(1:nx1,1:nx2,1:nx3,1:nbins) = dcmplx(0.0d0,0.0d0)

    do l = 1,nbins

      k1 = bins(l,1) - bins(1,1) + 1
      k2 = bins(l,2) - bins(1,1) + 1

      write(*,'(A,i5.5,A,i5.5)') ' --> Working on bins : ', k1, '.', k2

      do k = 1,nx3
        do j = 1,nx2
          do i = 1,nx1

            dsignal(1:binsize) = zone(m)%q(i,j,k,k1:k2)
            call compute_dft(mode_number,binsize,dsignal,coeff)
            zone(m)%q_hat(i,j,k,l) = coeff

          enddo
        enddo
      enddo
      
      if (.false.) then
      
        write(*,'(A)') ' --> Exporting Fourier modes ... '

        write(dft_file,'(a,i2.2,a5)') 'dft', l,'.cgns' 
      
        if (data_2D .eqv. .true.) then
          call create_file_CGNS(trim(dft_file),'2D')
          call write_fourier_modes_2D_CGNS(m,[nx1,nx2],l,zone(m)%q_hat(1:nx1,1:nx2,1,l),trim(working_var), &
                                   trim(output_path)//trim(dft_file))
          call write_link_CGNS(m,trim(dft_file),'./grid_mod_2D.cgns')
        endif
        
        if (data_3D .eqv. .true.) then
          call create_file_CGNS(trim(dft_file),'3D')
          call write_fourier_modes_3D_CGNS(m,[nx1,nx2,nx3],l,zone(m)%q_hat(1:nx1,1:nx2,1:nx3,l),trim(working_var), &
                                   trim(output_path)//trim(dft_file))
          call write_link_CGNS(m,trim(dft_file),'./grid_mod.cgns')
        endif

      endif

    enddo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    !
    !~~~ Exportar a matriz para o SVD
    if (flag .eqv. .true.) then

      call system('printf "\033[0;33m \n"')
      write(*,'(A)') ' Exporting Fourier filtered data for Harmonic Spectral POD ... '
      call system('printf "\033[0m \n"')    
    
      call system('mkdir -p output_fourier_modes/')
      
      !~~~ Exportando a matriz em BINARIO para fazer o SVD pelo Python
      ! Eu ignoro o ultimo ponto da malha pois a concatenacao fica automatica (eu acho)
      ! A concatenacao eh necessaria em casos que a malha eh quebrada em diversas zonas
      ! No caso de uma zona apenas, eliminar o ultimo ponto tambem eh necessario pois ele eh o ponto 1 replicado
!      write(pod_file,'(a,i2.2,a4)') 'output_fourier_modes/fourier_matrix', l,'.dat'
!      write(123) nx1-1, nx2, nx3, nbins !XXX ignorar o ultimo facilita a concatenacao
!      do k = 1,nx3
!        do j = 1,nx2
!          do i = 1,nx1-1                !XXX ignorar o ultimo facilita a concatenacao
!            write(123) zone(m)%q_hat(i,j,k,:)*dsqrt(zone(m)%area(i,j))
!          enddo
!        enddo
!      enddo
!      idummy = 0
!      write(123) idummy
!      close(123)

      open(124,file='output_fourier_modes/header.dat',position='append')
      write(124,*) m, nx1-1, nx2, nx3, nbins !XXX ignorar o ultimo facilita a concatenacao
      !write(124,*) m, nx1, nx2, nx3, nbins
      close(124)

      write(pod_file,'(A,i2.2,A4)') 'output_fourier_modes/fourier_matrix_', m,'.dat'
      open(123,file=trim(pod_file))
      do k = 1,nx3
        do j = 1,nx2
          do i = 1,nx1-1 !XXX ignorar o ultimo facilita a concatenacao
          !do i = 1,nx1
            write(123,'(9999e23.15)') zone(m)%q_hat(i,j,k,:)!*dsqrt(zone(m)%area(i,j))
          enddo
        enddo
      enddo
      !idummy = 0
      !write(123,*) idummy
      close(123)
      
    endif

  return

end subroutine compute_temporal_fourier_modes
!####################################################################################
