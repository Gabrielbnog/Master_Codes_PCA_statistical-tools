
!####################################################################################
subroutine compute_spatial_fourier_modes_mean

  use mod_field
  use mod_signal_process
  use mod_pod_modes, only : working_zone, nFourierModes
  use mod_CGNS  
  implicit none
  complex(kind=8) :: coeff
  
  character(len=3) :: char_mode
  character(len=4) :: char_zone
  
  integer(kind=4) :: i, j, k, m, n, t
  integer(kind=4) :: nx1, nx2
  
  real(kind=8), allocatable :: fourier_coeffs(:,:,:)
   
  m = working_zone
  
  nFourierModes = 16
  
  if (m .eq. 1) then
  
    do n = 1,nFourierModes
      
      write(char_mode,'(i3.3)') n
  
      open(666,file='file_'//char_mode//'.dat')
      write(666,'(A)') 'VARIABLES="X","Y","Q"'
      close(666)
      
    enddo
    
  endif    
  
  !write(*,'(A,i0,A)') ' Working with spanwise FOURIER POD for ', nFourierModes, ' modes !!!'


    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)

    binsize = kmax(m)
    allocate(dsignal(1:binsize))

    allocate(fourier_coeffs(1:nx1,1:nx2,1:nFourierModes))
             fourier_coeffs(1:nx1,1:nx2,1:nFourierModes) = -1234567.0d0
    
    do n = 1,nFourierModes

      mode_number = n
      
      do j = 1,nx2
        do i = 1,nx1

          dsignal(1:binsize) = zone(m)%qmean3D(i,j,1:kmax(m))
          call compute_dft(mode_number,binsize,dsignal,coeff)
          fourier_coeffs(i,j,n) = cdabs(coeff)

        enddo
      enddo
      
      write(char_mode,'(i3.3)') n
      write(char_zone,'(i4.4)') working_zone
      open(666,file='file_'//char_mode//'.dat',position='append')
      write(666,'(2(A,i0),A,f10.4)') 'ZONE T="'//char_zone//'",I=',nx1,',J=',nx2,',F=POINT, STRANDID=1, SOLUTIONTIME=', dble(n)   
      k = 1
      do j = 1,nx2
        do i = 1,nx1
          write(666,'(3e23.15)') zone(m)%x(i,j,k), zone(m)%y(i,j,k), fourier_coeffs(i,j,n)
        enddo
      enddo
      close(666)   
      
    enddo

  return

end subroutine compute_spatial_fourier_modes_mean
!####################################################################################
