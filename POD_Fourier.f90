!####################################################################################
subroutine POD_Fourier(idummy)

  use mod_pod_modes
  implicit none
  integer(kind=4), intent(in) :: idummy

  nFourierModes = 3

  ! Durr
  call compute_cell_volume(working_zone)
   
  ! Fourier transform in the spanwise direction
  call compute_spatial_fourier_modes

  ! The output is the correlation matrix
  call Correlation_Matrix_Fourier_POD

  return

end subroutine
!####################################################################################

!####################################################################################
subroutine compute_spatial_fourier_modes

  use mod_field
  use mod_signal_process
  use mod_pod_modes, only : working_zone, nFourierModes
  use mod_CGNS  
  implicit none
  complex(kind=8) :: coeff
  
  integer(kind=4) :: i, j, m, n, t
  integer(kind=4) :: nx1, nx2
   
  m = working_zone
  
  write(*,'(A,i0,A)') ' Working with spanwise FOURIER POD for ', nFourierModes, ' modes !!!'

    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)

    binsize = kmax(m)
    allocate(dsignal(1:binsize))

    allocate(zone(m)%q_hat(1:nx1,1:nx2,1:nFourierModes,1:nsnap))
    
    do n = 1,nFourierModes

      mode_number = n
      
      do t = 1,nsnap
        do j = 1,nx2
          do i = 1,nx1

            dsignal(1:binsize) = zone(m)%q(i,j,1:kmax(m),t)
            call compute_dft(mode_number,binsize,dsignal,coeff)
            zone(m)%q_hat(i,j,n,t) = coeff

          enddo
        enddo
      enddo
      
    enddo

  return

end subroutine compute_spatial_fourier_modes
!####################################################################################

!####################################################################################
subroutine Correlation_Matrix_Fourier_POD

  use mod_pod_modes
  use mod_field
  implicit none  
!  integer(kind=4) :: i, j, k, l, m
  integer(kind=4) :: k, m
  integer(kind=4) :: t1, t2
  
  integer(kind=4) :: nx1, nx2, nx3
  
  complex(kind=8), allocatable :: corr_matrix(:,:)
  complex(kind=8), allocatable :: q1(:)
  complex(kind=8), allocatable :: q2(:)
  real(kind=8), allocatable :: wk(:)
  
  character(len=5) :: char_zone

  call system('date')

  print*, 'Computing correlation matrix ...'

  allocate(corr_matrix(1:nsnap,1:nsnap))
  corr_matrix(1:nsnap,1:nsnap) = 0.0d0

  call system('mkdir -p ./SVD_files/')
  open(999,file='./SVD_files/nsnap')
  write(999,*) nsnap
  close(999)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !do m = 1,nzones
  m = working_zone
    
    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    nx3 = kmax(m)
    ALLOCATE(q1(nx1*nx2))
    ALLOCATE(q2(nx1*nx2))

    ALLOCATE(wk(nx1*nx2))
    wk = reshape(zone(m)%area(imin(m):imax(m),1:jmax(m)),[nx1*nx2])
    
   !XXX !$OMP PARALLEL DO PRIVATE(t1,t2,k,q1,q2) collapse(2) NUM_THREADS(20)
    do t1 = 1,nsnap
            
      do k = 1,nFourierModes
    
        q1(1:nx1*nx2) = reshape( zone(m)%q(imin(m):imax(m),1:jmax(m),k,t1), [nx1*nx2] )
      
        q1(:) = q1(:)*wk(:)

        do t2 = t1,nsnap
          
          q2(1:nx1*nx2) = reshape( zone(m)%q(imin(m):imax(m),1:jmax(m),k,t2), [nx1*nx2] )
          q2(1:nx1*nx2) = dconjg( q2(1:nx1*nx2) )

          corr_matrix(t1,t2) = dot_product(q1(:),q2(:))

        enddo
        
      enddo
                 
    enddo
   !XXX !$END PARALLEL

    DEALLOCATE(q1)
    DEALLOCATE(q2)
    DEALLOCATE(wk)

  !enddo

  !end of the correlation matrix construction

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  write(char_zone,'(a1,i4.4)') 'z', output_zone
!  exportando a matriz de correlacao
  if (.true.) open(1000, file='./SVD_files/corr_matrix_v_' // char_zone // 'R.dat')  !FIXME - Especficar o modo no nome do arquivo!
  if (.true.) open(1001, file='./SVD_files/corr_matrix_v_' // char_zone // 'I.dat')  !FIXME - Especficar o modo no nome do arquivo!
  
  do t1 = 1,nsnap

    !mirroring the correlation matrix
    !DIR$ NOVECTOR
    do t2 = 1,t1-1
      corr_matrix(t1,t2) = corr_matrix(t2,t1)
    enddo

  enddo

  do t1 = 1,nsnap

    write(1000,'(9999es21.13)')  dble( corr_matrix(t1,:) )
    write(1001,'(9999es21.13)') dimag( corr_matrix(t1,:) )

  enddo
  
  close(1000)
  close(1001)
  
  call system('date')  

  return

end subroutine Correlation_Matrix_Fourier_POD
!################################################################################################
