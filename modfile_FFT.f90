!##################################################################################################################
!
!   modfile_FFT.f90
!
!   Last modification:
!
!        Qui Mai 24 09:31:43 BRT 2018
!
!
!##################################################################################################################
module mod_signal_process
 
  !real(kind=8), parameter :: max_freq = 10000.0d0
  real(kind=8) :: dt

  real(kind=8), parameter :: overlap = 0.666666666666666d0
  
  !integer(kind=4) :: bin_size
  !integer(kind=4) :: nsnap

  character(len=5) :: char1, char2, char_freq
  
  integer(kind=4) :: nfreq, nyquist
  real(kind=8) :: period, delta, freq
  
  real(kind=8) mean
  integer(kind=4) :: nbins, binsize          !
  !integer(kind=4) :: k1, k2                  !aux variables
  integer(kind=4), allocatable :: bins(:,:)  !range of each bin
  
  real(kind=8) window_cte
  real(kind=8), allocatable :: window(:)

  !FFT
  complex(kind=8), allocatable :: zsignal(:)
  complex(kind=8), allocatable :: zdummy(:)  
  
  !DFT

  real(kind=8), allocatable :: dsignal(:)
  
!  type mode_select
!    !real(kind=8), allocatable, dimension(:,:,:,:) :: fft_mean
!    !real(kind=8), allocatable, dimension(:,:,:,:,:) :: fft_coeff
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: u_hat
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: v_hat
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: w_hat
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: p_hat
!  end type
!  type(mode_select), allocatable :: mode(:)
 
  integer(kind=4) :: nmodes, mode_number
  integer(kind=4), allocatable :: list_of_modes(:)

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~   
    subroutine compute_bins(nmin,skip,nmax,flag)
     
      implicit none
      integer(kind=4), intent(in) :: nmin
      integer(kind=4), intent(in) :: skip
      integer(kind=4), intent(in) :: nmax
      integer(kind=4) :: i, k1, k2 
      logical, intent(in) :: flag

      nbins = 0
      k1 = nmin
      k2 = k1+(binsize-1)

      if (k2 .gt. nmax) then
        write(*,*) ' The number of snapshots is smaller than the bin size...', k2, nmax
        write(*,*) ' The number of snapshots is smaller than the bin size...', k2, nmax
        write(*,*) ' The number of snapshots is smaller than the bin size...', k2, nmax
        write(*,*) ' The number of snapshots is smaller than the bin size...', k2, nmax
        write(*,*) ' The number of snapshots is smaller than the bin size...', k2, nmax
        stop
      endif

      if (flag .eqv. .true.) write(*,*) '  Proposed bins:'
      do while (k2 .le. nmax)
        k2 = k1 + (binsize-1)
        if (flag .eqv. .true.) write(*,'(3x,i2,A,2i6)') nbins+1,')',k1,k2
        if (k2 .le. nmax) then
          nbins = nbins + 1
          k1 = int(dble(nbins)*dble(binsize)*(1.0d0-overlap)) + nmin
        endif
      enddo

      allocate(bins(nbins,1:2))
      do i = 1,nbins
        bins(i,1) = int(dble(i-1)*dble(binsize)*(1.0d0-overlap)) + nmin
        bins(i,2) = bins(i,1) + (binsize-1)
      enddo

      write(*,*) '  Working with the bins:'
      do i = 1,nbins
        write(*,'(3x,i2,A,X,i5.5,A,i5.5)') i, ')', bins(i,1), '.', bins(i,2)
      enddo
      write(*,*) ' '

      return

    end subroutine
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    !Os inputs devem ser entre 1 e N
    !Essa rotina converte para a logica de 0 a N-1
    subroutine compute_dft(omega_index,N,input_signal,summ)

      implicit none
      integer(kind=4), intent(in) :: omega_index
      integer(kind=4), intent(in) :: N
      real(kind=8), intent(in) :: input_signal(0:N-1)
      complex(kind=8), intent(out) :: summ

      integer(kind=4) :: k1, k2

      real(kind=8) :: pi
      complex(kind=8) :: imag

      if (omega_index .eq. 0) stop 'The index of the frequency should be between "1" and "N"'

      pi = dacos(-1.0d0)  
      imag = dcmplx(0.0,1.0)

      summ = dcmplx(0.0d0,0.0d0)
      k1 = omega_index - 1
      do k2 = 0,N-1
        summ = summ + input_signal(k2)*cdexp(-imag*2.0d0*pi/dble(N))**(k1*k2)
      enddo
      summ = summ/dble(N)
      
      return

    end subroutine
    !~~~~~~~~~~~~~~~~~~~~~~~~~
       
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine tukey_wf(N,alfa,window)

    !https://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.signal.tukey.html#scipy.signal.tukey

      implicit none
      integer(kind=4), intent(in) :: N
      real(kind=8), intent(in) :: alfa
      real(kind=8), intent(out) :: window(1:N)
      
      integer(kind=4) :: i, iaux, i1, i2, i3, i4, i5, i6
      real(kind=8) :: daux
      real(kind=8), parameter :: pi = dacos(-1.0d0)  

      daux = alfa*dble(N)/2.0d0
      iaux = int(daux)
      
      i1 = 0
      i2 = iaux
      i3 = iaux
      i4 = int(dble(N-1)*(1.0d0-(alfa/2.0d0)))
      i5 = int(dble(N-1)*(1.0d0-(alfa/2.0d0)))
      i6 = N
      
      !DIR$ NOVECTOR
      do i = 0,N-1
        if (i .ge. i1 .and. i .lt. i2) then
          window(i+1) = 0.5d0*( 1.0d0 + dcos(pi*(dble(i)/daux - 1.0d0)) )
        endif
        if (i .ge. i3 .and. i .le. i4) then
          window(i+1) = 1.0d0
        endif
        if (i .gt. i5 .and. i .lt. i6) then 
          window(i+1) = 0.5d0*( 1.0d0 + dcos(pi*(dble(i)/daux - 2.0d0/alfa + 1.0d0)) )
        endif
      enddo

      call system('printf "\033[1;33m Using TUKEY window function for the dipole sources!!"')             
      write(*,'(A,f5.2)') '  - Tukey window function parameter = ', alfa
      call system('"\033[0m"')

      return

    end subroutine tukey_wf
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine hamming_wf(N,window)

    !The window with these particular coefficients was proposed by Richard W. Hamming.
    !The window is optimized to minimize the maximum (nearest) side lobe, giving it a height of about one-fifth that of the Hann window

    !https://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.signal.hamming.html#scipy.signal.hamming

      implicit none
      integer(kind=4), intent(in) :: N
      real(kind=8), intent(out) :: window(1:N)
      
      integer(kind=4) :: i
      real(kind=8), parameter :: alfa = 0.54d0
      real(kind=8), parameter :: beta = 0.46d0
    !  real(kind=8), parameter :: alfa = 25.0d0/46.0d0
    !  real(kind=8), parameter :: beta = 21.0d0/46.0d0  
      real(kind=8), parameter :: pi = dacos(-1.0d0)  
      
      do i = 0,N-1
        window(i+1) = alfa - beta*dcos(2.0d0*pi*dble(i)/dble(N)) 
      enddo

      call system('printf "\033[1;33m Using HAMMING window function for the dipole sources!! \033[0m"')             
      
      return

    end subroutine hamming_wf
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine retangular_wf(N,window)
    
      implicit none
      integer(kind=4), intent(in) :: N
      real(kind=8), intent(out) :: window(1:N)    
    
      integer(kind=4) :: i
          
      do i = 1,N-1
        window(i) = 1.0d0
      enddo
      window(N) = 0.0d0

      write(*,*) ' '      
      call system('printf "\033[0;32m WARNING: No window function used for the dipole sources!! \033[0m"') 
      write(*,*) ' '
      call system('printf "\033[0;32m WARNING: No window function used for the dipole sources!! \033[0m"')
      write(*,*) ' '
      
      return
      
    end subroutine retangular_wf
    !~~~~~~~~~~~~~~~~~~~~~~~~~      
    
end module mod_signal_process
!##################################################################################################################
