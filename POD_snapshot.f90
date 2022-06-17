subroutine POD_snapshot(flag)

  use mod_pod_modes
  implicit none
  character(len=*), intent(in) :: flag
  
  logical logical_temp
  logical logical_spat
  logical logical_rcns

  logical_temp = .false.
  logical_spat = .false.
  logical_rcns = .false.
  
  select case(flag)
    case ('temporal') 
      logical_temp = .true.
    case ('spatial')
      logical_spat = .true.
    case ('reconstrucao')
      logical_rcns = .true.
    case default
      stop 'Definir "temporal" ou "spatial" ou "reconstrucao" no parameters.in'    
  end select

  ! Compute correlation matrix
  if (logical_temp .eqv. .true.)  then
  
    ! Durr
    call compute_cell_volume(working_zone)
    
    ! The output is the correlation matrix
    call Correlation_Matrix_Classic_POD
  
  endif
  
  ! Compute spatial modes. It requires the pre-computation of the SVD
  if (logical_spat .eqv. .true.)  then
     
    ! The output are the first "nPODmodes" spatial modes
    call Spatial_Modes_Classic_POD(.true.)
  
  endif
  
  ! Reconstruct the flowfield based on the truncated SVD for the first "nPODmodes"
  if (logical_rcns .eqv. .true.)  then

    call Spatial_Modes_Classic_POD(.false.)
     
    call Reconstruct_Flow_Classic_POD
  
  endif  

  return

end subroutine POD_snapshot
!#####################################################################################

!#####################################################################################
subroutine Correlation_Matrix_Classic_POD

  use mod_CGNS, only : working_var
  use mod_pod_modes
  use mod_field
  implicit none
  integer(kind=4) :: j, k, m
  integer(kind=4) :: t1, t2
  
  integer(kind=4) :: nx1, nx2, nx3
  
  real(kind=8), allocatable :: corr_matrix(:,:,:)
  real(kind=8), allocatable :: q1(:)
  real(kind=8), allocatable :: q2(:)
  real(kind=8), allocatable :: wk(:)
  
  character(len=5) :: char_zone

  print*, 'Computing correlation matrix ...'

  call system('date')

  call system('mkdir -p ./SVD_files/')
  open(999,file='./SVD_files/nsnap')
  write(999,*) nsnap
  close(999)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~ Matriz de correlação
  m = working_zone
    
    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    nx3 = kmax(m)
  
    ALLOCATE(corr_matrix(1:nsnap,1:nsnap,0:nx3))
             corr_matrix(1:nsnap,1:nsnap,0:nx3) = 0.0d0
    
    ALLOCATE(q1(1:nx1))
             q1(1:nx1) = 0.0d0
             
    ALLOCATE(q2(1:nx1))
             q2(1:nx1) = 0.0d0

    ALLOCATE(wk(1:nx1))
    
    !$OMP PARALLEL DO PRIVATE(t1,t2,j,k,q1,q2,wk) collapse(2) NUM_THREADS(10) schedule(dynamic,50)
    do t1 = 1,nsnap
            
      do k = 1,nx3
       
        corr_matrix(t1,:,k) = 0.0d0 
        do j = 1,nx2

          if (data_2D .eqv. .true.) wk(1:nx1) = zone(m)%area(imin(m):imax(m),j)
          if (data_3D .eqv. .true.) wk(1:nx1) = zone(m)%volume(imin(m):imax(m),j)
    
          q1(1:nx1) = zone(m)%q(imin(m):imax(m),j,k,t1)
          q1(1:nx1) = q1(1:nx1)*wk(1:nx1)

          do t2 = t1,nsnap
            q2(1:nx1) = zone(m)%q(imin(m):imax(m),j,k,t2)
            corr_matrix(t1,t2,k) = corr_matrix(t1,t2,k) + dot_product(q1(1:nx1),q2(1:nx1))
          enddo
 
        enddo
          
      enddo
      
    enddo
    !$END PARALLEL

    DEALLOCATE(q1)
    DEALLOCATE(q2)
    DEALLOCATE(wk)
  
    do t1 = 1,nsnap  
      do t2 = t1,nsnap
        corr_matrix(t1,t2,0) = sum(corr_matrix(t1,t2,1:nx3))
      enddo
    enddo

  !end of the correlation matrix construction

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~ Exportando a matriz de correlacao
  write(char_zone,'(a1,i4.4)') 'z', output_zone
  if (trim(working_var) .eq.  'Pressure') open(1000, file='./SVD_files/corr_matrix_p_' // char_zone // '.dat')
  
  if (trim(working_var) .eq. 'VelocityX') open(1000, file='./SVD_files/corr_matrix_u_' // char_zone // '.dat')
  if (trim(working_var) .eq. 'VelocityY') open(1000, file='./SVD_files/corr_matrix_v_' // char_zone // '.dat')
  if (trim(working_var) .eq. 'VelocityZ') open(1000, file='./SVD_files/corr_matrix_w_' // char_zone // '.dat')
  
  if (trim(working_var) .eq. 'MomentumX') open(1000, file='./SVD_files/corr_matrix_ru_' // char_zone // '.dat')
  if (trim(working_var) .eq. 'MomentumY') open(1000, file='./SVD_files/corr_matrix_rv_' // char_zone // '.dat')  
  if (trim(working_var) .eq. 'MomentumZ') open(1000, file='./SVD_files/corr_matrix_rw_' // char_zone // '.dat')       
 
  do t1 = 1,nsnap

    !mirroring the correlation matrix
    !DIR$ NOVECTOR
    do t2 = 1,t1-1
      corr_matrix(t1,t2,0) = corr_matrix(t2,t1,0)
    enddo

  enddo

  do t1 = 1,nsnap
    write(1000,'(9999es21.13)') corr_matrix(t1,:,0)
  enddo
  close(1000)
  
  call system('date')  

  return

end subroutine Correlation_Matrix_Classic_POD
!#######################################################################################

!#######################################################################################
subroutine Spatial_Modes_Classic_POD(flag)

  use mod_pod_modes
  use mod_CGNS
  use mod_field
  implicit none
  logical, intent(in) :: flag
  integer(kind=4) :: i, j, k, m, n
  integer(kind=4) :: nx1!, nx2, nx3
  real(kind=8), allocatable :: vdummy(:)
  
  character(len=5) :: char_zone
  character(len=3) :: cdummy
  
  allocate(pod_zone(1:nzones))

  !~~~ Matriz sigma da SVD
  allocate(lambda(nsnap))
  open(10,file='./SVD_files/sigma.dat')
  do i = 1,nsnap
    read(10,*) lambda(i)
  enddo
  close(10)

  !~~~ Matriz V da SVD
  allocate(temporal_modes(1:nsnap,1:nPODmodes))
  allocate(vdummy(1:nsnap))
  open(30,file='./SVD_files/V_matrix.dat')
  do i = 1,nsnap
    read(30,*) vdummy(1:nsnap)
    temporal_modes(i,1:nPODmodes) = vdummy(1:nPODmodes)
  enddo
  close(30)

  !~~~ Cálculo dos modos espaciais
  write(*,*) 'Computing spatial modes ...'

  m = working_zone
  
    nx1 = imax(m)-imin(m)+1
    !nx2 = jmax(m)
    !nx3 = kmax(m)  

    ALLOCATE(pod_zone(m)%spatial_modes(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nPODmodes))
             pod_zone(m)%spatial_modes(:,:,:,:) = 0.0d0

    do n = 1,nPODmodes
      print*, n
      vdummy(1:nsnap) = temporal_modes(1:nsnap,n)
      do k = 1,kmax(m)
        !$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(SHARED) NUM_THREADS(10)
        do j = 1,jmax(m)
          do i = imin(m),imax(m)
            if (zone(m)%iblank(i,j) .ne. 10 .and. iblank .eqv. .true.) zone(m)%q(i,j,k,1:nsnap) = 0.0d0
            pod_zone(m)%spatial_modes(i,j,k,n) = dot_product(zone(m)%q(i,j,k,1:nsnap),vdummy(1:nsnap))/lambda(n)
          enddo
        enddo
        !$OMP END PARALLELDO
      enddo
    enddo
    deallocate(vdummy)
    
    !~~~~~~~~~~~~~~    
    
    if (flag .eqv. .false.) return

    !~~~~~~~~~~~~~~
    
    write(*,'(A)') ' Writing spatial modes ...'

    write(char_zone,'(a1,i4.4)') 'z', output_zone

    if (trim(working_var) .eq. 'MomentumX') open(99,file=trim(output_path)//'modes_'//char_zone//'_rU.dat')
    if (trim(working_var) .eq. 'MomentumY') open(99,file=trim(output_path)//'modes_'//char_zone//'_rV.dat')    
    if (trim(working_var) .eq. 'VelocityX') open(99,file=trim(output_path)//'modes_'//char_zone//'_U.dat')
    if (trim(working_var) .eq. 'VelocityY') open(99,file=trim(output_path)//'modes_'//char_zone//'_V.dat')    
    if (trim(working_var) .eq.  'Pressure') open(99,file=trim(output_path)//'modes_'//char_zone//'_p.dat')

    !~~~~~~~~~~~~~~
    !~~~ Header

      write(99,'(A)') 'VARIABLES='
      write(99,'(A)') '"X"'
      write(99,'(A)') '"Y"'
      if (data_3D .eqv. .true.) write(99,'(A)') '"Z"'
      if (iblank .eqv. .true.) write(99,'(A)') '"Iblank"'
      do n = 1,nPODmodes
        write(cdummy,'(a1,i2.2)') 'M', n
        write(99,'(A)') '"'//cdummy//'"'
      enddo
      write(99,'(3(A,i0),A)') 'ZONE T="'//char_zone//'",I=',nx1,',J=',jmax(m),',K=',kmax(m),',F=BLOCK'
    
    !~~~~~~~~~~~~~~    

    print*, 'X'
    do k = 1,kmax(m); do j = 1,jmax(m); do i = imin(m),imax(m)
      write(99,'(f13.7)') zone(m)%x(i,j,k)
    enddo; enddo; enddo

    print*, 'Y'
    do k = 1,kmax(m); do j = 1,jmax(m); do i = imin(m),imax(m)
      write(99,'(f13.7)') zone(m)%y(i,j,k)
    enddo; enddo; enddo

    if (data_3D .eqv. .true.) then
      print*, 'Z'
      do k = 1,kmax(m); do j = 1,jmax(m); do i = imin(m),imax(m)
        write(99,'(f13.7)') zone(m)%z(i,j,k)
      enddo; enddo; enddo
    endif
    
    if (iblank .eqv. .true.) then
      print*, 'Iblank'
      do k = 1,kmax(m); do j = 1,jmax(m); do i = imin(m),imax(m)
        write(99,'(i0)') zone(m)%iblank(i,j)
      enddo; enddo; enddo
    endif

    do n = 1,nPODmodes
      print*, n
      do k = 1,kmax(m); do j = 1,jmax(m); do i = imin(m),imax(m)
        write(99,'(e15.7)') pod_zone(m)%spatial_modes(i,j,k,n)
      enddo; enddo; enddo
    enddo
    
    !~~~~~~~~~~~~~~

  return

end subroutine Spatial_Modes_Classic_POD
!#######################################################################################

!#######################################################################################
subroutine Reconstruct_Flow_Classic_POD

  use mod_pod_modes
  use mod_CGNS
  use mod_field
  implicit none
  integer(kind=4) :: i, j, k, m, n, t
  integer(kind=4) :: nx1, nx2, nx3
  
  ! character(len=5) :: char_zone
  character(len=32) :: output_file
  
  logical mean_2D
  logical mean_3D
  
  !~~~ Reconstrução do escoamento
  
  write(*,*) 'Reconstructing flow field ...'

  m = working_zone
  
    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    nx3 = kmax(m)

    if (allocated(zone(m)%q) .eqv. .true.) print*, ' DEALLOCATE "q" ...'    
    if (allocated(zone(m)%q) .eqv. .true.) DEALLOCATE(zone(m)%q)

    if (allocated(zone(m)%r) .eqv. .true.) print*, ' DEALLOCATE "r" ...'
    if (allocated(zone(m)%r) .eqv. .true.) DEALLOCATE(zone(m)%r)

    if (allocated(zone(m)%u) .eqv. .true.) print*, ' DEALLOCATE "u" ...'
    if (allocated(zone(m)%u) .eqv. .true.) DEALLOCATE(zone(m)%u)
    
    if (allocated(zone(m)%v) .eqv. .true.) print*, ' DEALLOCATE "v" ...'
    if (allocated(zone(m)%v) .eqv. .true.) DEALLOCATE(zone(m)%v)
    
    if (allocated(zone(m)%w) .eqv. .true.) print*, ' DEALLOCATE "w" ...'
    if (allocated(zone(m)%w) .eqv. .true.) DEALLOCATE(zone(m)%w)

    ALLOCATE(pod_zone(m)%reconstructed(imin(m):imax(m),1:jmax(m),1:kmax(m),1))
             pod_zone(m)%reconstructed(:,:,:,:) = 0.0d0

    print*, ' '
    print*, ' Mean : '
    mean_2D = .false.
    mean_3D = .false.
    if (allocated(zone(m)%qmean2D) .eqv. .true.) mean_2D = .true.
    if (allocated(zone(m)%qmean3D) .eqv. .true.) mean_3D = .true.
    print*, mean_2D, mean_3D
    print*, ' '

    do t = nsnap,nsnap
    
      print*, t
      
      do k = 1,kmax(m)
        do j = 1,jmax(m)
          do i = imin(m),imax(m)
          
            !if (mean_2D .eqv. .true.) pod_zone(m)%reconstructed(i,j,k,1) = 0.0d0
            !if (mean_3D .eqv. .true.) pod_zone(m)%reconstructed(i,j,k,1) = 0.0d0
            if (mean_2D .eqv. .true.) pod_zone(m)%reconstructed(i,j,k,1) = zone(m)%qmean2D(i,j)
            if (mean_3D .eqv. .true.) pod_zone(m)%reconstructed(i,j,k,1) = zone(m)%qmean3D(i,j,k)
            do n = 1,nPODmodes
              pod_zone(m)%reconstructed(i,j,k,1) = pod_zone(m)%reconstructed(i,j,k,1) + &
                                                         temporal_modes(t,n)*lambda(n)*pod_zone(m)%spatial_modes(i,j,k,n)
            enddo
            
          enddo
        enddo
      enddo
    
      write(*,'(A,i0)') ' Writing reconstructed solution file for snapshot ', t
     
      !write(output_file,'(A,i4.4,A)') './rcns', t, '_f.cgns'  
      write(output_file,'(A,i4.4,A)') './rcns', t, '.cgns'  
      call create_file_CGNS(trim(output_file),'3D')
      call write_full_soln_3D_CGNS(m,[nx1,nx2,nx3],pod_zone(m)%reconstructed(:,:,:,1),trim(working_var),trim(output_file))
      call write_link_CGNS(m,trim(output_file),'./grid_mod.cgns')

    enddo

    !~~~~~~~~~~~~~~

  return

end subroutine Reconstruct_Flow_Classic_POD
!#######################################################################################
