!####################################################################################
program xPOD3D

  !use mod_CGNS!, only : working_var
  use mod_field  !, only : time, idxi, idxr, idxf, nsnap, output_path
  use mod_pod_modes
  use mod_signal_process
  
  implicit none
  integer(kind=4)     :: m !Working zone
  ! Set-up of the code
  call data_setup
  
  ! Change the path where is the files will be outputed
  if (.TRUE.) call change_path(trim(output_path))
  !Reading the grid data for the working zone only !
  call read_grid
  
  !call compute_cell_area(working_zone)
  ! Reads the variables "rho", "u", "v" and "w" based on "working_var"
  if (.TRUE.) call read_soln_ruvw
  
  ! Convert from momentum to velocity depending on the working variables
  if (.FALSE.) call convert_momentum_to_velocity
      
  ! Read variables from .dat file in ASCII format 
  ! Put values in nx, ny, nz, nsnap, imin, imax, jmin, jmax, kmax 
  ! The mean can be considered as 3D or 2D (averaging in "Z" axis as well)
  ! The mean data is computed for the variable "q", which can be anything..
  ! To compute the mean for 'n' different variables, this routine should be called 'n' times
  if (.FALSE.) call Reynolds_Decomposition(.false.)   !arquivos 3D apenas
  
  if (.FALSE.) call Double_Decomposition(.false.)     !arquivos 2D ou 3D...
  
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !  -> The steps below are not required for POD and its variants ...
    !   They are optional and may be used for further analysis of the data
    if (.TRUE.) then
      ! To compute the aerodynamic coefficients (Cp, Cl, Cd, Cf)
      if (.false.) then
        call compute_aerodynamic_coefficients
      endif
      if (.TRUE.) then
        
        m = working_zone
        !Dimensoes do trecho a ser considerado para a zona CGNS
        x1 = imin(m)
        x2 = imax(m)
        y1 = jmin(m) 
        y2 = jmax(m)
        z1 = kmin(m) 
        z2 = kmax(m)
        nt = nsnap
                
        allocate(dqsidx(x1:x2,y1:y2,z1:z2),dqsidy(x1:x2,y1:y2,z1:z2),dqsidz(x1:x2,y1:y2,z1:z2))                                                                      
        allocate(detadx(x1:x2,y1:y2,z1:z2),detady(x1:x2,y1:y2,z1:z2),detadz(x1:x2,y1:y2,z1:z2)) 
        allocate(dphidx(x1:x2,y1:y2,z1:z2),dphidy(x1:x2,y1:y2,z1:z2),dphidz(x1:x2,y1:y2,z1:z2)) 
        allocate(dxdqsi(x1:x2,y1:y2,z1:z2),dydqsi(x1:x2,y1:y2,z1:z2),jacobian(x1:x2,y1:y2,z1:z2),yplus(y2))

        call metric(zone(m)%x(x1:x2,y1:y2,z1:z2),zone(m)%y(x1:x2,y1:y2,z1:z2),zone(m)%z(x1:x2,y1:y2,z1:z2),&
        dqsidx(x1:x2,y1:y2,z1:z2),dqsidy(x1:x2,y1:y2,z1:z2),dqsidz(x1:x2,y1:y2,z1:z2),&
        detadx(x1:x2,y1:y2,z1:z2),detady(x1:x2,y1:y2,z1:z2),detadz(x1:x2,y1:y2,z1:z2),&
        dphidx(x1:x2,y1:y2,z1:z2),dphidy(x1:x2,y1:y2,z1:z2),dphidz(x1:x2,y1:y2,z1:z2),&
        dxdqsi(x1:x2,y1:y2,z1:z2),dydqsi(x1:x2,y1:y2,z1:z2),jacobian(x1:x2,y1:y2,z1:z2),x2-x1+1,y2,z2)
        
        call velocity_profile(x2-x1+1,y2,z2,nt,zone(m)%x(x1:x2,y1:y2,z1:z2),zone(m)%y(x1:x2,y1:y2,z1:z2),Re,Ma, &                                                                                                                                       
         zone(m)%u(x1:x2,1:y2,1:z2,1:nt),zone(m)%v(x1:x2,1:y2,1:z2,1:nt),zone(m)%r(x1:x2,1:y2,1:z2,1:nt), &
         utau,mu,dens,yplus(1:y2)) 
        
        ! call Lumley(x2-x1+1,y2,z2,nt,zone(m)%./x  u(x1:x2,1:y2,1:z2,1:nt),zone(m)%v(x1:x2,1:y2,1:z2,1:nt),& 
        !  zone(m)%w(x1:x2,1:y2,1:z2,1:nt),yplus(1:y2),detadx(x1:x2,y1:y2,z1:z2),detady(x1:x2,y1:y2,z1:z2),&
        !  dqsidx(x1:x2,y1:y2,z1:z2),dqsidy(x1:x2,y1:y2,z1:z2))        
         
        ! call ReynoldsStress_profile(x2-x1+1,y2,z2,nt,zone(m)%x(x1:x2,y1:y2,z1:z2), &                                           
        !  zone(m)%u(x1:x2,1:y2,1:z2,1:nt),zone(m)%v(x1:x2,1:y2,1:z2,1:nt),zone(m)%w(x1:x2,1:y2,1:z2,1:nt), &                                                     
        !  yplus(1:y2),utau)     
         
        call probe_along_Z(x2-x1+1,y2,z2,nt,zone(m)%u(x1:x2,1:y2,1:z2,1:nt),zone(m)%v(x1:x2,1:y2,1:z2,1:nt),&
         zone(m)%w(x1:x2,1:y2,1:z2,1:nt),detadx(x1:x2,y1:y2,z1:z2),detady(x1:x2,y1:y2,z1:z2),&
         dqsidx(x1:x2,y1:y2,z1:z2),dqsidy(x1:x2,y1:y2,z1:z2),position,yplus(1:y2))
         
        !  print*, ' ' 
        !  print*, 'Computing TKE'
        !  print*, ' '

        ! call BalanceTKE(yplus(1:y2),mu,utau,dens,x2-x1+1,y2,z2,nt,zone(m)%u(x1:x2,1:y2,1:z2,1:nt),zone(m)%v(x1:x2,1:y2,1:z2,1:nt),&
        !  zone(m)%w(x1:x2,1:y2,1:z2,1:nt),zone(m)%p(x1:x2,1:y2,1:z2,1:nt),&
        !  dqsidx(x1:x2,y1:y2,z1:z2),dqsidy(x1:x2,y1:y2,z1:z2),dqsidz(x1:x2,y1:y2,z1:z2),&
        !  detadx(x1:x2,y1:y2,z1:z2),detady(x1:x2,y1:y2,z1:z2),detadz(x1:x2,y1:y2,z1:z2),&
        !  dphidx(x1:x2,y1:y2,z1:z2),dphidy(x1:x2,y1:y2,z1:z2),dphidz(x1:x2,y1:y2,z1:z2))
                           
      endif 
      ! OPTIONAL: 
      ! To get the temporal fluctuation data at a single point (may be specified several times)
      if (.FALSE.) then

        if (working_zone .eq. 1) call probe(112,100,1,working_zone,nt,dt,'./')
        if (working_zone .eq. 9) call probe(114,426,1,working_zone,nt,dt,'./')

      endif

      ! OPTIONAL: 
      ! To get the temporal fluctuation data at a circle
      ! Not to be used in the Bhaskaran code !!
      if (.FALSE.) then

        call circular_probe(500,1,2,nt,dt,'./')

      endif

      ! OPTIONAL: 
      ! To get the temporal fluctuation data at a line along Z 
      ! It may be specified several times...
      if (.FALSE.) then

        write(*,*) ' Exporting probe data'
        if (working_zone .eq. 1) call probe_along_Z(1,35,working_zone,nt,dt,'./')
      
      endif

    endif

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 
    !  -> Proper Orthogonal Decomposition using the "Classic" Snapshot Method
    ! 
    !   It is computationally necessary to compute the correlation matrix
    !   Then, it may be modified to Sieber's Spectral POD in another code
    !
    if (.FALSE.) then
      print*, ''
      print*, trim(POD_operation)
      print*, trim(POD_operation)
      print*, trim(POD_operation)        
      call POD_snapshot(trim(POD_operation))
    endif
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 
    !  -> "Spatial Fourier" Proper Orthogonal Decomposition using the "Classic" Snapshot Method
    ! 
    !   It is necessary to compute the correlation matrix. Note that it is a COMPLEX matrix now!!
    !   The SVD can be computed for each fourier mode separately
    !
    if (.FALSE.) call POD_Fourier(nmodes)
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 
    !  -> Harmonic "Towne's" Spectral POD
    !
    !   In this case, it is NOT NECESSARY to compute the Correlation matrix.
    !   Instead, the temporal fourier modes are directly exported.
    !   Then, the SVD is computed for the "temporal fouried filtered data".
    !   This way, the "temporal modes -- associated to the phase shift" and the spatial modes are directly computed
    !
    !   The paper from Towne ( )
    !     says to work with the correlation matrix. However, the results are the same and less programming is necessary!!
    !
    if (.FALSE.) call POD_harmonic

end program
!####################################################################################

!####################################################################################
subroutine data_setup

  use mod_CGNS
  use mod_field
  use mod_pod_modes
  use mod_signal_process
  implicit none  
  integer(kind=4) :: m
  logical batch_script

  integer(kind=4) :: aux(3)

  integer(kind=4)  :: idummy
  character(len=1) :: cdummy 
  ! character(len=16) :: var_dummy

  namelist /WORK_PAR/ working_zone, output_zone, working_var, start_index_i, final_index_i, output_path
  namelist /GRID_PAR/ grid_name, path_to_grid, iblank
  namelist /SOLN_PAR/ idxi, idxf, idxr, path_to_soln, Re, Ma, dt, fresult,position
  namelist  /POD_PAR/ nPODmodes, POD_operation
  namelist /SPOD_PAR/ mode_number, binsize

  open(610,file='parameters.in',status='old')
    read(610,WORK_PAR)
    read(610,GRID_PAR)
    read(610,SOLN_PAR)
    read(610, POD_PAR)
    read(610,SPOD_PAR)
  close(610)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  call system('printf "\033[0;36m \n"')

  batch_script = .false.
  inquire(file='working_zone.in',exist=batch_script)
  if (batch_script .eqv. .true.) then
  
    open(610,file='working_zone.in')
    read(610,*) working_zone
    close(610)
    output_zone = working_zone
    write(*,'(A,i2)') ' Forcing new zone :: ', working_zone
    write(*,'(A,i2)') ' Forcing new zone :: ', working_zone
    
    call system('rm -f working_zone.in')

  endif

  m = working_zone
  
  call system('printf "\033[0m \n"')  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  call system('printf "\033[0;37m \n"')
  
  batch_script = .false.
  inquire(file='working_var.in',exist=batch_script)
  if (batch_script .eqv. .true.) then
 
    open(610,file='working_var.in')
    read(610,*) working_var
    close(610)

    write(*,'(A,A)') ' Forcing new variable :: ', working_var
    write(*,'(A,A)') ' Forcing new variable :: ', working_var
    
    call system('rm -f working_var.in')

  endif
 
  call system('printf "\033[0m \n"')
  
  if (trim(working_var) .eq. "") stop ' Juquinha! Especifique uma variavel!! '
  print*, ' ', trim(working_var)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
!  dt = dt*dble(fresult)*dble(idxr)
! 
!  allocate(time(1:nsnap))
!  time(1:nsnap) = 0.0d0
!  do i = 1,nsnap
!    time(i) = dble(i-1)*dt
!    print*, time(i)
!  enddo
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~ Malha
  idummy = len(trim(path_to_grid))
  cdummy = path_to_grid(idummy:idummy)
  if (cdummy .ne. '/') then
    path_to_grid(idummy+1:idummy+1) = '/'
  endif

  write(*,'(A)') ' Grid file:'
  !write(*,'(A)') trim(grid_name)
  write(CGNS_gridname,'(A)') trim(path_to_grid)//trim(grid_name)
  write(*,'(x,A)') trim(CGNS_gridname)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~ Solucao
  idummy = len(trim(path_to_soln))
  cdummy = path_to_soln(idummy:idummy)
  if (cdummy .ne. '/') then
    path_to_soln(idummy+1:idummy+1) = '/'
  endif 
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~ Output
  idummy = len(trim(output_path))
  cdummy = output_path(idummy:idummy)
  if (cdummy .ne. '/') then
    output_path(idummy+1:idummy+1) = '/'
  endif   
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  call zoneinfo(trim(CGNS_gridname),nzones)
  print*, ' Number of zones = ', nzones
  if (working_zone .gt. nzones) stop ' --> working_zone .gt. nzones <-- '
 
  allocate(zone(working_zone:working_zone))

  m = working_zone
  call read_size_CGNS(m,trim(CGNS_gridname),aux)
  zone(m)%nx = aux(1)
  zone(m)%ny = aux(2)
  zone(m)%nz = aux(3)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~ Final tweaks in the setup...   

  !~~~ In case anything goes wrong...

  call system('printf "\033[0;31m \n"')

  if (zone(m)%nz .lt. 1) then
    write(*,'(A)') ' FUDEU: nz < 1 !!!'
    call system('printf "\033[0m \n"')
    stop 
  endif
  
  allocate(imin(working_zone:working_zone))
  allocate(imax(working_zone:working_zone))
  allocate(jmin(working_zone:working_zone))
  allocate(jmax(working_zone:working_zone))
  allocate(kmin(working_zone:working_zone))
  allocate(kmax(working_zone:working_zone))
  
  imin(m) = start_index_i
  imax(m) = final_index_i
    
  jmin(:) = 1 !deixar sempre 1...
  jmax(:) = 0 !pode mudar aqui ...
 
  !~~~ 2D case
  if (zone(m)%nz .eq. 1) then

    write(*,'(A,i0)') ' --> Os dados sao  bi-dimensionais (2D) : ', zone(m)%nz

    data_2D = .true.
    data_3D = .false.

    kmin(:) = 1 !deixar sempre 1...
    kmax(:) = 1 !deixar sempre 1... 

  endif
  
  !~~~ 3D case
  if (zone(m)%nz .gt. 1) then      

    write(*,'(A,i0)') ' --> Os dados sao tri-dimensionais (3D) : ', zone(m)%nz
    
    data_2D = .false.
    data_3D = .true.

    kmin(:) = 1 !deixar sempre 1 ...
    kmax(:) = 0 !deixar sempre 0 ...
    
    write(*,'(A)') ' --> Check the 3D setup in "data_setup" in "main.f90" ...'
    allocate(dz(nzones))  
    dz = [1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]

  endif  
 
  call system('printf "\033[0m \n"')
  
  return
  
end subroutine
!####################################################################################

!####################################################################################
subroutine change_path(path)

  implicit none
  character(len=*), intent(in) :: path
  character(len=132) aux_path
  character(len=132) current_path
  
  integer(kind=4)  :: idummy
  character(len=1) :: cdummy

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Changing work path
  aux_path = trim(path)
  
  idummy = len(trim(aux_path))
  cdummy = aux_path(idummy:idummy)
  if (cdummy .ne. '/') then
    aux_path(idummy+1:idummy+1) = '/'
  endif  
  
  call system('mkdir -p ' // trim(aux_path))
  
  call chdir(trim(aux_path))
  call getcwd(current_path)
  write(*,'(3A)') ' Changing working directory to ', trim(current_path), ' ...'

  return

end subroutine
!####################################################################################

