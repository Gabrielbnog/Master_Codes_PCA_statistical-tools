!####################################################################################
program xPOD3D

  use mod_CGNS, only : working_var
  use mod_field, only : idxi, idxr, idxf, nsnap, output_path
  use mod_pod_modes
  use mod_signal_process
  implicit none
 
  ! Set-up of the code
  call data_setup
  
  call read_grid
  
  ! Change the path where is the files will be outputed
  call change_path(trim(output_path))
  
  call read_mean_soln
      
  call compute_spatial_fourier_modes_mean
 
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
  character(len=16) :: var_dummy

  namelist /WORK_PAR/ working_zone, output_zone, working_var, start_index_i, final_index_i, output_path
  namelist /GRID_PAR/ grid_name, path_to_grid, iblank
  namelist /SOLN_PAR/ path_to_soln

  open(610,file='parameters.in',status='old')
    read(610,WORK_PAR)
    read(610,GRID_PAR)
    read(610,SOLN_PAR)
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
 
!    print*, ' (1)'
 
    open(610,file='working_var.in')
    read(610,*) working_var
    close(610)
!    print*, ' (2)'
!    working_var = trim(var_dummy)
!    print*, ' (3)'

    write(*,'(A,A)') ' Forcing new variable :: ', working_var
    write(*,'(A,A)') ' Forcing new variable :: ', working_var
    
    call system('rm -f working_var.in')

  endif
 
  call system('printf "\033[0m \n"')
  
  !~~~ 

  if (trim(working_var) .eq. "") stop ' Juquinha! Especifique uma variavel!! '
  print*, ' ', trim(working_var)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  gamma = 1.4d0
  gas_cte = 287.0d0
  Twall = 2.5d0 ! = (1.0d0/(gamma-1.0d0))
  
  dt = dt*dble(fresult)*dble(idxr)
  
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

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
