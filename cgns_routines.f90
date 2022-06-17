
  !FIXME : o acesso a zonas já existentes está cagado...
  !FIXME : comparar "ngrid" com "isize" aumenta robustez!

! Jean Marques
! 08/08/2016
!
! Tulio R.
! Qua Set 19 10:11:09 -03 2018

!****************************************************************************
!
!     LIST OF SUBROUTINES :   
!
!
!
!
!
!

!############################################################################
module mod_CGNS

  character(len=16)  :: working_var

  character(len=132) :: CGNS_filename
  character(len=132) :: CGNS_gridname
  character(len=132) :: CGNS_solnname

  integer(kind=4)    :: index_grid, index_soln, index_file, index_mean
  integer(kind=4)    :: index_base, index_zone, index_node, index_flow, index_mode
  integer(kind=4)    :: index_field, index_coord
  
  !integer(kind=4)    :: index_hole
 
  integer(kind=4)    :: ier
  integer(kind=4)    :: location
  integer(kind=4)    :: cell_dim, phys_dim

  integer(kind=4), allocatable :: isize(:,:)

  character(len=16)  :: basename
  character(len=32)  :: zonename
  character(len=32)  :: solnname
  character(len=32)  :: fieldname
  character(len=5)   :: citer
  character(len=2)   :: cbin

  character(len=5)   :: cmode
  
  logical :: CGNS_file_exists

end module
!############################################################################

!------------------------------------------------------------------------------------
! Reads Number of Zones
subroutine zoneinfo(char_filename,nzones)

  use cgns
  use mod_CGNS
  implicit none
  character(len=*), intent(in) :: char_filename
  integer(kind=4), intent(out) :: nzones

  CGNS_filename = trim(char_filename)
  inquire(file=trim(CGNS_filename),exist=CGNS_file_exists)
    if (CGNS_file_exists .eqv. .false.) then
      print*, trim(CGNS_filename)
      stop ' O arquivo CGNS não existe!!!'
    endif

  !... open CGNS file
  call cg_open_f(CGNS_filename,CG_MODE_READ,index_file,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  !... read zone
  index_base = 1
  call cg_nzones_f(index_file,index_base,nzones,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_file,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  return

end subroutine zoneinfo
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Reads NX, NY, NZ dimensions on each zone
subroutine read_size_CGNS(zone_number,char_filename,ngrid)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4),intent(in)   :: zone_number
  character(len=*), intent(in) :: char_filename
  integer(kind=4), intent(out) :: ngrid(3)

  integer(kind=4) :: m

  m = zone_number

  CGNS_filename = trim(char_filename)
  inquire(file=trim(CGNS_filename),exist=CGNS_file_exists)
    if (CGNS_file_exists .eqv. .false.) then
      write(*,'(A)') trim(CGNS_filename)
      stop ' O arquivo CGNS não existe!!! ' 
    endif

  !... open CGNS file
  call cg_open_f(trim(CGNS_filename),CG_MODE_READ,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  !... read dimensions, nx, ny and nz
  index_zone = m
  call cg_zone_read_f(index_grid,index_base,index_zone,zonename,isize,ier)

  ngrid(1) = isize(1,1)
  ngrid(2) = isize(2,1)
  if (cell_dim .gt. 2) then
    ngrid(3) = isize(3,1)
  else
    ngrid(3) = 1
  endif

  print*, ngrid(1)
  print*, ngrid(2)
  print*, ngrid(3)

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine read_size_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine create_file_CGNS(char_filename,aux_dim)

  use cgns
  use mod_CGNS
  implicit none
  
  character(len=*), intent(in) :: char_filename
  character(len=*), intent(in) :: aux_dim
  
  logical :: check
  
  CGNS_solnname = char_filename
  
  check = .false.
  
  inquire(file=CGNS_solnname,exist=CGNS_file_exists)
  if (CGNS_file_exists .eqv. .false.) then

    !open solution file
    call cg_open_f(trim(CGNS_solnname),CG_MODE_WRITE,index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    !create base in the recently opened solution file
    if (aux_dim .eq. '2D' .or. aux_dim .eq. '2d') then
      call cg_base_write_f(index_file,'Base',2,2,index_base,ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
      check = .true.
    endif

    if (aux_dim .eq. '3D' .or. aux_dim .eq. '3d') then
      call cg_base_write_f(index_file,'Base',3,3,index_base,ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
      check = .true.
    endif
    
    if (check .eqv. .false.) stop ' Tem cagada na chamada do create_file_CGNS ... '

    !close the file
    call cg_close_f(index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    write(*,'(A)') trim(char_filename) // ' already exists!'

  endif

  return

end subroutine create_file_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_link_CGNS(zone_number,char_filename,char_pathname)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  character(len=*), intent(in) :: char_filename
  character(len=*), intent(in) :: char_pathname
  integer(kind=4) :: m
  !integer(kind=4) :: nx1, nx2, nx3

  character(len=132) linkpath

  write(*,'(A)') '  -> Linking file ' // trim(char_filename) // ' to ' // trim(char_pathname)
  inquire(file=trim(char_pathname),exist=CGNS_file_exists)
    if (CGNS_file_exists .eqv. .false.) then
      print*, trim(char_pathname)
      stop ' O arquivo CGNS não existe!!!'
    endif
  
  !open solution file
  CGNS_solnname = trim(char_filename)
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  m = zone_number
    
    write(zonename,'(a4,i4.4)') 'Zone', m

    linkpath = '/Base/' // trim(zonename)

    call cg_gopath_f(index_file,trim(linkpath),ier)
      if (ier .ne. CG_OK) call cg_error_print_f

    linkpath = trim(linkpath) // '/GridCoordinates'

    call cg_link_write_f('GridCoordinates',trim(char_pathname),trim(linkpath),ier)
      if (ier .ne. CG_OK) call cg_error_print_f
      
    write(*,'(A,i0)') 'Linking zone ', m

  call cg_close_f(index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  return

end subroutine write_link_CGNS
!------------------------------------------------------------------------------------

!      |||||||||||||||     |||||||||||||       ||||||||||||  |||||||||||||||
!     |||||||||||||||||   |||||||||||||||      ||||||||||||  |||||||||||||||||
!     |||||||             ||||||    ||||||        ||||||     ||||||     |||||||
!     ||||||              ||||||     ||||||       ||||||     ||||||       ||||||
!     ||||||   ||||||||   ||||||    ||||||        ||||||     ||||||       ||||||
!     ||||||   |||||||||  |||||||||||||||         ||||||     ||||||       ||||||
!     ||||||      ||||||  |||||||||||||||         ||||||     ||||||       ||||||
!     |||||||     ||||||  ||||||     ||||||       ||||||     ||||||     |||||||
!     ||||||||||||||||||  ||||||      ||||||   ||||||||||||  |||||||||||||||||
!      ||||||||||||||||   ||||||       ||||||  ||||||||||||  |||||||||||||||

!------------------------------------------------------------------------------------
subroutine read_grid_2D_CGNS(zone_number,ngrid,char_filename,xcoord,ycoord)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)   :: zone_number
  integer(kind=4), intent(in)   :: ngrid(1:2)
  character(len=*), intent(in)  :: char_filename  
  real(kind=8), intent(out) :: xcoord(1:ngrid(1),1:ngrid(2))
  real(kind=8), intent(out) :: ycoord(1:ngrid(1),1:ngrid(2))
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2!, nx3
  integer(kind=4) :: ijk_min(2), ijk_max(2)

  CGNS_filename = trim(char_filename)

  m = zone_number

  call cg_open_f(trim(CGNS_filename),CG_MODE_READ,index_grid,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

    index_zone = m
    call cg_zone_read_f(index_grid,index_base,index_zone,zonename,isize,ier)

    nx1 = ngrid(1)
    nx2 = ngrid(2)

    ijk_min(1) = 1
    ijk_min(2) = 1

    ijk_max(1) = nx1
    ijk_max(2) = nx2

    call cg_coord_read_f(index_grid,index_base,index_zone,'CoordinateX',RealDouble,ijk_min(1:2),ijk_max(1:2),xcoord(:,:),ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
    call cg_coord_read_f(index_grid,index_base,index_zone,'CoordinateY',RealDouble,ijk_min(1:2),ijk_max(1:2),ycoord(:,:),ier)

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

end subroutine read_grid_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_partial_grid_2D_CGNS(zone_number,ngrid,xcoord,ycoord,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:2)
  real(kind=8), intent(in)     :: xcoord(1:ngrid(1),1:ngrid(2))
  real(kind=8), intent(in)     :: ycoord(1:ngrid(1),1:ngrid(2))
  character(len=*), intent(in) :: char_filename
  integer(kind=4)  :: m
  


  

  
  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

    call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
    if (ier .eq. CG_NODE_NOT_FOUND) then

      write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

      isize(1,1) = ngrid(1)
      isize(2,1) = ngrid(2)

      isize(1,2) = isize(1,1) - 1
      isize(2,2) = isize(2,1) - 1

      isize(1,3) = 0
      isize(2,3) = 0
      
      call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Structured,index_zone,ier)
        if (ier .ne. CG_OK) call cg_error_exit_f

      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateX',xcoord(:,:),index_coord,ier )
        if (ier .ne. CG_OK) call cg_error_exit_f
      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateY',ycoord(:,:),index_coord,ier )
        if (ier .ne. CG_OK) call cg_error_exit_f

!      !~~~~~~~~~~~~~
!      call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/ZoneGridConnectivity/', ier)
!        if (ier .eq. CG_NODE_NOT_FOUND) call cg_error_print_f
!      
!      npnts = 20
!      allocate(pnts(2,npnts))
!      pnts(1,1:20) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 /)
!      pnts(2,1:20) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)      
!      call  cg_hole_write_f(index_grid, index_base, index_zone, 'yabadabadoo', Vertex, PointList, 1, npnts, pnts, index_hole, ier)
!      
!      print*, ier

    else
    
      write(*,'(A,i0,A)') 'The grid coordinates for zone ', m, ' already exist! Skipping ...'

    endif

    call cg_close_f(index_grid,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

  return

end subroutine write_partial_grid_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine read_grid_3D_CGNS(zone_number,ngrid,char_filename,xcoord,ycoord,zcoord)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4),intent(in)   :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:3)  
  character(len=*), intent(in) :: char_filename
  real(kind=8), intent(out) :: xcoord(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  real(kind=8), intent(out) :: ycoord(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  real(kind=8), intent(out) :: zcoord(1:ngrid(1),1:ngrid(2),1:ngrid(3))  
  integer(kind=4) :: ijk_min(3), ijk_max(3)
  integer(kind=4) :: m

  CGNS_filename = trim(char_filename)

  call cg_open_f(trim(CGNS_filename),CG_MODE_READ,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number
  
    call cg_zone_read_f(index_grid,index_base,m,zonename,isize,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    ijk_min(1) = 1
    ijk_min(2) = 1
    ijk_min(3) = 1

    ijk_max(1) = ngrid(1)
    ijk_max(2) = ngrid(2)
    ijk_max(3) = ngrid(3)

    call cg_coord_read_f(index_grid,index_base,m,'CoordinateX',RealDouble,ijk_min(1:3),ijk_max(1:3),xcoord(:,:,:),ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
    call cg_coord_read_f(index_grid,index_base,m,'CoordinateY',RealDouble,ijk_min(1:3),ijk_max(1:3),ycoord(:,:,:),ier)
    call cg_coord_read_f(index_grid,index_base,m,'CoordinateZ',RealDouble,ijk_min(1:3),ijk_max(1:3),zcoord(:,:,:),ier)

  call cg_close_f(index_grid,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)
  
  return

end subroutine read_grid_3D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_partial_grid_3D_CGNS(zone_number,ngrid,xcoord,ycoord,zcoord,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:3)
  real(kind=8), intent(in)     :: xcoord(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  real(kind=8), intent(in)     :: ycoord(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  real(kind=8), intent(in)     :: zcoord(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  character(len=*), intent(in) :: char_filename
  integer(kind=4)  :: m
  

  
  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

    call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
      if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1,1) = ngrid(1)
    isize(2,1) = ngrid(2)
    isize(3,1) = ngrid(3)

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1
    isize(3,2) = isize(3,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0
    
    call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateX',xcoord(:,:,:),index_coord,ier )
      if (ier .ne. CG_OK) call cg_error_exit_f
    call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateY',ycoord(:,:,:),index_coord,ier )
    call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateZ',zcoord(:,:,:),index_coord,ier )

  else
  
    write(*,'(A,i0,A)') 'The grid coordinates for zone ', m, ' already exist! Skipping ...'

  endif

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

  return

end subroutine write_partial_grid_3D_CGNS
!------------------------------------------------------------------------------------

!      /|||||||||||||||||   |||||||||||||||||||   ||||||                 //////||||       //////
!     |||||||||||||||||||   |||||||||||||||||||   ||||||                //////|||||      //////
!     ||||||                ||||||       ||||||   ||||||               //////||||||     //////
!     ||||||                ||||||       ||||||   ||||||              ////// ||||||    //////
!     |||||||||||||||||||   ||||||       ||||||   ||||||             //////  ||||||   ////// 
!     |||||||||||||||||||   ||||||       ||||||   ||||||            //////   ||||||  //////
!                  ||||||   ||||||       ||||||   ||||||           //////    |||||| //////
!                  ||||||   ||||||       ||||||   ||||||          //////     ||||||//////
!     |||||||||||||||||||   |||||||||||||||||||   |||||||||||||  //////      |||||//////    
!     |||||||||||||||||/    |||||||||||||||||||   ||||||||||||| //////       ||||////// 
!
!------------------------------------------------------------------------------------
subroutine write_full_soln_2D_CGNS(zone_number,ngrid,solution,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:2)
  real(kind=8), intent(in)     :: solution(1:ngrid(1),1:ngrid(2))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2

  character(len=128) :: cdummy

  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)

  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1,1) = nx1
    isize(2,1) = nx2

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0

    call cg_zone_write_f(index_soln,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    index_zone = m !só funciona quando todas as zonas tiverem sido escritas..
    
    !call system('printf "\033[0;31m \n"')   
    !PRINT*, ' '
    !print*, 'MIGUE NA SUBROTINA "write_full_soln_2D_CGNS" !!!! '
    !PRINT*, ' '
    !call system('printf "\033[0m \n"') 

  endif

  call cg_zone_read_f(index_soln,index_base,index_zone,cdummy,isize,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/FlowSolution/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(index_soln,index_base,index_zone,'FlowSolution',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    index_flow = 1
    !call cg_sol_info_f(index_soln, index_base, index_zone, index_flow, solnname, location, ier)
    !  if (ier .ne. CG_OK) call cg_error_print_f
    !  if (trim(solnname) .ne. "Mean") stop ' A media nao se chama "Mean"... Tem alguma treta!!'

  endif

  call cg_field_write_f(index_soln,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name),solution(1:nx1,1:nx2),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_full_soln_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_full_soln_3D_CGNS(zone_number,ngrid,qdata,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none

  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:3)
  real(kind=8), intent(in)     :: qdata(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4)              :: m
  integer(kind=4)              :: nx1, nx2, nx3

  character(len=128) :: cdummy

  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

  write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)
  nx3 = ngrid(3)

  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(CGNS_solnname)

    isize(1,1) = nx1
    isize(2,1) = nx2
    isize(3,1) = nx3

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1
    isize(3,2) = isize(3,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

    call cg_zone_write_f(index_soln,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    index_zone = m !só funciona quando todas as zonas tiverem sido escritas..
    
    !call system('printf "\033[0;31m \n"')   
    !PRINT*, ' '
    !print*, 'MIGUE NA SUBROTINA "write_full_soln_3D_CGNS" !!!! '
    !PRINT*, ' '
    !call system('printf "\033[0m \n"') 

  endif

  call cg_zone_read_f(index_soln,index_base,index_zone,cdummy,isize,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  !call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/FlowSolution/', ier)
  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/teste/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    !call cg_sol_write_f(index_soln,index_base,index_zone,'FlowSolution',Vertex,index_flow,ier)
    call cg_sol_write_f(index_soln,index_base,index_zone,'teste',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    index_flow = 1
    call cg_sol_info_f(index_soln, index_base, index_zone, index_flow, solnname, location, ier)
      if (ier .ne. CG_OK) call cg_error_print_f
      if (trim(solnname) .ne. "teste") stop ' A solucao nao se chama "teste"... Tem alguma treta!!'

  endif

  call cg_field_write_f(index_soln,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name),qdata,index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    !if (ier .ne. CG_OK) call cg_error_print_f

  call cg_close_f(index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_full_soln_3D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine read_partial_soln_CGNS(zone_number,ijk_min,ijk_max,var_name,char_filename,output)

  use cgns
  use mod_CGNS
  !use mod_field
  implicit none
  
  integer(kind=4), intent(in) :: zone_number
  integer(kind=4), intent(in) :: ijk_min(3)
  integer(kind=4), intent(in) :: ijk_max(3)
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  real(kind=8), intent(out) :: output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2),ijk_min(3):ijk_max(3))

  integer(kind=4) :: m
  
  CGNS_solnname = trim(char_filename)
  
  write(*,'(1a1,A,$)') char(13), trim(CGNS_solnname)

  call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  index_base = 1
  call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
    if (ier .ne. CG_OK) call cg_error_exit_f

  allocate(isize(cell_dim,3))

  m = zone_number  

    index_zone = m
    call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
    
    index_flow = 1
    call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier)
        
!    call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/FlowSolution_t/', ier)
    call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
      if (ier .ne. CG_OK) call cg_error_exit_f  

    call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:3),ijk_max(1:3),&
                         output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2),ijk_min(3):ijk_max(3)),ier)
      if (ier .ne. CG_OK) call cg_error_exit_f  

  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

end subroutine read_partial_soln_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine read_partial_INT_soln_CGNS(zone_number,ijk_min,ijk_max,var_name,char_filename,output)

  use cgns
  use mod_CGNS
  !use mod_field
  implicit none
  
  integer(kind=4), intent(in) :: zone_number
  integer(kind=4), intent(in) :: ijk_min(3)
  integer(kind=4), intent(in) :: ijk_max(3)
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4), intent(out) :: output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2),ijk_min(3):ijk_max(3))

  integer(kind=4) :: m
  
  CGNS_solnname = trim(char_filename)
  
  write(*,'(1a1,A,$)') char(13), trim(CGNS_solnname)

  call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  index_base = 1
  call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
    if (ier .ne. CG_OK) call cg_error_exit_f

  allocate(isize(cell_dim,3))

  m = zone_number  

    index_zone = m
    call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
    
    index_flow = 1
    call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier)
        
!    call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/FlowSolution_t/', ier)
    call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
      if (ier .ne. CG_OK) call cg_error_exit_f  

    if (trim(var_name) .eq. "Iblank") then
      call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),Integer,ijk_min(1:3),ijk_max(1:3),&
                           output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2),ijk_min(3):ijk_max(3)),ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
    else
      stop ' Juquinha -- Iblank read soln'
    endif    

  call cg_close_f(index_soln,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

end subroutine read_partial_INT_soln_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
!
!             //////||||       //////||||  |||||||||||||          //////||||          //////||||       //////
!            //////|||||      //////|||||  |||||||||||||         //////|||||         //////|||||      //////
!           //////||||||     //////||||||  ||||||               //////||||||        //////||||||     //////
!          ////// ||||||    ////// ||||||  ||||||              ////// ||||||       ////// ||||||    //////
!         //////  ||||||   //////  ||||||  ||||||||||         //////  ||||||      //////  ||||||   ////// 
!        //////   ||||||  //////   ||||||  ||||||||||        //////   ||||||     //////   ||||||  //////
!       //////    |||||| //////    ||||||  ||||||           //////||||||||||    //////    |||||| //////
!      //////     ||||||//////     ||||||  ||||||          //////|||||||||||   //////     ||||||//////
!     //////      |||||//////      ||||||  |||||||||||||  //////      ||||||  //////      |||||//////
!    //////       ||||//////       ||||||  ||||||||||||| //////       |||||| //////       ||||//////
!
!------------------------------------------------------------------------------------
subroutine write_mean_soln_2D_CGNS(zone_number,ngrid,mean,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:2)
  real(kind=8), intent(in)     :: mean(1:ngrid(1),1:ngrid(2))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2!, nx3

  character(len=128) :: cdummy

  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname), CG_MODE_MODIFY, index_mean, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_mean, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)

  call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1,1) = nx1
    isize(2,1) = nx2

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0

    call cg_zone_write_f(index_mean, index_base, trim(zonename), isize, Structured, index_zone, ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    index_zone = m !só funciona quando todas as zonas tiverem sido escritas..
!    !FIXME FIXME
!    index_zone = output_zone
    call system('printf "\033[0;31m \n"')   
    PRINT*, ' '
    print*, 'MIGUE NA SUBROTINA "write_mean_soln_2D_CGNS" !!!! '
    PRINT*, ' '
    call system('printf "\033[0m \n"') 

  endif

  call cg_zone_read_f(index_mean,index_base,index_zone,cdummy,isize,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/Mean/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(index_mean,index_base,index_zone,'Mean',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    index_flow = 1
    !call cg_sol_info_f(index_mean, index_base, index_zone, index_flow, solnname, location, ier)
    !  if (ier .ne. CG_OK) call cg_error_print_f
    !  if (trim(solnname) .ne. "Mean") stop ' A media nao se chama "Mean"... Tem alguma treta!!'

  endif

  call cg_field_write_f(index_mean,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name),mean(1:nx1,1:nx2),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_mean_soln_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine read_mean_soln_2D_CGNS(zone_number,ijk_min,ijk_max,var_name,char_filename,mean)

  use cgns
  use mod_CGNS
  implicit none
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ijk_min(1:2)
  integer(kind=4), intent(in)  :: ijk_max(1:2)
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  real(kind=8), intent(out)    :: mean(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2))
  integer(kind=4) :: m

  CGNS_solnname = trim(char_filename)

  call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_mean,ier)
  if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_mean, index_base, basename, cell_dim, phys_dim, ier)
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
    if (ier .ne. CG_OK) call cg_error_exit_f

  allocate(isize(cell_dim,3))

  m = zone_number  

    index_zone = 1
    call cg_zone_read_f(index_mean,index_base,index_zone,zonename,isize,ier)
    
    call cg_gopath_f(index_mean, '/Base/'//trim(zonename), ier)
      if (ier .ne. CG_OK) call cg_error_exit_f 
          
    call cg_sol_info_f(index_mean,index_base,index_zone,index_flow,solnname,location,ier)
      !print*, 'solnname = ', trim(solnname)
      !if (trim(solnname) .ne. "Mean")  stop ' A solução nao se chama "Mean"... '
        
    call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/Mean/', ier)
      if (ier .ne. CG_OK) call cg_error_exit_f 

  write(*,*) trim(zonename)
  write(*,*) ijk_min(1:2)
  write(*,*) ijk_max(1:2)

  index_flow = 1
  !call cg_sol_info_f(index_mean, index_base, index_zone, index_flow, solnname, location, ier)
  !  if (ier .ne. CG_OK) call cg_error_print_f
  !  if (trim(solnname) .ne. "Mean") stop ' A media nao se chama "Mean"... Tem alguma treta!!'

  call cg_field_read_f(index_mean,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:2),ijk_max(1:2),mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

end subroutine read_mean_soln_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_mean_soln_3D_CGNS(zone_number,ngrid,mean,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:3)
  real(kind=8), intent(in)     :: mean(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2, nx3

  character(len=128) :: cdummy

  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_mean, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)
  nx3 = ngrid(3)

  call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1,1) = nx1
    isize(2,1) = nx2
    isize(3,1) = nx3

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1
    isize(3,2) = isize(3,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

    call cg_zone_write_f(index_mean,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    index_zone = m !só funciona quando todas as zonas tiverem sido escritas..

    !call system('printf "\033[0;31m \n"')   
    !PRINT*, ' '
    !print*, 'MIGUE NA SUBROTINA "write_mean_soln_3D_CGNS" !!!! '
    !PRINT*, ' '
    !call system('printf "\033[0m \n"') 

  endif

  call cg_zone_read_f(index_mean,index_base,index_zone,cdummy,isize,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/Mean/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(index_mean,index_base,index_zone,'Mean',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    index_flow = 1
    !call cg_sol_info_f(index_mean, index_base, index_zone, index_flow, solnname, location, ier)
    !  if (ier .ne. CG_OK) call cg_error_print_f
    !  if (trim(solnname) .ne. "Mean") stop ' A media nao se chama "Mean"... Tem alguma treta!!'

  endif

  call cg_field_write_f(index_mean,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name),mean(1:nx1,1:nx2,1:nx3),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_mean_soln_3D_CGNS
!------------------------------------------------------------------------------------






!------------------------------------------------------------------------------------
subroutine write_fourier_modes_2D_CGNS(zone_number,ngrid,ibin,dft_mode,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:2)
  integer(kind=4), intent(in)  :: ibin
  complex(kind=8), intent(in)  :: dft_mode(1:ngrid(1),1:ngrid(2))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2


  
  character(len=6)   :: caux

  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_mode,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_mode, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number
  write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)

  call cg_gopath_f(index_mode, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4)') '  -> Creating zone ', m

    isize(1,1) = nx1
    isize(2,1) = nx2

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0

    call cg_zone_write_f(index_mode,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    !index_zone = m só funciona quando todas as zonas tiverem sido escritas..
    !FIXME FIXME
!    index_zone = 1
    index_zone = m
    call system('printf "\033[0;31m \n"')   
    PRINT*, ' '
    PRINT*, index_zone, ibin
    PRINT*, 'MIGUE NA SUBROTINA "write_fourier_modes_2D_CGNS" !!!! '
    PRINT*, ' '
    call system('printf "\033[0m \n"') 

  endif
  
!  call cg_zone_read_f(index_mode,index_base,index_zone,cdummy,isize,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

  !~~~ Parte real
  write(caux,'(a4,i2.2)') 'Real', ibin
  call cg_sol_write_f(index_mode,index_base,index_zone,caux,Vertex,index_flow,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_field_write_f(index_mode,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name)//'_R', dble(dft_mode(1:nx1,1:nx2)),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  !~~~

!  !~~~ Parte imaginária
!  write(caux,'(a4,i2.2)') 'Imag', ibin
!  call cg_sol_write_f(index_mode,index_base,index_zone,caux,Vertex,index_flow,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!  call cg_field_write_f(index_mode,index_base,index_zone,index_flow, &
!                        RealDouble,trim(var_name),dimag(dft_mode(1:nx1,1:nx2)),index_field,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f
!  !~~~
  
!  !~~~ Absoluto
!  write(caux,'(a4,i2.2)') 'Magn', ibin
!  call cg_sol_write_f(index_mode,index_base,index_zone,caux,Vertex,index_flow,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f
!
!  call cg_field_write_f(index_mode,index_base,index_zone,index_flow, &
!                        RealDouble,trim(var_name)//'_A',cdabs(dft_mode(1:nx1,1:nx2)),index_field,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f
!  !~~~  
  
  call cg_close_f(index_mode,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_fourier_modes_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_fourier_modes_3D_CGNS_dummy(zone_number,ngrid,ibin,dft_mode,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:3)
  integer(kind=4), intent(in)  :: ibin
  complex(kind=8), intent(in)  :: dft_mode(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2, nx3


  
  character(len=6)   :: caux

  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_mode,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_mode, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number
  write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)
  nx3 = ngrid(3)

  call cg_gopath_f(index_mode, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4)') '  -> Creating zone ', m

    isize(1,1) = nx1
    isize(2,1) = nx2
    isize(3,1) = nx3

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1
    isize(3,2) = isize(3,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

    call cg_zone_write_f(index_mode,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    !index_zone = m só funciona quando todas as zonas tiverem sido escritas..
    !FIXME FIXME
!    index_zone = 1
    index_zone = m
    call system('printf "\033[0;31m \n"')   
!    PRINT*, ' '
    PRINT*, index_zone, ibin
    PRINT*, 'MIGUE NA SUBROTINA "write_fourier_modes_3D_CGNS" !!!! '
!    PRINT*, ' '
    call system('printf "\033[0m \n"')     

  endif

!  call cg_zone_read_f(index_mode,index_base,index_zone,cdummy,isize,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

  !~~~ Parte real
  write(caux,'(a4,i2.2)') 'Real', ibin
  call cg_sol_write_f(index_mode,index_base,index_zone,caux,Vertex,index_flow,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_field_write_f(index_mode,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name), dble(dft_mode(1:nx1,1:nx2,1:nx3)),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  !~~~

!  !~~~ Parte imaginária
!  write(caux,'(a4,i2.2)') 'Imag', ibin
!  call cg_sol_write_f(index_mode,index_base,index_zone,caux,Vertex,index_flow,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!  call cg_field_write_f(index_mode,index_base,index_zone,index_flow, &
!                        RealDouble,trim(var_name),dimag(dft_mode(1:nx1,1:nx2,1:nx3)),index_field,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f
!  !~~~
  
  !~~~ Absoluto
  write(caux,'(a4,i2.2)') 'Magn', ibin
  call cg_sol_write_f(index_mode,index_base,index_zone,caux,Vertex,index_flow,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_field_write_f(index_mode,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name),cdabs(dft_mode(1:nx1,1:nx2,1:nx3)),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  !~~~
  
  call cg_close_f(index_mode,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_fourier_modes_3D_CGNS_dummy
!------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------
subroutine write_fourier_modes_3D_CGNS(zone_number,ngrid,ibin,dft_mode,var_name,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:3)
  integer(kind=4), intent(in)  :: ibin
  complex(kind=8), intent(in)  :: dft_mode(1:ngrid(1),1:ngrid(2),1:ngrid(3))
  character(len=*), intent(in) :: var_name
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2, nx3

  character(len=128) :: cdummy
  


  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_mean, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)
  nx3 = ngrid(3)

  call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1,1) = nx1
    isize(2,1) = nx2
    isize(3,1) = nx3

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1
    isize(3,2) = isize(3,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

    call cg_zone_write_f(index_mean,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    !index_zone = m só funciona quando todas as zonas tiverem sido escritas..
    !FIXME FIXME
    index_zone = m
    call system('printf "\033[0;31m \n"')   
!    PRINT*, ' '
    PRINT*, index_zone, ibin
    PRINT*, 'MIGUE NA SUBROTINA "write_fourier_modes_3D_CGNS" !!!! '
!    PRINT*, ' '
    call system('printf "\033[0m \n"') 

  endif

  call cg_zone_read_f(index_mean,index_base,index_zone,cdummy,isize,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(index_mean, '/Base/'//trim(zonename)//'/teste/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(index_mean,index_base,index_zone,'teste',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    index_flow = 1
    !call cg_sol_info_f(index_mean, index_base, index_zone, index_flow, solnname, location, ier)
    !  if (ier .ne. CG_OK) call cg_error_print_f
    !  if (trim(solnname) .ne. "Mean") stop ' A media nao se chama "Mean"... Tem alguma treta!!'

  endif

  call cg_field_write_f(index_mean,index_base,index_zone,index_flow, &
                        RealDouble,trim(var_name)//'_R',real(dft_mode(1:nx1,1:nx2,1:nx3),8),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_mean,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_fourier_modes_3D_CGNS
!------------------------------------------------------------------------------------













!  !------------------------------------------------------------------------------------
!  subroutine write_link_IBLANK_CGNS(zone_number,char_filename,char_pathname)

!    use cgns
!    use mod_CGNS
!    implicit none
!    
!    integer(kind=4), intent(in)  :: zone_number
!    character(len=*), intent(in) :: char_filename
!    character(len=*), intent(in) :: char_pathname
!    integer(kind=4) :: m
!    !integer(kind=4) :: nx1, nx2, nx3

!    character(len=132) linkpath

!    write(*,'(A)') '  -> Linking file ' // trim(char_filename) // ' to ' // trim(char_pathname)

!    !open solution file
!    CGNS_solnname = trim(char_filename)
!    call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_file,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    m = zone_number
!      
!      write(zonename,'(a4,i4.4)') 'Zone', m

!      !FIXME
!      linkpath = '/Base/' // trim(zonename) // '/teste'

!      call cg_gopath_f(index_file,trim(linkpath),ier)
!        if (ier .ne. CG_OK) call cg_error_print_f

!      linkpath = trim(linkpath) // '/Iblank'

!      call cg_link_write_f('dummy',trim(char_pathname),trim(linkpath),ier)
!        if (ier .ne. CG_OK) call cg_error_print_f
!        
!      write(*,'(A,i0)') 'Linking zone ', m

!    call cg_close_f(index_file,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    return

!  end subroutine write_link_IBLANK_CGNS
!  !------------------------------------------------------------------------------------



!  !------------------------------------------------------------------------------------
!  subroutine write_grid_2D_CGNS(zone_number,ngrid,xcoord,ycoord,char_filename)

!    use cgns
!    use mod_CGNS
!    !use mod_field
!    implicit none
!    
!    integer(kind=4), intent(in)  :: zone_number
!    integer(kind=4), intent(in)  :: ngrid(1:2)
!    real(kind=8), intent(in)     :: xcoord(1:ngrid(1),1:ngrid(2))
!    real(kind=8), intent(in)     :: ycoord(1:ngrid(1),1:ngrid(2))
!    character(len=*), intent(in) :: char_filename
!    integer(kind=4)  :: m
!    integer(kind=4)  :: nx1, nx2
!    
!    character(len=128) :: cdummy
!    
!    m = zone_number

!    !... open CGNS file
!    CGNS_filename = trim(char_filename)
!    call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    index_base = 1
!    call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
!      if (ier .ne. CG_OK) call cg_error_print_f
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

!    allocate(isize(cell_dim,3))
!    
!    write(zonename,'(a4,i4.4)') 'Zone', m

!    call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
!    if (ier .eq. CG_NODE_NOT_FOUND) then

!      write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)
!    
!      nx1 = ngrid(1)
!      nx2 = ngrid(2)

!      isize(1,1) = nx1
!      isize(2,1) = nx2

!      isize(1,2) = isize(1,1) - 1
!      isize(2,2) = isize(2,1) - 1

!      isize(1,3) = 0
!      isize(2,3) = 0
!      
!      call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Structured,index_zone,ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateX',xcoord(1:nx1,1:nx2),index_coord,ier )
!        if (ier .ne. CG_OK) call cg_error_exit_f
!      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateY',ycoord(1:nx1,1:nx2),index_coord,ier )

!    else
!    
!      write(*,'(A,i0,A)') 'The grid coordinates for zone ', m, ' already exist! Skipping ...'

!    endif

!    call cg_close_f(index_grid,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!    return

!  end subroutine write_grid_2D_CGNS
!  !------------------------------------------------------------------------------------



!  !------------------------------------------------------------------------------------
!  subroutine write_grid_3D_CGNS(zone_number,ngrid,char_filename)

!    use cgns
!    use mod_CGNS
!    use mod_field
!    implicit none
!    
!    integer(kind=4), intent(in)  :: zone_number
!    integer(kind=4), intent(in)  :: ngrid(1:3)
!    character(len=*), intent(in) :: char_filename
!    integer(kind=4)  :: m
!    integer(kind=4)  :: nx1, nx2, nx3
!    
!    character(len=128) :: cdummy
!    
!    !... open CGNS file
!    CGNS_filename = trim(char_filename)
!    call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    index_base = 1
!    call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
!      if (ier .ne. CG_OK) call cg_error_print_f
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

!    allocate(isize(cell_dim,3))

!    m = zone_number

!      write(zonename,'(a4,i4.4)') 'Zone', m

!      call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
!        if (ier .eq. CG_NODE_NOT_FOUND) then

!      write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)
!    
!      nx1 = ngrid(1)
!      nx2 = ngrid(2)
!      nx3 = ngrid(3)

!      isize(1,1) = nx1
!      isize(2,1) = nx2
!      isize(3,1) = nx3

!      isize(1,2) = isize(1,1) - 1
!      isize(2,2) = isize(2,1) - 1
!      isize(3,2) = isize(3,1) - 1

!      isize(1,3) = 0
!      isize(2,3) = 0
!      isize(3,3) = 0
!      
!      call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Structured,index_zone,ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateX',zone(m)%x(:,:,:),index_coord,ier )
!        if (ier .ne. CG_OK) call cg_error_exit_f
!      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateY',zone(m)%y(:,:,:),index_coord,ier )
!      call cg_coord_write_f(index_grid,index_base,index_zone,RealDouble,'CoordinateZ',zone(m)%z(:,:,:),index_coord,ier )

!    else
!    
!      write(*,'(A,i0,A)') 'The grid coordinates for zone ', m, ' already exist! Skipping ...'

!    endif

!    call cg_close_f(index_grid,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!    return

!  end subroutine write_grid_3D_CGNS
!  !------------------------------------------------------------------------------------



!  !------------------------------------------------------------------------------------
!  subroutine read_full_soln_2D_CGNS(step,niter,zone_number,var_name,char_filename,output)

!    use cgns
!    use mod_CGNS
!    use mod_field
!    implicit none
!    
!    integer(kind=4), intent(in)  :: step
!    integer(kind=4), intent(in)  :: niter
!    integer(kind=4), intent(in)  :: zone_number
!    character(len=*), intent(in) :: var_name
!    character(len=*), intent(in) :: char_filename
!    real(kind=8), intent(out)    :: output(1:zone(zone_number)%nx,1:zone(zone_number)%ny,1:zone(zone_number)%nz)

!    integer(kind=4) :: m
!    integer(kind=4) :: nx1, nx2, nx3
!    integer(kind=4) :: ijk_min(3), ijk_max(3)
!    
!    CGNS_solnname = trim(char_filename)
!    
!    write(*,'(1a1,A,$)') char(13), trim(CGNS_solnname)

!    call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
!      if (ier .ne. CG_OK) call cg_error_print_f

!    index_base = 1
!    call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    allocate(isize(cell_dim,3))

!    m = zone_number  

!      index_zone = m
!      call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
!      
!      index_flow = 1
!      call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier)
!          
!      call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f  

!      nx1 = zone(m)%nx
!      nx2 = zone(m)%ny
!      nx3 = zone(m)%nz

!      if (nx3 .gt. 1) stop '  Resultado esta tridimensional... Deveria ser bidimensional!!'

!      ijk_min(1) = 1
!      ijk_min(2) = 1
!      !ijk_min(3) = 1

!      ijk_max(1) = nx1
!      ijk_max(2) = nx2
!      !ijk_max(3) = nx3

!      call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:2),ijk_max(1:2),output,ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!    call cg_close_f(index_soln,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!  end subroutine read_full_soln_2D_CGNS
!  !------------------------------------------------------------------------------------

!  !------------------------------------------------------------------------------------
!  subroutine read_slice_soln_3D_CGNS(niter,zone_number,zslice,var_name,char_filename,output)

!    use cgns
!    use mod_CGNS
!    use mod_field
!    implicit none
!    
!    integer(kind=4), intent(in) :: niter
!    integer(kind=4), intent(in) :: zone_number
!    integer(kind=4), intent(in) :: zslice
!    character(len=*), intent(in) :: var_name
!    character(len=*), intent(in) :: char_filename
!    real(kind=8), intent(out) :: output(1:zone(zone_number)%nx,1:zone(zone_number)%ny,zslice:zslice)

!    integer(kind=4) :: m
!    integer(kind=4) :: nx1, nx2, nx3
!    integer(kind=4) :: ijk_min(3), ijk_max(3)

!    !FIXME: essa subrotina tem que ser melhorada
!    !FIXME: essa subrotina tem que ser melhorada
!    !FIXME: essa subrotina tem que ser melhorada
!    !FIXME: essa subrotina tem que ser melhorada
!    !FIXME: essa subrotina tem que ser melhorada

!    CGNS_solnname = trim(char_filename)
!    write(*,'(1a1,2A,$)') char(13), '   |~~ ', trim(CGNS_solnname)

!    call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    index_base = 1
!    call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
!      if (ier .ne. CG_OK) call cg_error_print_f
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

!    allocate(isize(cell_dim,3))

!    ! call cg_gopath_f(index_soln, '/Base/TimeIterValues', ier)
!    ! if (ier .ne. CG_OK) call cg_error_exit_f
!    ! ! Read time instant
!    ! call cg_array_read_f(1,tmp,ier)
!    ! if (ier .ne. CG_OK) call cg_error_exit_f

!    m = zone_number

!      index_zone = m
!      call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!      nx1 = zone(m)%nx
!      nx2 = zone(m)%ny
!      nx3 = zone(m)%nz

!      if (nx3 .eq. 1) stop '  Resultado esta bidimensional... Deveria ser tridimensional!!'

!      ijk_min(1) = 1
!      ijk_min(2) = 1
!      ijk_min(3) = zslice

!      ijk_max(1) = nx1
!      ijk_max(2) = nx2
!      ijk_max(3) = zslice

!  !    write(*,*) trim(zonename)
!  !    write(*,*) ' ijk_min(1:3) = ', ijk_min(1:3)
!  !    write(*,*) ' ijk_max(1:3) = ', ijk_max(1:3)

!      !FIXME
!      !index_flow = 1
!      !call cg_sol_info_f(index_soln, index_base, index_zone, index_flow, solnname, location, ier)
!      !  if (ier .ne. CG_OK) call cg_error_exit_f

!      call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:3),ijk_max(1:3),&
!                           output(1:ijk_max(1),1:ijk_max(2),ijk_min(3):ijk_max(3)),ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!    call cg_close_f(index_soln,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!  end subroutine read_slice_soln_3D_CGNS
!  !------------------------------------------------------------------------------------



!  !------------------------------------------------------------------------------------
!  subroutine read_full_soln_3D_CGNS(zone_number,ijk_min,ijk_max,var_name,char_filename,output)

!    use cgns
!    use mod_CGNS
!    !use mod_field
!    implicit none
!    
!    integer(kind=4), intent(in) :: zone_number
!    character(len=*), intent(in) :: var_name
!    character(len=*), intent(in) :: char_filename
!    real(kind=8), intent(out) :: output(1:zone(zone_number)%nx,1:zone(zone_number)%ny,1:zone(zone_number)%nz)

!    integer(kind=4) :: m
!    integer(kind=4) :: nx1, nx2, nx3
!    integer(kind=4) :: ijk_min(3), ijk_max(3)
!    
!    CGNS_solnname = trim(char_filename)
!    
!    write(*,'(1a1,A,$)') char(13), trim(CGNS_solnname)

!    call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!    index_base = 1
!    call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    allocate(isize(cell_dim,3))

!    m = zone_number  

!      index_zone = m
!      call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
!      
!      index_flow = 1
!      call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier)
!          
!  !    call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/FlowSolution_t/', ier)
!      call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f   

!      nx1 = zone(m)%nx
!      nx2 = zone(m)%ny
!      nx3 = zone(m)%nz

!      if (nx3 .eq. 1) stop '  Resultado esta bidimensional... Deveria ser tridimensional!!'

!      ijk_min(1) = 1
!      ijk_min(2) = 1
!      ijk_min(3) = 1

!      ijk_max(1) = nx1
!      ijk_max(2) = nx2
!      ijk_max(3) = nx3

!      call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:3),ijk_max(1:3),output,ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!    call cg_close_f(index_soln,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!  end subroutine read_full_soln_3D_CGNS
!  !------------------------------------------------------------------------------------

!  !------------------------------------------------------------------------------------
!  subroutine read_partial_soln_3D_CGNS(zone_number,ijk_min,ijk_max,var_name,char_filename,output)

!    use cgns
!    use mod_CGNS
!    !use mod_field
!    implicit none
!    
!    integer(kind=4), intent(in) :: zone_number
!    integer(kind=4), intent(in) :: ijk_min(3)
!    integer(kind=4), intent(in) :: ijk_max(3)
!    character(len=*), intent(in) :: var_name
!    character(len=*), intent(in) :: char_filename
!    real(kind=8), intent(out) :: output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2),ijk_min(3):ijk_max(3))

!    integer(kind=4) :: m
!    
!    CGNS_solnname = trim(char_filename)
!    
!    write(*,'(1a1,A,$)') char(13), trim(CGNS_solnname)

!    call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
!      if (ier .ne. CG_OK) call cg_error_print_f

!    index_base = 1
!    call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    allocate(isize(cell_dim,3))

!    m = zone_number  

!      index_zone = m
!      call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
!      
!      index_flow = 1
!      call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier)
!          
!  !    call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/FlowSolution_t/', ier)
!      call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f  

!      call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:3),ijk_max(1:3),&
!                           output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2),ijk_min(3):ijk_max(3)),ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!    call cg_close_f(index_soln,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!  end subroutine read_partial_soln_3D_CGNS
!  !------------------------------------------------------------------------------------

!  !------------------------------------------------------------------------------------
!  subroutine read_partial_soln_2D_CGNS(zone_number,ijk_min,ijk_max,var_name,char_filename,output)

!    use cgns
!    use mod_CGNS
!    implicit none
!    integer(kind=4), intent(in) :: zone_number
!    integer(kind=4), intent(in) :: ijk_min(1:2)
!    integer(kind=4), intent(in) :: ijk_max(1:2)
!    character(len=*), intent(in) :: var_name
!    character(len=*), intent(in) :: char_filename
!    real(kind=8), intent(out) :: output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2))

!    integer(kind=4) :: m
!    integer(kind=4) :: nx1, nx2, nx3
!    
!    CGNS_solnname = trim(char_filename)
!    
!    write(*,'(1a1,A,$)') char(13), trim(CGNS_solnname)

!    call cg_open_f(trim(CGNS_solnname),CG_MODE_READ,index_soln,ier)
!      if (ier .ne. CG_OK) call cg_error_print_f

!    index_base = 1
!    call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
!      if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '
!      if (ier .ne. CG_OK) call cg_error_exit_f

!    allocate(isize(cell_dim,3))

!    m = zone_number  

!      index_zone = m
!      call cg_zone_read_f(index_soln,index_base,index_zone,zonename,isize,ier)
!      
!      index_flow = 1
!      call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier)

!      call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f  

!      call cg_field_read_f(index_soln,index_base,index_zone,index_flow,trim(var_name),RealDouble,ijk_min(1:2),ijk_max(1:2), &
!                                              output(ijk_min(1):ijk_max(1),ijk_min(2):ijk_max(2)),ier)
!        if (ier .ne. CG_OK) call cg_error_exit_f

!    call cg_close_f(index_soln,ier)
!    if (ier .ne. CG_OK) call cg_error_exit_f

!    DEALLOCATE(isize)

!  end subroutine read_partial_soln_2D_CGNS
!  !------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------
subroutine write_iblank_2D_CGNS(zone_number,ngrid,iblank,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: ngrid(1:2)
  integer(kind=4), intent(in)  :: iblank(1:ngrid(1),1:ngrid(2))
  character(len=*), intent(in) :: char_filename
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2

  character(len=128) :: cdummy


  CGNS_solnname = trim(char_filename)

  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_MODIFY,index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_soln, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  m = zone_number

  write(zonename,'(a4,i4.4)') 'Zone', m

  nx1 = ngrid(1)
  nx2 = ngrid(2)

  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/', ier)
  
  index_zone = 1

  call cg_zone_read_f(index_soln,index_base,index_zone,cdummy,isize,ier)
  if (ier .ne. CG_OK) call cg_error_print_f

  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1,1) = nx1
    isize(2,1) = nx2

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0

    call cg_zone_write_f(index_soln,index_base,trim(zonename),isize,Structured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      
  else
  
    !a reordenação das zonas pode complicar um cado.. Tem que pensar no jeito mais robusto possível!!
    index_zone = m !só funciona quando todas as zonas tiverem sido escritas..
    !FIXME FIXME
    !index_zone = 1
    
    !index_flow = 1
    !index_field = 1
    !call cg_field_info_f(index_soln,index_base, index_zone, index_flow, index_field, idummy, fieldname, ier)
    !  print*, trim(fieldname)
    !  if (ier .ne. CG_OK) call cg_error_exit_f 
    
    index_flow = 1 
    call cg_sol_info_f(index_soln,index_base,index_zone,index_flow,solnname,location,ier) 
      print*, trim(solnname)
      if (ier .ne. CG_OK) call cg_error_exit_f 
    
    call system('printf "\033[0;31m \n"')   
    PRINT*, ' '
    print*, 'MIGUE NA SUBROTINA "write_full_soln_2D_CGNS" !!!! '
    PRINT*, ' '
    call system('printf "\033[0m \n"') 

  endif

  call cg_zone_read_f(index_soln,index_base,index_zone,cdummy,isize,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(index_soln, '/Base/'//trim(zonename)//'/'//trim(solnname)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(index_soln,index_base,index_zone,trim(solnname),Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    index_flow = 1
    !call cg_sol_info_f(index_soln, index_base, index_zone, index_flow, solnname, location, ier)
    !  if (ier .ne. CG_OK) call cg_error_print_f
    !  if (trim(solnname) .ne. "Mean") stop ' A media nao se chama "Mean"... Tem alguma treta!!'

  endif

  call cg_field_write_f(index_soln,index_base,index_zone,index_flow, &
                           Integer,'iblk',iblank(1:nx1,1:nx2),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_soln,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_iblank_2D_CGNS
!------------------------------------------------------------------------------------
