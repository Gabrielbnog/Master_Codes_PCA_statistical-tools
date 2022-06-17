!##################################################################################
subroutine extract_grid

  !extrair a malha de superficie do aerofolio
  !eventualmente pode-se mudar para outra posição caso seja necessário...

  use mod_CGNS
  use mod_field
  use mod_pod_modes, only : working_zone
  implicit none
  integer(kind=4) i, m
  
  m = working_zone
  
  write(*,'(A)') ' Exporting grid line ... '

  do i = imin(m),imax(m)-1
    
    write(80,'(3f20.14)') zone(m)%x(i,1,1), zone(m)%y(i,1,1), 0.0d0
    
  enddo

  return

end subroutine extract_grid
!##################################################################################

!##################################################################################
!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine export_soln(zone_number, idx_beg, idx_skip, idx_end, export_grid)

  use mod_field
  use mod_CGNS
  implicit none
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: idx_beg, idx_skip, idx_end
  logical, intent(in)          :: export_grid
  ! character(len=*), intent(in) :: new_grid
  integer(kind=4) :: m, t
  integer(kind=4) :: i, file_snap
  
  integer(kind=4) :: nx1, nx2, nx3
  
  m = zone_number
  
  nx1 = imax(m)-imin(m)+1
  nx2 = jmax(m)
  nx3 = kmax(m)
  
  t = 0
  do i = idx_beg,idx_end,idx_skip

    t = t + 1

    file_snap = (idx_beg-1)*idxr*idx_skip + (t-1)*idxr*idx_skip + idxi

    write(CGNS_solnname,'(A,A,i6.6,A5)') trim(output_path), 'soln', file_snap, '.cgns'
    
    if (data_2D .eqv. .true.) then
      call create_file_CGNS(trim(CGNS_solnname),'2D')
      call write_full_soln_2D_CGNS(m,[nx1,nx2],zone(m)%q(1:nx1,1:nx2,1,t),trim(working_var),trim(CGNS_solnname))
      call write_link_CGNS(m,trim(CGNS_solnname),trim(output_path)//'/grid_mod_2D.cgns')
    endif

    if (data_3D .eqv. .true.) then
      call create_file_CGNS(trim(CGNS_solnname),'3D')
      call write_full_soln_3D_CGNS(m,[nx1,nx2,nx3],zone(m)%q(1:nx1,1:nx2,1:nx3,t),trim(working_var),trim(CGNS_solnname))
      call write_link_CGNS(m,trim(CGNS_solnname),trim(output_path)//'/grid_mod.cgns')
    endif
    
  enddo
  write(*,*) ''
  
  if (export_grid .eqv. .false.) return
  
  !TODO eventuar adicionar a possibilidade de criar um novo arquivo de malha só para o chunk de interesse
  
  return

end subroutine
!####################################################################################

!####################################################################################
!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine export_ASCII_grid(zone_number)

  use mod_field
  use mod_CGNS
  use mod_pod_modes, only : istride, jstride, kstride
  implicit none
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4) :: i, j, k, m
  
  character(len=5) :: char_zone
  
  integer(kind=4) :: nx1, nx2, nx3
  
  m = zone_number
  write(char_zone,'(a1,i4.4)') 'z', m
  
  nx1 = 0
  do i = imin(m),imax(m),istride
    nx1 = nx1 + 1
  enddo

  nx2 = 0
  do i = 1,jmax(m),jstride
    nx2 = nx2 + 1
  enddo
  
  nx3 = 0
  do i = 1,kmax(m),kstride
    nx3 = nx3 + 1
  enddo

  !open(666,file='ascii_grid_'//char_zone//'.dat')
  open(666,file='grid_'//char_zone//'.dat')
!  write(666,'(A)') 'FILETYPE=GRID'  
  
  if (nx3 .eq. 1) then
  
!    write(666,'(A)') 'VARIABLES="X","Y"'
!    write(666,'(3(A,i0),A)') 'ZONE T="'//char_zone//'",I=',nx1,',J=',nx2,',F=POINT'     
    k = 1
    do j = 1,jmax(m),jstride
      do i = imin(m),imax(m),istride
        write(666,'(2e23.15)') zone(m)%x(i,j,k), zone(m)%y(i,j,k)
      enddo
    enddo
    
  else
  
    write(666,'(A)') 'VARIABLES="X","Y","Z"'
    write(666,'(3(A,i0),A)') 'ZONE T="'//char_zone//'",I=',nx1,',J=',nx2,',K=',nx3,',F=BLOCK'
    do k = 1,kmax(m)
      do j = 1,jmax(m),jstride
        do i = imin(m),imax(m),istride
          write(666,'(e23.15)') zone(m)%x(i,j,k)
        enddo
      enddo
    enddo
    
    do k = 1,kmax(m)
      do j = 1,jmax(m),jstride
        do i = imin(m),imax(m),istride
          write(666,'(e23.15)') zone(m)%y(i,j,k)
        enddo
      enddo
    enddo
    
    do k = 1,kmax(m)
      do j = 1,jmax(m),jstride
        do i = imin(m),imax(m),istride
          write(666,'(e23.15)') zone(m)%z(i,j,k)
        enddo
      enddo
    enddo        
    
  endif
  
  close(666)
  
  call system('preplot grid_'//char_zone//'.dat')
  
  return

end subroutine
!####################################################################################

!####################################################################################
!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine export_ASCII_soln(zone_number, idx_beg, idx_skip, idx_end)

  use mod_field
  use mod_CGNS
  use mod_pod_modes, only : istride, jstride
  implicit none
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: idx_beg, idx_skip, idx_end
  integer(kind=4) :: i, j, k, m, n, t, file_snap
  
  character(len=5) :: char_zone
  character(len=64) :: filename
  
  integer(kind=4) :: nx1, nx2, nx3
  
  m = zone_number
  write(char_zone,'(a1,i4.4)') 'z', m
  
  nx1 = 0
  do i = imin(m),imax(m),istride
    nx1 = nx1 + 1
  enddo

  nx2 = 0
  do i = 1,jmax(m),jstride
    nx2 = nx2 + 1
  enddo
  
  nx3 = 1
 
  t = 0
  do n = idx_beg,idx_end,idx_skip

    t = t + 1
    print*, t

    file_snap = (idx_beg-1)*idxr*idx_skip + (t-1)*idxr*idx_skip + idxi

    !if (trim(working_var) .eq. 'VelocityX') write(filename,'(a6,i6.6,a6)') 'ascii_', file_snap, '_U.dat'
    !if (trim(working_var) .eq. 'VelocityY') write(filename,'(a6,i6.6,a6)') 'ascii_', file_snap, '_V.dat'    

    if (trim(working_var) .eq. 'VelocityX') write(filename,'(i6.6,a6)') file_snap, '_U.dat'
    if (trim(working_var) .eq. 'VelocityY') write(filename,'(i6.6,a6)') file_snap, '_V.dat' 
    open(666,file=trim(filename))
!    write(666,'(A)') 'FILETYPE=SOLUTION'
!    if (trim(working_var) .eq. 'VelocityX') write(666,'(A)') 'VARIABLES="U"'
!    if (trim(working_var) .eq. 'VelocityY') write(666,'(A)') 'VARIABLES="V"'
!    write(666,'(3(A,i0),A,f10.3)') 'ZONE T="'//char_zone//'",I=',nx1,',J=',nx2,',K=',nx3,',F=POINT,SOLUTIONTIME=',dble(t)
    k = 1
    do j = 1,jmax(m),jstride
      do i = imin(m),imax(m),istride
        write(666,'(e23.15)') zone(m)%q(i,j,k,t)
      enddo
    enddo
    close(666)
    
  enddo
  write(*,*) ''
  
  return

end subroutine
!####################################################################################

!####################################################################################
subroutine initial_condition

  use mod_CGNS
  use mod_field
  implicit none
  integer(kind=4) :: i, j, k, m
  integer(kind=4) :: nx1, nx2, nx3
  
  character(len=5) :: char_zone
  
  m = output_zone
  
    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    nx3 = kmax(m)    

    write(char_zone,'(a1,i4.4)') 'z', output_zone
    if (trim(working_var) .eq. 'VelocityX') then
      open(99,file=trim(output_path)//'initial_'//char_zone//'_U.dat')
      write(99,'(A)') 'VARIABLES="X","Y","U"'      
    endif
    if (trim(working_var) .eq. 'VelocityY') then
      open(99,file=trim(output_path)//'initial_'//char_zone//'_V.dat')    
      write(99,'(A)') 'VARIABLES="X","Y","V"'
    endif

    write(99,'(3(A,i0),A)') 'ZONE T="'//char_zone//'",I=',nx1,',J=',jmax(m),',K=',kmax(m),',F=POINT'                          

    do k = 1,kmax(m)
      do j = 1,jmax(m)
        do i = imin(m),imax(m)
          write(99,'(9999e23.15)') zone(m)%x(i,j,k), zone(m)%y(i,j,k), zone(m)%q(i,j,k,1)
        enddo
      enddo
    enddo


end subroutine
!####################################################################################
