!####################################################################################
!read grid data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_grid

  use mod_field
  use mod_CGNS
  use mod_pod_modes, only : working_zone
  implicit none
  integer(kind=4) :: m
  integer(kind=4) :: temp(3)

  integer(kind=4) :: nx1, nx2, nx3 

  write(*,'(A)') ' ' 
  write(*,'(A)') ' Reading grid data ...' 
  write(*,'(A)') ' ' 

  m = working_zone  

  write(*,'(A,i2,3(A5,i4))') ' Zone ', m, 'nx =', zone(m)%nx, 'ny =', zone(m)%ny, 'nz =', zone(m)%nz

  !~~~ Fazendo uma troca de valores aqui pq eu leio a malha toda SEMPRE...
  temp(1) = imax(m)
  temp(2) = jmax(m)
  temp(3) = kmax(m)

  imax(m) = zone(m)%nx
  jmax(m) = zone(m)%ny
  kmax(m) = zone(m)%nz 

  !~~~ A parte da malha eu aloco para td mundo mesmo.
  allocate(zone(m)%x(1:imax(m),1:jmax(m),1:kmax(m)))
  allocate(zone(m)%y(1:imax(m),1:jmax(m),1:kmax(m)))
  allocate(zone(m)%z(1:imax(m),1:jmax(m),1:kmax(m)))
  
  print*, shape( zone(m)%x )

  if (data_2D .eqv. .true.) then

    !~~ A leitura da malha assume imin = 1. Ou seja, lê a zona toda!!
    call read_grid_2D_CGNS( m,[zone(m)%nx,zone(m)%ny],trim(CGNS_gridname), &
          zone(m)%x(1:imax(m),1:jmax(m),1), zone(m)%y(1:imax(m),1:jmax(m),1) )
    
    !~~ Exportando a malha...   
    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    call create_file_CGNS(trim(output_path)//'/grid_mod_2D.cgns','2D')
    call write_partial_grid_2D_CGNS( output_zone,[nx1,nx2], &
          zone(m)%x(imin(m):imax(m),1:jmax(m),1), zone(m)%y(imin(m):imax(m),1:jmax(m),1), trim(output_path)//'/grid_mod_2D.cgns' )
  
  endif
  
  if (data_3D .eqv. .true.) then

    !~~ A leitura da malha é feita para a zona toda!!
    call read_grid_3D_CGNS( m,[zone(m)%nx,zone(m)%ny,zone(m)%nz],trim(CGNS_gridname), &
          zone(m)%x(1:imax(m),1:jmax(m),1:kmax(m)), zone(m)%y(1:imax(m),1:jmax(m),1:kmax(m)), zone(m)%z(1:imax(m),1:jmax(m),1:kmax(m)) )
 
    !~~ Mas eu só exporto um chunk da malha... 
    imax(m) = temp(1) 
    jmax(m) = temp(2)
    kmax(m) = temp(3)
     
    if (imax(m) .eq. 0) imax(m) = zone(m)%nx
    if (jmax(m) .eq. 0) jmax(m) = zone(m)%ny
    if (kmax(m) .eq. 0) kmax(m) = zone(m)%nz
    
    print*, ' MIN : ', m, imin(m), jmin(m), kmin(m)
    print*, ' MAX : ', m, imax(m), jmax(m), kmax(m)

    nx1 = imax(m)-imin(m)+1
    nx2 = jmax(m)
    nx3 = kmax(m)
    call create_file_CGNS(trim(output_path)//'/grid.cgns','3D')
    call write_partial_grid_3D_CGNS( output_zone,[nx1,nx2,nx3], &
          zone(m)%x(imin(m):imax(m),1:jmax(m),1:kmax(m)), zone(m)%y(imin(m):imax(m),1:jmax(m),1:kmax(m)), &
          zone(m)%z(imin(m):imax(m),1:jmax(m),1:kmax(m)), trim(output_path)//'/grid.cgns' )
      
  endif
  
  write(*,'(A)') ' Grid ... OK!'
  write(*,'(A)') ' '
 
  return
  
end subroutine
!####################################################################################

!####################################################################################
!read solution data for ALL zones
!in order to read specific zones, one should modify the routines in "cgns_routines.f90"
subroutine read_soln_ruvw()

  use mod_field
  use mod_CGNS
  use mod_pod_modes, only : working_zone
  use mod_signal_process, only : dt
  implicit none
  integer(kind=4) :: i, m, t
  integer(kind=4) :: idummy

  logical :: new_file_version

  character(len=16) :: cdummy
  
  nsnap = (idxf - idxi)/idxr + 1
  idxf = idxi + (nsnap-1)*idxr
  
  write(*,*) 'Solution files:'
  write(*,*) idxi, idxf, idxr
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  dt = dt*dble(fresult)*dble(idxr)
 
  allocate(time(1:nsnap))
  time(1:nsnap) = 0.0d0
  do i = 1,nsnap
    time(i) = dble(i-1)*dt
  enddo  

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  m = working_zone

    if (iblank .eqv. .true.) then
      allocate(zone(m)%iblank(imin(m):imax(m),1:jmax(m)))
               zone(m)%iblank(imin(m):imax(m),1:jmax(m)) = 0
    endif
  
    if (trim(working_var) .eq.  'Pressure' .or. trim(working_var) .eq.   'Density' .or. &
        trim(working_var) .eq. 'MomentumX' .or. trim(working_var) .eq. 'MomentumY' .or. trim(working_var) .eq. 'MomentumZ') then
      write(*,'(A,i0,A)') ' Reading ' // trim(working_var) // ' data for ', nsnap, ' snapshots ...'
      allocate(zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
    endif
    
    if (trim(working_var) .eq. 'VelocityX' .or. trim(working_var) .eq. 'VelocityY' .or. trim(working_var) .eq. 'VelocityZ') then
      write(*,'(A,i0,A)') ' Reading the "Momentum/Density" data for ' // trim(working_var) // ' for ', nsnap, ' snapshots ...' 
      allocate(zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
      allocate(zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
    endif    
    
    if (trim(working_var) .eq. "ReynoldsStress") then
      write(*,'(A,i0,A)') ' Reading the all the data for ', nsnap, ' snapshots ...' 
      allocate(zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
      allocate(zone(m)%u(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%u(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
      allocate(zone(m)%v(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%v(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
      allocate(zone(m)%w(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%w(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
      allocate(zone(m)%p(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
               zone(m)%p(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
!      allocate(zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
!               zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap) = 0.0d0
    endif
  
    !---- reading the solution files
    t = 0
    soln : do idummy = idxi,idxf,idxr

      t = t + 1
      ! !Note: .cgnc =  5 caracteres!, logo A5 em i6.6 --> numero de digitos do qout
      ! if (data_2D .eqv. .true.) write(CGNS_solnname,'(A,A,i6.6,A8)') trim(path_to_soln), 'qout', idummy, '_2D.cgns'
      ! !Dados CGNS do HD completo
      if (data_2D .eqv. .true.) write(CGNS_solnname,'(A,A,i6.6,A5)') trim(path_to_soln), 'qout', idummy, '.cgns'
      if (data_3D .eqv. .true.) write(CGNS_solnname,'(A,A,i6.6,A5)') trim(path_to_soln), 'qout', idummy, '.cgns'
      !Reconstrução
      !if (data_3D .eqv. .true.) write(CGNS_solnname,'(A,A,i4.4,A5)') trim(path_to_soln), 'rcns', idummy, '.cgns'
      !Dados CGNS do HD reduzido
      ! if (data_2D .eqv. .true.) write(CGNS_solnname,'(A,A,i6.6,A8)') trim(path_to_soln), 'qout', idummy, '_mod.cgns'
      ! if (data_3D .eqv. .true.) write(CGNS_solnname,'(A,A,i6.6,A9)') trim(path_to_soln), 'qout', idummy, '_mod.cgns'
      
      ! ! Checa se o arquivo esta no padrão novo (Tulio, Brener e Miotto) ou no padrão antigo (William e Jean)
      ! inquire(file=trim(CGNS_solnname),exist=new_file_version)
      ! if (new_file_version .eqv. .true.) then
        
      !   if (t .eq. 1) print*,  '(1) Arquivo do William !!'

      !   if (data_2D .eqv. .true.) write(CGNS_solnname,'(A,A,i6.6,A8)') trim(path_to_soln), 'qout', idummy, '_2D.cgns'
      !   if (data_3D .eqv. .true.) write(CGNS_solnname,'(A,A,i4.4,A5)') trim(path_to_soln), 'qout', idummy, '.cgns'
        
      ! endif

      !~~~~~~
      ! Lê o iblank se for necessário
      if (t .eq. 1 .and. iblank .eqv. .true.) then
        write(*,'(A)') ' '
        write(*,'(A)') ' Reading iblank ...'
        call read_partial_INT_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),1], &
              'Iblank',trim(CGNS_solnname),zone(m)%iblank(imin(m):imax(m),1:jmax(m)))
        write(*,'(A)') ' Iblank OK!'

        if (trim(working_var) .eq. "Grid") then
          print*, ' Reading grid and Iblank only...'      
          exit soln
        endif

      endif               
      
      !~~~~~~
      if (trim(working_var) .ne. "ReynoldsStress" .and. working_var(1:8) .ne. "Velocity" .and. trim(working_var) .ne. "Grid" ) then
      
        ! Lê uma variável existente no arquivo CGNS versão nova (6 digitos, e momento)
        if (new_file_version .eqv. .true.) then
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)], &
                        trim(working_var),trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
        endif

        ! ! Lê uma variável existente no arquivo CGNS na versão antiga (4 digitos e velocidade)
        ! if (new_file_version .eqv. .true.) then
          
        !   if (t .eq. 1) then
        !     if (trim(working_var) .eq. "MomentumX") cdummy = 'VelocityX'
        !     if (trim(working_var) .eq. "MomentumY") cdummy = 'VelocityY'
        !     if (trim(working_var) .eq. "MomentumZ") cdummy = 'VelocityZ'
        !     if (trim(working_var) .eq.  "Pressure") cdummy = 'Pressure'
        !     print*,  '(2) Arquivo do William !!', trim(cdummy)
        !   endif
          
        !   call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)], &
        !                 trim(cdummy),trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

        ! endif

      endif

      !Bloco que geralmente eh utilizado, comentar caso for rodar o do william
       
      if (new_file_version .eqv. .true.) then

        if ( working_var(1:8) .eq. "Velocity"  ) &
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],  'Density',trim(CGNS_solnname),zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
          
        if (trim(working_var) .eq. "VelocityX" ) &
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'MomentumX',trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

        if (trim(working_var) .eq. "VelocityY" ) &
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'MomentumY',trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

        if (trim(working_var) .eq. "VelocityZ" ) &
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'MomentumZ',trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))  
        
      endif
      
      
    !  if (new_file_version .eqv. .true.) then
       
    !    if (t .eq. 1) print*,  '(3) Arquivo do William !!'

    !    if (trim(working_var) .eq. "Density"  ) &
    !      call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],  'Density',trim(CGNS_solnname),zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

    !    if (trim(working_var) .eq. "VelocityX" ) &
    !      call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'VelocityX',trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

    !    if (trim(working_var) .eq. "VelocityY" ) &
    !      call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'VelocityY',trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

    !    if (trim(working_var) .eq. "VelocityZ" ) &
    !      call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'VelocityZ',trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t)) 

    !      print*, m, imin(m), imax(m)

    !      call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)], &
    !                    trim(working_var),trim(CGNS_solnname),zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
     
    !  endif
        
      !~~~~~~
      if (trim(working_var) .eq. "ReynoldsStress") then

        if (new_file_version .eqv. .true.) then  !Bloco utilizado para as arquivos que não são do William

          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],  'Density',trim(CGNS_solnname),zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t))      
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'MomentumX',trim(CGNS_solnname),zone(m)%u(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'MomentumY',trim(CGNS_solnname),zone(m)%v(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'MomentumZ',trim(CGNS_solnname),zone(m)%w(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
          call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)], 'Pressure',trim(CGNS_solnname),zone(m)%p(imin(m):imax(m),1:jmax(m),1:kmax(m),t))    

        endif
        
        ! if (new_file_version .eqv. .true.) then
        
        !   if (t .eq. 1) print*,  '(4) Arquivo do William !!'

        !  zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = 1.0d0

        !   call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],  'Density',trim(CGNS_solnname),zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
        !   call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'VelocityX',trim(CGNS_solnname),zone(m)%u(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
        !   call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'VelocityY',trim(CGNS_solnname),zone(m)%v(imin(m):imax(m),1:jmax(m),1:kmax(m),t))
        !   call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)],'VelocityZ',trim(CGNS_solnname),zone(m)%w(imin(m):imax(m),1:jmax(m),1:kmax(m),t))  
        !   call read_partial_soln_CGNS(m,[imin(m),1,1],[imax(m),jmax(m),kmax(m)], 'Pressure',trim(CGNS_solnname),zone(m)%p(imin(m):imax(m),1:jmax(m),1:kmax(m),t))

        ! endif
      
        endif   
    
    enddo soln

  ! !XXX Migue
  !     call system('printf "\033[0;32m \n"')
  !     print*, ' --> Migue :: Estou alocando "q" pq eu uso no calculo das flutuacoes!' 
  !     print*, ' --> Migue :: Tentar pensar numa forma mais organizada ...' 
  !     call system('printf "\033[0m \n"') 
    print*, ' '
          
  if (allocated(zone(m)%q) .eqv. .false.) allocate(zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),1:nsnap))
       
!  print*, ' '
!  print*, ' cdummy = ', trim(cdummy)
!  print*, ' cdummy = ', trim(cdummy)
!  print*, ' cdummy = ', trim(cdummy)
!  print*, ' cdummy = ', trim(cdummy)
!  print*, ' cdummy = ', trim(cdummy)
!  print*, ' '

  return

end subroutine read_soln_ruvw
!####################################################################################

!####################################################################################
subroutine link_iblank

  use mod_field
  use mod_CGNS
  use mod_pod_modes, only : working_zone
  implicit none
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2  
  
  m = working_zone
  
  nx1 = imax(m)-imin(m)+1
  nx2 = jmax(m)  

  call create_file_CGNS('./iblank.cgns','2D')
  call write_iblank_2D_CGNS(m,[nx1,nx2],zone(m)%iblank(imin(m):imax(m),1:jmax(m)),'./iblank.cgns')
  call write_link_CGNS(m,'./iblank.cgns','./grid_mod_2D.cgns')
 
  return

end subroutine link_iblank
!####################################################################################

!####################################################################################
subroutine convert_momentum_to_velocity

  use mod_field
  use mod_CGNS
  use mod_pod_modes, only : working_zone
  implicit none
  integer(kind=4) :: m, t

  m = working_zone

    if (trim(working_var) .eq. "VelocityX" .or. trim(working_var) .eq. "VelocityY" .or. trim(working_var) .eq. "VelocityZ" .or. &
        trim(working_var) .eq. "ReynoldsStress") then

    do t = 1,nsnap
    
      write(*,'(1a1,A,i0,$)') char(13), ' Converting from momentum to velocity the snapshot number ', t

      if (trim(working_var) .eq. "VelocityX" ) &
        zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t) / zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t)
      if (trim(working_var) .eq. "VelocityY" ) &
        zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t) / zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t)
      if (trim(working_var) .eq. "VelocityZ" ) &
        zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = zone(m)%q(imin(m):imax(m),1:jmax(m),1:kmax(m),t) / zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t)
        
      if (trim(working_var) .eq. "ReynoldsStress" ) then
        zone(m)%u(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = &
                                 zone(m)%u(imin(m):imax(m),1:jmax(m),1:kmax(m),t) / zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t)
        zone(m)%v(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = &
                                 zone(m)%v(imin(m):imax(m),1:jmax(m),1:kmax(m),t) / zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t)
        zone(m)%w(imin(m):imax(m),1:jmax(m),1:kmax(m),t) = &
                                 zone(m)%w(imin(m):imax(m),1:jmax(m),1:kmax(m),t) / zone(m)%r(imin(m):imax(m),1:jmax(m),1:kmax(m),t)
      endif

    enddo
  
  endif

  write(*,*) ' '

  return

end subroutine convert_momentum_to_velocity
!####################################################################################
