!#################################################################################
subroutine compute_cell_volume(zone_number)
 
  use mod_field
  use mod_pod_modes
  implicit none
  integer(kind=4), intent(in) :: zone_number  
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: m

  m = zone_number

  call compute_cell_area(m)

  if (iblank .eqv. .false.) return

  if (iblank .eqv.  .true.) then
  
    write(*,'(A)') ' Working with Iblank ...' 
  
    !2D case
    if (data_2D .eqv. .true.) then

      do j = 1,jmax(m)
        do i = imin(m),imax(m)
          if (zone(m)%iblank(i,j) .ne. 10) zone(m)%area(i,j) = 0.0d0
        enddo
      enddo
      
    endif
    
    !3D case
    if (data_3D .eqv. .true.) then
    
      write(*,'(A)') ' Computing cell volume ...' 

      allocate(zone(m)%volume(imin(m):imax(m),1:jmax(m)))
               zone(m)%volume(imin(m):imax(m),1:jmax(m)) = 0.0d0

      do j = 1,jmax(m)
        do i = imin(m),imax(m)
          if (zone(m)%iblank(i,j) .eq. 10) zone(m)%volume(i,j) = zone(m)%area(i,j)*dz(m)
          if (zone(m)%iblank(i,j) .ne. 10) zone(m)%volume(i,j) = 0.0d0
        enddo
      enddo
      
    endif
    
  endif

  return

end subroutine compute_cell_volume
!=================================================================================

!=================================================================================
subroutine compute_cell_area(zone_number)
 
  ! the area is for a ficticious cell-vertex element
  ! the "center" is the node from the finite difference scheme
  ! the total area is the sum of portions of the "cell center" surrounding elements

  use mod_field
  implicit none
  integer(kind=4), intent(in) :: zone_number
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: m
  integer(kind=4) :: nx1, nx2
  
  real(kind=8) xaux(4), yaux(4), area

  write(*,'(A)') ' Computing cell area ...'   

  m = zone_number
  
    nx1 = imax(m)
    nx2 = jmax(m)

    allocate(zone(m)%area(1:nx1,1:nx2))
    zone(m)%area(:,:) = 0.0d0

    do j = 1,nx2-1
  
      do i = 1,nx1-1
      
        xaux(1) = zone(m)%x( i , j ,1)
        xaux(2) = zone(m)%x(i+1, j ,1)
        xaux(3) = zone(m)%x(i+1,j+1,1)
        xaux(4) = zone(m)%x( i ,j+1,1)
        
        yaux(1) = zone(m)%y( i , j ,1)
        yaux(2) = zone(m)%y(i+1, j ,1)
        yaux(3) = zone(m)%y(i+1,j+1,1)
        yaux(4) = zone(m)%y( i ,j+1,1)
        
        call compute_quad_area(xaux(1:4),yaux(1:4),area)       
        
        zone(m)%area( i , j ) = zone(m)%area( i , j ) + area*0.25d0
        zone(m)%area(i+1, j ) = zone(m)%area(i+1, j ) + area*0.25d0
        zone(m)%area(i+1,j+1) = zone(m)%area(i+1,j+1) + area*0.25d0
        zone(m)%area( i ,j+1) = zone(m)%area( i ,j+1) + area*0.25d0

      enddo
      
    enddo
 
    if (nzones .eq. 2 .and. m .eq. 1) then
    
      !~~~ Periodic boundary conditions
      !
      !  This is valid when the implicit region is a single zone
   
        call system('printf "\033[0;32m \n"')   
      write(*,'(A,i0)') ' --> Using periodic boundary conditions for zone ', m
      write(*,'(/,A )') ' --> This is usually valid when the implicit region is a single zone !'
      write(*,'(  A )') ' --> Also, make sure that the gap is closed !!! '
        call system('printf "\033[0m \n"')

      ! Corrigindo o gap
      zone(m)%area(nx1,:) = zone(m)%area(nx1,:) + zone(m)%area(1,:)
      zone(m)%area( 1 ,:) = zone(m)%area(nx1,:)

    endif

    ! Primeiro e ultimo elemento na direção normal à parede 
    zone(m)%area(:, 1 ) = 2.0d0*zone(m)%area(:, 1 )
    zone(m)%area(:,nx2) = 2.0d0*zone(m)%area(:,nx2)    
    
  !enddo
  
!    do j = 1,jmax(m)
!      do i = imin(m),imax(m)
!        !TODO mudar isso aqui pra quando não existir Iblank (Nagarajan)
!        if (zone(m)%iblank(i,j) .ne. 10) zone(m)%area(i,j) = 0.0d0
!      enddo
!    enddo  

!  open(99,file='/home/tulio/Desktop/fort.dat',position='append')
!  !open(99,file='/home/tulio/Desktop/fort.dat')
!  !write(99,'(A)') 'VARIABLES="X","Y","S","Iblank"'
!  write(99,'(A)') 'VARIABLES="X","Y","S"'
!  !do m = 1,nzones

!    nx1 = zone(m)%nx
!    nx2 = zone(m)%ny

!    write(99,'(A,i0,A,i0,A)') 'ZONE T="",I=',nx1,',J=',nx2,',F=POINT'

!    do j = 1,nx2
!     do i = 1,nx1
!       !write(99,*) zone(m)%x(i,j,1),zone(m)%y(i,j,1),zone(m)%area(i,j),zone(m)%iblank(i,j)
!       write(99,*) zone(m)%x(i,j,1),zone(m)%y(i,j,1),zone(m)%area(i,j)
!     enddo
!    enddo

!  !enddo

  return

end subroutine compute_cell_area
!=================================================================================

!=================================================================================
!
!   Area of a quadrilateral element according to the right-hand rule:
!
!                    A = index 1
!                    B = index 2
!                    C = index 3
!                    D = index 4
!
!              D4____________________ C3 
!               /                    |
!              /                     |
!             /                      |
!            /                       |
!           /________________________|
!         A1                         B2
!       
!
!   https://en.wikipedia.org/wiki/Quadrilateral#Vector_formulas
!
subroutine compute_quad_area(x,y,area)

  implicit none
  real(kind=8), intent(in)  :: x(4), y(4)
  real(kind=8), intent(out) :: area
  real(kind=8) x1,x2,y1,y2

  x1 = x(3) - x(1)
  y1 = y(3) - y(1)
  
  x2 = x(4) - x(2)
  y2 = y(4) - y(2)

  area = 0.5*dabs(x1*y2 - x2*y1)
  
  return

end subroutine
!#################################################################################
