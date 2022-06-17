!#####################################################################################
module mod_field

  real(kind=8) :: Re, Ma, Reinv, gas_cte, gamma, pinf, Twall, mu

  logical :: data_3D
  logical :: data_2D
  
  logical :: iblank
  !logical :: use_r, use_u, use_v, use_w, use_q  

  character(len=132) :: path_to_grid
  character(len=132) :: path_to_soln
  character(len=132) :: output_path
  character(len=132) :: grid_name

  integer(kind=4) idxi     ! index initial
  integer(kind=4) idxf     ! index final
  integer(kind=4) idxr     ! index rate
  integer(kind=4) nsnap
  integer(kind=4) :: fresult
  
  integer(kind=4) :: nzones
  integer(kind=4) :: output_zone
  
!  logical         :: logical_primitive

  integer(kind=4) :: start_index_i, final_index_i
  integer(kind=4) :: start_index_j, final_index_j  
  integer(kind=4), allocatable, dimension(:) :: imin, imax
  integer(kind=4), allocatable, dimension(:) :: jmin, jmax
  integer(kind=4), allocatable, dimension(:) :: kmin, kmax
  
  real(kind=8), allocatable :: dz(:)

  type zone_vars

    real(kind=8), allocatable, dimension(:,:,:,:) :: q

    real(kind=8), allocatable, dimension(:,:)     :: qmean2D
!    real(kind=8), allocatable, dimension(:,:)     :: qstdd2D
        
    real(kind=8), allocatable, dimension(:,:,:)   :: qmean3D

    real(kind=8), allocatable, dimension(:,:,:,:) :: r
    real(kind=8), allocatable, dimension(:,:,:,:) :: u
    real(kind=8), allocatable, dimension(:,:,:,:) :: v
    real(kind=8), allocatable, dimension(:,:,:,:) :: w
    real(kind=8), allocatable, dimension(:,:,:,:) :: p

!    real(kind=8), allocatable, dimension(:,:,:,:) :: kinetic
!    real(kind=8), allocatable, dimension(:,:,:,:,:) :: quadrupole

    complex(kind=8), allocatable, dimension(:,:,:,:) :: q_hat
    
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: u_hat
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: v_hat
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: w_hat
!    complex(kind=8), allocatable, dimension(:,:,:,:) :: p_hat    
    
    real(kind=8), allocatable, dimension(:,:,:) :: x, y, z
    
    integer(kind=4), allocatable :: iblank(:,:)    
    real(kind=8), allocatable :: area(:,:)
    real(kind=8), allocatable :: volume(:,:)
    
    integer(kind=4) :: nx, ny, nz
!    integer(kind=4) :: nx1, nx2, nx3
    
  end type zone_vars
  type(zone_vars), allocatable :: zone(:)

end module
!#####################################################################################
module mod_pod_modes

  integer(kind=4) :: working_zone
  integer(kind=4) :: zone_min, zone_max

  !~~~~~ sparsity!
  integer(kind=4), parameter :: istride = 2
  integer(kind=4), parameter :: jstride = 2
  integer(kind=4), parameter :: kstride = 2

  !~~~~~ 
  integer(kind=4) nPODmodes
  integer(kind=4) nFourierModes

  real(kind=8), allocatable :: lambda(:)
  real(kind=8), allocatable :: temporal_modes(:,:)

  type pod_zones
    real(kind=8), allocatable :: spatial_modes(:,:,:,:)
    real(kind=8), allocatable :: reconstructed(:,:,:,:)
  end type pod_zones
  type(pod_zones), allocatable :: pod_zone(:)
  
  character(len=250) :: path_to_output_matrix

  character(len=16) :: POD_operation

end module
!#####################################################################################
