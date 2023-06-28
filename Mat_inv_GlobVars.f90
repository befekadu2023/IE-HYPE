MODULE Mat_inv_GlobVars

double precision, dimension(:), allocatable ::  uxx   !real matrix (dim_jacob)
double precision, dimension(:), allocatable ::  uxx_SoilP   !real matrix (dim_jacob_SoilP)
double precision, dimension(:,:), allocatable ::  JacobMatrx  !real matrix (dim_jacob x dim_jacob)
double precision, dimension(:,:), allocatable ::  JacobMatrx_SoilP  !real matrix (dim_jacob_SoilP x dim_jacob_SoilP)
  double precision, dimension(:,:), allocatable ::  JacobMatrxCopy
  double precision, dimension(:,:), allocatable ::  JacobMatrxCopy_SoilP  !copy of matrix JacobMatrx
  double precision, dimension(:,:), allocatable ::  inv_Jacob
  double precision, dimension(:,:), allocatable ::  inv_Jacob_SoilP   !real matrix (dim_jacob_SoilP x dim_jacob_SoilP)
  integer, allocatable ::  row_permutat(:)  !integer vector (dim_jacob)
  integer, allocatable ::  row_permutat_SoilP(:)  !integer vector (dim_jacob_SoilP)
  double precision, dimension(:), allocatable ::  matrx_multip
  double precision, dimension(:), allocatable ::  matrx_multip_SoilP   !real matrix (dim_jacob)
  double precision, dimension(:), allocatable ::  S_new, S_old,S_t0	!!real matrix (dim_jacob)
  double precision, dimension(:), allocatable ::  S_new_SoilP, S_old_SoilP,S_t0_SoilP	!!real matrix (dim_jacob_SoilP)
  double precision, dimension(:), allocatable :: F_numerator, F_numerator_SoilP
  double precision :: dt
  integer :: dim_jacob,CODE_mtrx
  integer :: dim2_ode
  


  
  END MODULE Mat_inv_GlobVars