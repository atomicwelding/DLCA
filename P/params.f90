module params
  implicit none

  ! filepath to save the data
  character(len=*), parameter :: filepath = "data/L05.dat"

  ! linear dimension of the box
  ! volume => L^3
  integer, parameter :: L = 05


  ! parameters handling the mobility of the clusters
  ! alpha ≈ 1/d where d is the fractal dimension
  real, parameter :: alpha = 0.55


  ! number of phi to use
  integer, parameter :: Npts = 40
  
  ! phi in [min_phi, max_phi]
  real, parameter :: MIN_PHI = 0.0 , MAX_PHI = 0.3

  ! number of runs per phi
  ! we are averaging over the runs
  integer, parameter :: runs = 9600

  ! stopping criterion
  ! maximum of consecutive steps without anything happening
  ! in practice, L^3 is a good choice
  integer, parameter :: MAX_STEPS_WITHOUT_AGGREGATION = 10000

end module params
