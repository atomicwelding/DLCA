module params
  implicit none

  ! filepath to save the data
  character(len=*), parameter :: filepath = "data/L70.dat"

  ! linear dimension of the box
  ! volume => L^3
  integer, parameter :: L = 70


  ! parameters handling the mobility of the clusters
  ! alpha ≈ 1/d where d is the fractal dimension
  real, parameter :: alpha = 0.55


  ! number of phi to use
  integer, parameter :: Npts = 100
  
  ! phi in [min_phi, max_phi]
  real, parameter :: MIN_PHI = 0.1 , MAX_PHI = 0.

  ! number of runs per phi
  ! we are averaging over the runs
  integer, parameter :: runs = 1

  ! stopping criterion
  ! maximum of consecutive steps without anything happening
  ! in practice, L^3 is a good choice
  integer, parameter :: MAX_STEPS_WITHOUT_AGGREGATION = L*L*L

end module params
