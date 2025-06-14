module params
  implicit none

  ! filepath to save the data
  character(len=*), parameter :: filepath = "L20LLL.dat"

  ! linear dimension of the box
  ! volume => L^3
  integer, parameter :: L = 20


  ! parameters handling the mobility of the clusters
  ! alpha ≈ 1/d where d is the fractal dimension
  real, parameter :: alpha = 0.55


  ! number of phi to use
  integer, parameter :: Npts = 40
  
  ! phi in [min_phi, max_phi]
  real, parameter :: MIN_PHI = 0.0 , MAX_PHI = 0.3

  ! number of runs per phi
  ! we are averaging over the runs
  integer, parameter :: runs = 48

  ! stopping criterion
  ! maximum of consecutive steps without anything happening
  ! in practice, L^3 is a good choice
  integer, parameter :: MAX_STEPS_WITHOUT_AGGREGATION = 8000 ! L*L*L

end module params
