module types
  use params
  implicit none
  type :: Particle ! represents a monomer
     integer :: id ! unique id
     integer :: x, y, z ! coordinates in the box
     integer :: m = 1   ! mass of the particle 
     integer :: cluster_id ! id of the cluster it's in
  end type Particle

  type FaceTouch ! utility type, to check for percolation
     logical :: min_x = .false.
     logical :: max_x = .false.
     logical :: min_y = .false.
     logical :: max_y = .false.
     logical :: min_z = .false.
     logical :: max_z = .false.
  end type FaceTouch


  type :: SimulationState
     integer :: N = 0
     integer, dimension(0:L-1, 0:L-1, 0:L-1) :: grid ! each site holds either -1 or id

     type(Particle), allocatable :: particles(:)
     integer,        allocatable :: cluster_size(:)

     integer :: no_growth_steps = 0 ! used to check if nothing happened

     ! flags
     logical :: hasPercolated  = .false.
     logical :: hasEnded       = .false.

     logical :: hasAdvanced    = .false.

  end type SimulationState
end module types
