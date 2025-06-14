! main program
! provided by weld
program DLCA
  use types
  use params
  use subroutines
  implicit none

  type(SimulationState) :: S
  real :: phi = min_phi
  integer :: idx_phi, idx_run
  integer :: actives
  integer :: prev_max_cluster
  real    :: tgel_this, tgel_total
  integer :: number_percolate
  real :: Pgel
  integer :: i

  real :: t = 0.0

  integer :: Nc = 0 

  integer :: seed_size
  integer, allocatable :: seed(:)
  integer :: base_seed

  integer :: fd, io
  character(len=100) :: errmsg
  
  open(fd, file=filepath, status="replace", iostat=io, iomsg=errmsg)
  if (io /= 0) then
     print *, "Error opening file:", trim(errmsg)
     stop 1
  end if

  ! file header
  write(fd, *) "@ L = ", L
  write(fd, *) "@ MAX_STEPS_WITHOUT_AGGREGATION = ", MAX_STEPS_WITHOUT_AGGREGATION
  write(fd, *) "@ t, cluster size" 
  call flush(fd)


  ! infer the number of monomers based on density phi
  phi = min_phi
  S%N = int(phi * L**3)

  allocate(S%particles(0:S%N - 1))
  allocate(S%cluster_size(0:S%N - 1))

  print *, "Running phi =", phi, "=> N =", S%N

  ! important to setup the seed
  ! if not, results are always the same
  call random_seed(size = seed_size)
  allocate(seed(seed_size))

  S%grid = -1
  S%cluster_size = 0
  do i = 0, S%N - 1
     S%particles(i)%cluster_id = -1
  end do
  call init_particles(S)

  S%no_growth_steps = 0
  
  do while ( t < MAX_TIME )
     call trial(S)
     Nc = count_active_clusters(S)
     t = t + 1/dble(Nc)
     if( mod(t, 1.0) == 0 ) then
        print *, "t: ", t
        write(fd, '(F7.2)', advance='no') t
        do i = 0, S%N - 1
           write(fd, '(1X,I6)', advance='no') S%cluster_size(i)
        end do
        write(fd, *)
     end if
  end do

  deallocate(seed)
  deallocate(S%particles)
  deallocate(S%cluster_size)

  close(fd)

end program DLCA
