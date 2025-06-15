! main program
! provided by weld
program DLCA
  use types
  use params
  use subroutines
  implicit none

  type(SimulationState) :: S
  integer :: t = 0, MAX_TIME = 20000
  real :: phi = min_phi
  integer :: idx_phi, idx_run
  integer :: actives
  integer :: prev_max_cluster
  real    :: tgel_this, tgel_total
  integer :: number_percolate
  real :: Pgel
 integer :: i

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

  call write_frame(S, 0)
  do while (S%hasEnded .NEQV. .true.)
     t = t + 1 
     call trial(S)
     if(mod(t, 100) == 0) then
       call write_frame(S, t)
     end if
  end do

  deallocate(seed)
  deallocate(S%particles)
  deallocate(S%cluster_size)

  close(fd)



  
contains
  subroutine write_frame(S, current_frame)
    type(SimulationState), intent(inout) :: S
    integer, intent(in) :: current_frame
    integer :: i
    character, parameter :: symbol = 'N'

    write(fd, *) S%N
    write(fd, *) "Frame ", current_frame
    do i = 0, S%N-1
       write(fd, *) "N", real(S%particles(i)%x), real(S%particles(i)%y), real(S%particles(i)%z)
    end do
  end subroutine write_frame
end program DLCA
