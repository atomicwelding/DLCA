! main program
! provided by weld
program DLCA
  use types
  use params
  use subroutines
  implicit none

  type(SimulationState) :: S
  real, dimension(0:Npts-1) :: phis
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

  ! build the list of phis
  do idx_phi = 0, Npts - 1
     phis(idx_phi) = MIN_PHI +  real(idx_phi) * (MAX_PHI - MIN_PHI)/ real(Npts - 1)
  end do

  open(fd, file=filepath, status="replace", iostat=io, iomsg=errmsg)
  if (io /= 0) then
     print *, "Error opening file:", trim(errmsg)
     stop 1
  end if

  ! file header
  write(fd, *) "@ L = ", L
  write(fd, *) "@ MAX_STEPS_WITHOUT_AGGREGATION = ", MAX_STEPS_WITHOUT_AGGREGATION
  write(fd, *) "@ runs = ", runs
  write(fd, *) "@ phi0, P, tgel_avg, Pgel_avg"
  call flush(fd)


  ! loop over each phi value, sequentially
  do idx_phi = 0, Npts - 1

     ! infer the number of monomers based on density phi
     S%N = int(phis(idx_phi) * L**3)

     ! If N < L, percolation is impossible, write trivial statement 
     if (S%N < L) then
        write(fd, "(F8.6, 4(F12.6))")     &
             phis(idx_phi),               &
             0.0,                         & ! P = 0
             -1.0,                         & ! tgel_avg = -1
             0.0                           ! Pgel_avg = 0
        call flush(fd)
        cycle
     end if

     allocate(S%particles(0:S%N - 1))
     allocate(S%cluster_size(0:S%N - 1))

     print *, "Running phi =", phis(idx_phi), "=> N =", S%N

     tgel_total = 0.0
     number_percolate = 0
     Pgel = 0.0


     ! important to setup the seed
     ! if not, results are always the same
     call random_seed(size = seed_size)
     allocate(seed(seed_size))


     ! parallelize this runs loop with OpenMP
     ! it reates an overhead when only few runs on low densities
     ! please prefer to use a single thread in this case
     !$OMP PARALLEL DO default(shared)                                      &
     !$OMP   private(idx_run, base_seed, tgel_this, prev_max_cluster, actives, i) &
     !$OMP   firstprivate(S, seed)                                            &
     !$OMP   reduction(+: number_percolate, tgel_total, Pgel)
     do idx_run = 1, runs


        ! chatgpt gave me a method to put more entropy into
        ! my random number generator ;D
        call system_clock(count=base_seed)
        base_seed = base_seed + idx_run + 1000 * idx_phi
        seed = base_seed

        call random_seed(put = seed)

        ! reset state for this run
        S%hasEnded      = .false.
        S%hasPercolated = .false.
        S%grid = -1
        S%cluster_size = 0
        do i = 0, S%N - 1
           S%particles(i)%cluster_id = -1
        end do

        call init_particles(S)

        tgel_this = 0.0
        S%no_growth_steps = 0
        prev_max_cluster = 1

        do while (.not. S%hasEnded)
           actives = count_active_clusters(S)
           tgel_this = tgel_this + 1.0 / real(actives)
           call trial(S)

           if (maxval(S%cluster_size) <= prev_max_cluster) then
              S%no_growth_steps = S%no_growth_steps + 1
           else
              S%no_growth_steps = 0
              prev_max_cluster = maxval(S%cluster_size)
           end if
        end do

        if (S%hasPercolated) then
           number_percolate = number_percolate + 1
           tgel_total = tgel_total + tgel_this
           Pgel = Pgel + real(maxval(S%cluster_size)) / real(S%N)
        end if

     end do
     !$OMP END PARALLEL DO

     deallocate(seed)

     ! write out
     if (number_percolate == 0) then
        write(fd,*) phis(idx_phi), 0.0, 0.0, -1.0
     else
        write(fd,*) phis(idx_phi),                       & 
             real(number_percolate)/real(runs),           &
             tgel_total/real(number_percolate),           &
             Pgel/real(number_percolate)
     end if
     call flush(fd)

     deallocate(S%particles)
     deallocate(S%cluster_size)

  end do  ! end loop over phis

  close(fd)

end program DLCA
