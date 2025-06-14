module subroutines
  use types
  use params
  implicit none
contains
  subroutine init_particles(S)
    type(SimulationState), intent(inout) :: S
    integer :: i, x, y, z
    real    :: r

    do i = 0, S%N - 1
10     call random_number(r)
       x = int(r * L)
       call random_number(r)
       y = int(r * L)
       call random_number(r)
       z = int(r * L)
       if (S%grid(x, y, z) /= -1) goto 10  ! keep trying until we find an empty cell

       S%particles(i)%id         = i
       S%particles(i)%x          = x
       S%particles(i)%y          = y
       S%particles(i)%z          = z
       S%particles(i)%cluster_id = i
       S%cluster_size(i)         = 1
       S%grid(x, y, z)           = i
    end do
  end subroutine init_particles


  ! count_active_clusters: return how many cluster IDs have size > 0
  integer function count_active_clusters(S) result(res)
    type(SimulationState), intent(in) :: S
    integer :: i

    res = 0
    do i = 0, S%N - 1
       if (S%cluster_size(i) > 0) res = res + 1
    end do
  end function count_active_clusters


  ! flood_fill_all:
  ! reassign each particle’s cluster_id by breadth‐first search
  ! on S%grid; also marks visited so each cluster gets a new ID.
  subroutine flood_fill_all(S)
    type(SimulationState), intent(inout) :: S

    integer :: i, cid
    logical, dimension(0:S%N-1) :: visited
    integer, dimension(:), allocatable :: queue
    integer :: qstart, qend
    integer :: current_id, neighbor_id
    integer :: x, y, z, nx, ny, nz, j
    integer, dimension(6,3) :: dirs


   
    ! the six neighbor directions: ±x, ±y, ±z
    dirs = reshape([ &
         1,  0,  0,  &
         -1,  0,  0,  &
         0,  1,  0,  &
         0, -1,  0,  &
         0,  0,  1,  &
         0,  0, -1   &
         ], [6,3])

    ! if no particles, there's nothing to do
    if (S%N <= 0) return

    visited = .false.
    cid = 0
    allocate(queue(S%N))

    do i = 0, S%N - 1
       if (.not. visited(i)) then
          qstart = 1
          qend   = 1
          queue(qstart) = i
          visited(i) = .true.
          S%particles(i)%cluster_id = cid

          do while (qstart <= qend)
             current_id = queue(qstart)
             qstart = qstart + 1

             x = S%particles(current_id)%x
             y = S%particles(current_id)%y
             z = S%particles(current_id)%z

             do j = 1, 6
                nx = x + dirs(j,1)
                ny = y + dirs(j,2)
                nz = z + dirs(j,3)

                if (nx <  0     .or. nx >= L) cycle
                if (ny <  0     .or. ny >= L) cycle
                if (nz <  0     .or. nz >= L) cycle

                neighbor_id = S%grid(nx, ny, nz)
                if (neighbor_id /= -1) then
                   if (.not. visited(neighbor_id)) then
                      visited(neighbor_id) = .true.
                      S%particles(neighbor_id)%cluster_id = cid
                      qend = qend + 1
                      queue(qend) = neighbor_id
                   end if
                end if

             end do
          end do

          cid = cid + 1
       end if
    end do

    deallocate(queue)
  end subroutine flood_fill_all



  ! update_cluster_sizes:
  ! given updated particle‐>cluster_id, recompute sizes of every clusters
  subroutine update_cluster_sizes(S)
    type(SimulationState), intent(inout) :: S
    integer :: i

    

    S%cluster_size = 0
    do i = 0, S%N - 1
       S%cluster_size(S%particles(i)%cluster_id) = &
            S%cluster_size(S%particles(i)%cluster_id) + 1
    end do
  end subroutine update_cluster_sizes


  ! checkEndSimulation:
  ! check if simulation has ended, with 3 criteria,
  ! 1 - if all particles belong to a unique cluster
  ! 2 - if the system has percolated, i.e. a cluster spans over two
  ! opposite faces
  ! 3 - if it exceeds a maximum nb of steps without doing anything
  subroutine checkEndSimulation(S)
    type(SimulationState), intent(inout) :: S

    integer :: i, cid
    type(FaceTouch), dimension(0:S%N-1) :: touches

    ! clear touch flags
    touches(:)%min_x = .false.
    touches(:)%max_x = .false.
    touches(:)%min_y = .false.
    touches(:)%max_y = .false.
    touches(:)%min_z = .false.
    touches(:)%max_z = .false.

    ! 1st criterion
    if (count_active_clusters(S) == 1) then
       S%hasEnded = .true.
       return
    end if

    ! mark which faces are touched by each particle cluster_id:
    do i = 0, S%N - 1
       cid = S%particles(i)%cluster_id
       if (S%particles(i)%x == 0 )    touches(cid)%min_x = .true.
       if (S%particles(i)%x == L-1)  touches(cid)%max_x = .true.
       if (S%particles(i)%y == 0 )    touches(cid)%min_y = .true.
       if (S%particles(i)%y == L-1)  touches(cid)%max_y = .true.
       if (S%particles(i)%z == 0 )    touches(cid)%min_z = .true.
       if (S%particles(i)%z == L-1)  touches(cid)%max_z = .true.
    end do

    ! 2nd criterion
    do cid = 0, S%N - 1
       if ((touches(cid)%min_x .and. touches(cid)%max_x) .or. &
            (touches(cid)%min_y .and. touches(cid)%max_y) .or. &
            (touches(cid)%min_z .and. touches(cid)%max_z)) then
          S%hasPercolated = .true.
          S%hasEnded       = .true.
          return
       end if
    end do

    ! 3rd criterion
    if (S%no_growth_steps >= MAX_STEPS_WITHOUT_AGGREGATION) then
       S%hasEnded = .true.
       return 
    end if
  end subroutine checkEndSimulation


  ! trial:
  ! pick a random cluster (with probability ∝ 1/m^-α)
  ! attempt move,
  ! then flood_fill_all, update_cluster_sizes, checkEndSimulation
  subroutine trial(S)
    implicit none

    type(SimulationState), intent(inout) :: S

    integer :: num_clusters, i, cid, d, cluster_idx
    integer :: dx, dy, dz, j, nx, ny, nz, neighbor_id
    integer :: x_old, y_old, z_old
    integer :: members_count

    integer, allocatable :: active_clusters(:)
    integer, allocatable :: members(:)

    real    :: r, prob
    logical :: can_move, hasAggregated
    integer, dimension(6,3) :: dirs

    S%hasAdvanced = .false.

    dirs = reshape([ &
         1,  0,  0,  &
         -1,  0,  0,  &
         0,  1,  0,  &
         0, -1,  0,  &
         0,  0,  1,  &
         0,  0, -1   &
         ], [6,3])

    ! list active clusters
    num_clusters = count_active_clusters(S)
    allocate(active_clusters(num_clusters))
    cluster_idx = 0
    do i = 0, S%N-1
       if (S%cluster_size(i) > 0) then
          cluster_idx = cluster_idx + 1
          active_clusters(cluster_idx) = i
       end if
    end do

    ! select a random cluster
    call random_number(r)
    cluster_idx = min(int(r * real(num_clusters)), num_clusters-1)
    cid = active_clusters(cluster_idx+1)
    deallocate(active_clusters)

    ! choose the direction randomly
    call random_number(r)
    d = merge(-1,1,r<0.5)
    dx = 0; dy = 0; dz = 0
    call random_number(r)
    select case (int(r*3))
    case (0); dx = d
    case (1); dy = d
    case (2); dz = d
    end select

    ! can it move?
    can_move = .true.
    do i = 0, S%N-1
       if (S%particles(i)%cluster_id == cid) then
          nx = S%particles(i)%x + dx
          ny = S%particles(i)%y + dy
          nz = S%particles(i)%z + dz
          if (nx<0 .or. nx>=L .or. ny<0 .or. ny>=L .or. nz<0 .or. nz>=L) then
             can_move = .false.; exit
          end if
          if (S%grid(nx,ny,nz) /= -1) then
             if (S%particles(S%grid(nx,ny,nz))%cluster_id /= cid) then
                can_move = .false.; exit
             end if
          end if
       end if
    end do
    if (.not. can_move) return

    ! MC acceptance
    prob = real(S%cluster_size(cid))**(-alpha)
    call random_number(r)
    if (r < prob) then
       S%hasAdvanced = .true.
    else
       return
    end if

    ! detect if some cluster will merge
    ! in the better case, it avoids doing the flood fill which
    ! is computationally costly
    hasAggregated = .false.
    do i = 0, S%N-1
       if (S%particles(i)%cluster_id == cid) then
          do j = 1,6
             nx = S%particles(i)%x + dx + dirs(j,1)
             ny = S%particles(i)%y + dy + dirs(j,2)
             nz = S%particles(i)%z + dz + dirs(j,3)
             if (nx<0 .or. nx>=L .or. ny<0 .or. ny>=L .or. nz<0 .or. nz>=L) cycle
             neighbor_id = S%grid(nx,ny,nz)
             if (neighbor_id /= -1 .and. S%particles(neighbor_id)%cluster_id /= cid) then
                hasAggregated = .true.; exit
             end if
          end do
          if (hasAggregated) exit
       end if
    end do

    ! count members of the cluster
    members_count = 0
    do i = 0, S%N-1
       if (S%particles(i)%cluster_id == cid) members_count = members_count + 1
    end do
    allocate(members(members_count))

    j = 0
    do i = 0, S%N-1
       if (S%particles(i)%cluster_id == cid) then
          j = j + 1
          members(j) = i
       end if
    end do

    ! free the grid befove the move
    do j = 1, members_count
       i = members(j)
       x_old = S%particles(i)%x
       y_old = S%particles(i)%y
       z_old = S%particles(i)%z
       S%grid(x_old,y_old,z_old) = -1
    end do

    ! move 
    do j = 1, members_count
       i = members(j)
       S%particles(i)%x = S%particles(i)%x + dx
       S%particles(i)%y = S%particles(i)%y + dy
       S%particles(i)%z = S%particles(i)%z + dz
    end do

    ! reinject in the grid
    do j = 1, members_count
       i = members(j)
       x_old = S%particles(i)%x
       y_old = S%particles(i)%y
       z_old = S%particles(i)%z
       S%grid(x_old,y_old,z_old) = i
    end do
    deallocate(members)

    if(hasAggregated) then
       call flood_fill_all(S)
       call update_cluster_sizes(S)
    end if
    
    call checkEndSimulation(S)
  end subroutine trial
end module subroutines
