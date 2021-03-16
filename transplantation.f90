! Clonal simulation of transplantation scenarios

! The simulation can be compiled in two different versions depending on the file
! with initial conditions to be used.

! To compile the simulation, see compile.sh for examples how to invoke
! gfortran (required).

! To run the simulation, provide the parameters as command line arguments:
!   transplantation [lambda] [eta] [gamma] [t0] [mu] [r] [L1] [L2] [output time step]
!       [total simulation time] [no. of realisations] [output path (relative)]

program transplantation
  implicit none

	real(kind=8),     parameter                           :: pi2 = 6.283185307179586
	real(kind=8),     parameter                           :: h = 0.
	real(kind=8),     parameter                           :: t_init = 2.

  ! paths
	character(255),                       parameter       :: ref_directory = './results/'
#ifdef gfr
	character(255),                       parameter       :: init_file = './data/init_gfra1.csv'
#endif
#ifdef ngn
	character(255),                       parameter       :: init_file = './data/init_ngn3.csv'
#endif

	integer,                              parameter       :: max_cs = 500

  integer,                              parameter       :: num_syncytia = 1024

  integer,                              parameter       :: ALL = 0
  integer,                              parameter       :: UNITS = 1, CELLS = 2

	! functions
  integer                                               :: rmod, unitstep, seed

  ! parameters
	real(kind=8)                                          :: lambda, eta, gamma, mu, t0, dt, sim_time
  integer                                               :: r_mig, L1, L2, nr, m

	! derived parameters
  integer                                               :: num_rps, num_rps_max, eq_rp, num_sites, num_syncytia_2, rp_min
	real(kind=8)                                          :: gamma_dyn

	! initial conditions
	integer,          dimension(1:2000,1:17)              :: init_cond
	integer,          dimension(1:17)                     :: ic
	integer																								:: io, n_ics

  ! fields
  integer,          dimension(:,:),     allocatable     :: lt
	integer                                               :: prog

	! static arrays
	real(kind=8),     dimension(:),       allocatable     :: larger_zero, larger_one, up_to_half

	! observables
	integer,          dimension(:),       allocatable     :: cmp, c_alive, c_prog_size, clones_sampled
	integer,          dimension(:,:),     allocatable     :: c_size, c_cmp, csd
	real(kind=8),     dimension(:,:,:),   allocatable     :: avg_csd
  real(kind=8),     dimension(:,:),     allocatable     :: avg_size, avg_cmp
	real(kind=8),     dimension(:),       allocatable     :: avg_alive, avg_prog_size

	integer,          dimension(:),       allocatable     :: n_diff
	integer,          dimension(:,:),     allocatable     :: avg_n_diff

  real(kind=8),     dimension(:),       allocatable     :: times

  ! computation variables
  integer																								:: unit_type, n_cells, cs, last_bp, orig_length

  real(kind=8),     dimension(:,:,:),   allocatable     :: prop
	real(kind=8)                                          :: prop_prog, norm
  real(kind=8)                                          :: total_prop, t, time_step, reaction, p_sum, die
  integer																								:: frag_point, temp_lt, temp_id, inv
  logical																								:: searching, temp_pm

	! counters
  integer                                               :: i, a, s, k, r, n, x1, x2, y1, y2, c, uc, rp, nt

	! file names
	character(255)                                        :: directory, directory_r, filename, filename_frame, title

	! command line arguments
	integer,                              parameter       :: req_num_args = 13
  integer                                               :: num_args
  character(255),   dimension(1:req_num_args)           :: cl_args

  real(kind=8)													    			   		:: start_time

	! read command line arguments
  num_args=iargc()
	if (num_args.ne.req_num_args) then
		write(*,'(a)') 'Wrong number of arguments. Terminated.'
		call exit()
	end if
  do i=1,num_args
    call getarg(i, cl_args(i))
  end do

	read(cl_args(1),*) lambda
	read(cl_args(2),*) eta
	read(cl_args(3),*) gamma
	read(cl_args(4),*) t0
	read(cl_args(5),*) mu

	read(cl_args(6),*) r_mig

  read(cl_args(7),*) L1
  read(cl_args(8),*) L2

  read(cl_args(9),*) dt
  read(cl_args(10),*) sim_time
  read(cl_args(11),*) nr

  read(cl_args(12),*) seed

	title = trim(cl_args(13))

	! read initial conditions
	init_cond = 0
	open(unit=1, file=trim(init_file), access='sequential', form='formatted')
		io = 0
		i = 0
		do while (io == 0)
			i = i + 1
			read(1, *, iostat=io) init_cond(i,:)
		end do
		n_ics = i - 1
	close(1)

	! create directory
	directory = trim(ref_directory) // trim(title) // '/'
  call system('mkdir -p ' // trim(directory))

  ! derived parameters
	num_rps_max = nint(sim_time / dt)
	rp_min = floor(t_init / dt)
	num_sites = L1 * L2
	num_syncytia_2 = floor(real(num_syncytia, 8) / 2.)

	! allocate dynamic arrays
	allocate(lt(1:L1,1:L2))
	allocate(prop(1:3,1:L1,1:L2))

	! allocate observables
	allocate(cmp(1:num_syncytia), c_size(rp_min:num_rps_max, 1:2), c_cmp(rp_min:num_rps_max, 1:num_syncytia))
	allocate(c_alive(rp_min:num_rps_max), c_prog_size(rp_min:num_rps_max), clones_sampled(rp_min:num_rps_max))
	allocate(avg_csd(rp_min:num_rps_max, 0:max_cs, 1:2), avg_alive(rp_min:num_rps_max))
	allocate(avg_size(rp_min:num_rps_max, 1:2), avg_cmp(rp_min:num_rps_max, 1:num_syncytia), avg_prog_size(rp_min:num_rps_max))

	allocate(larger_zero(0:num_syncytia), larger_one(0:num_syncytia), up_to_half(0:num_syncytia))

	! initalize static arrays
	larger_zero = 1.
	larger_zero(0) = 0.
	larger_one = 1.
	larger_one(0:1) = 0.
	! note that up_to_half(0) = 0. because of how it is used
	up_to_half = 0.
	up_to_half(1:num_syncytia_2) = 1.

	! initialize random number generator
	call init_random_seed(seed)

	! initialize global observables
	avg_csd = 0.
	avg_size = 0.
	avg_cmp = 0.
	avg_alive = 0.
	avg_prog_size = 0.
	clones_sampled = 0.

	! iterate over realisations
	do r = 1, nr
    write(*,'(a,f5.2,a,i7)',advance='no') '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', &
      100. * real(r) / real(nr), '%    ', r
		num_rps = num_rps_max
		clones_sampled(:num_rps) = clones_sampled(:num_rps) + 1

		! initialise timers and observables
		t = t_init
		rp = -1

		c_cmp = 0
		c_size = 0
		c_alive = 0
		c_alive(:num_rps) = 1
		c_prog_size = 0

		! === physical initial condition ===
		! choose random initial condition
		call random_number(die)
		i = ceiling(die * n_ics)
		ic = init_cond(i,:)

		lt = 0
		prog = 0

		! start site
		x1 = 1
		x2 = 1
		do s = 1, 16
			do i = 1, ic(s)
				searching = .true.
				do while (searching)
					! is lattice site vacant?
					if (lt(x1, x2) == 0) then
						! populate free site with a syncytium of the corresponding length
						lt(x1, x2) = s
						searching = .false.
					else
						! search for adjacent free site
						call random_number(die)
						x1 = rmod(x1 + nint(2. * die - 1.), L1)
						call random_number(die)
						x2 = rmod(x2 + nint(2. * die - 1.), L2)
					end if
				end do
			end do
		end do
		! set progeny
		prog = ic(17)

		! === time evolution ===
		do while (rp <= num_rps)
			! data analysis
			rp = ceiling(t / dt)

			! determine composition of the transplanted clone
			cmp = 0
			do x1 = 1, L1
				do x2 = 1, L2
					s = lt(x1, x2)
					if (s > 0) then
						cmp(s) = cmp(s) + 1
					end if
				end do
			end do
			do k = rp, num_rps
				c_cmp(k, :) = cmp(:)
			end do

			! determine the number of cells
			n_cells = 0.
			do s = 1, num_syncytia
				n_cells = n_cells + s * cmp(s)
			end do
			c_size(rp:num_rps, CELLS) = n_cells

			! determine the number of units
			c_size(rp:num_rps, UNITS) = sum(cmp)

			if (n_cells == 0) then
				c_alive(rp:num_rps_max) = 0
			end if

			! record size of progeny
			c_prog_size(rp:num_rps) = prog

			! === propensities ===
			! time-dependent rates
			if (t < t0) then
				gamma_dyn = gamma
			else
				gamma_dyn = 0.
			end if

			! partial propensities
			prop = 0
			do x1 = 1, L1
				do x2 = 1, L2
					! incomplete division
					prop(1, x1, x2) = up_to_half(lt(x1, x2)) * lambda
					! fragmentation: proportional to the number of bridges
					prop(2, x1, x2) = larger_one(lt(x1, x2)) * (lt(x1, x2) - 1) * eta
					! loss
					prop(3, x1, x2) = larger_zero(lt(x1, x2)) * gamma_dyn
				end do
			end do
			! division of progenitors
			prop_prog = prog * mu

			! === advance time and find reaction ===
      total_prop = sum(prop) + prop_prog
			if (total_prop == 0.) exit

			call random_number(time_step)
			time_step = - (1. / total_prop) * log(time_step)
      t = t + time_step

			call random_number(reaction)
			reaction = total_prop * reaction

			! === reactions ===
			p_sum = 0.
			searching = .true.

			! process 4: progenitor division
			p_sum = p_sum + prop_prog
			if (reaction < p_sum .and. searching) then
				! generate progeny
				prog = prog + 1
				searching = .false.
			end if

			x1 = 0
			do while (x1 < L1 .and. searching)
				x1 = x1 + 1
				x2 = 0
				do while (x2 < L2 .and. searching)
					x2 = x2 + 1

					! process 1: incomplete division
					p_sum = p_sum + prop(1, x1, x2)
					if (reaction < p_sum .and. searching) then
						lt(x1, x2) = lt(x1, x2) * 2

						searching = .false.
						exit
					end if

					! process 2: fragmentation
					p_sum = p_sum + prop(2, x1, x2)
					if (reaction < p_sum .and. searching) then
						! remember the original state of the lattice site
						last_bp = 0
						orig_length = lt(x1, x2)

						! loop through all bridges
						do i = 1, orig_length - 1
							! determine whether bridge breaks
							call random_number(die)
							if (die < 0.5) then
								! truncate original unit
								lt(x1, x2) = orig_length - i

								! choose neighbouring lattice site for invasion
								call random_neighbor(r_mig, x1, x2, L1, L2, y1, y2)

								! track differentiated cells
								prog = prog + lt(y1, y2)

								! invade neighboring site
								lt(y1, y2) = i - last_bp
								! track last breakpoint
								last_bp = i
							end if
						end do

						searching = .false.
						exit
					end if

					! process 3: loss
					p_sum = p_sum + prop(3, x1, x2)
					if (reaction < p_sum .and. searching) then
						! kill
						lt(x1, x2) = 0

						searching = .false.
						exit
					end if

				end do
				if (.not.searching) exit
			end do

		end do

		! === add clone ===
		do unit_type = 1, 2
			do k = rp_min, num_rps
				cs = c_size(k, unit_type)
				if (cs <= max_cs) then
					avg_csd(k, cs, unit_type) = avg_csd(k, cs, unit_type) + 1
				end if
			end do
		end do

		avg_size = avg_size + c_size
		avg_cmp = avg_cmp + c_cmp
		avg_alive = avg_alive + c_alive
		avg_prog_size = avg_prog_size + c_prog_size

	! end realisations
	end do

	! === data analysis ===
	! normalise global quantities
	do k = rp_min, num_rps_max
		if (clones_sampled(k) > 0.) then
			norm = 1. / clones_sampled(k)
		else
			norm = 0.
		end if

		avg_csd(k,:,:) = avg_csd(k,:,:) * norm
		avg_size(k,:) = avg_size(k,:) * norm
		avg_cmp(k,:) = avg_cmp(k,:) * norm
		avg_alive(k) = avg_alive(k) * norm
		avg_prog_size(k) = avg_prog_size(k) * norm
	end do

	! write parameter file
  open(unit=1, file=trim(directory) // 'parameters.csv')
		write(1, '(a20,f17.6)') 'lambda, ',		          lambda
		write(1, '(a20,f17.6)') 'eta, ',			          eta
		write(1, '(a20,f17.6)') 'gamma, ',		          gamma
		write(1, '(a20,f17.6)') 't0, ',		              t0
		write(1, '(a20,f17.6)') 'mu, ',		              mu

		write(1, '(a20,i15)')   'L1, ',			            L1
		write(1, '(a20,i15)')   'L2, ',	                L2

		write(1, '(a20,f17.6)') 'dt, ',			    	      dt
		write(1, '(a20,f17.6)') 'simulation_time, ',	  sim_time
		write(1, '(a20,i15)')   'num_realizations, ',   nr
		write(1, '(a20,i15)')   'seed, ',               seed

		write(1, '(a20,i15)')   'num_syncytia, ',       num_syncytia
  close(1)

	! --- data storage ---

	open(unit=1, file=trim(directory) // 'csd_units.csv')
		write(1,'(a)', advance='no') 'time, '
		do m = 0, max_cs - 1
			write(1,'(i2,a1)', advance='no') m, ','
		end do
		write(1,'(i2)') max_cs

		do k = rp_min, num_rps_max
			write(1,'(f17.6,a1)', advance='no') k * dt, ','
			do m = 0, max_cs - 1
				write(1, '(f15.10,a1)', advance='no') avg_csd(k, m, 1), ','
			end do
			write(1, '(f15.10)') avg_csd(k, max_cs, 1)
		end do
	close(1)

	open(unit=1, file=trim(directory) // 'csd_cells.csv')
		write(1,'(a)', advance='no') 'time, '
		do m = 0, max_cs - 1
			write(1,'(i2,a1)', advance='no') m, ','
		end do
		write(1,'(i2)') max_cs

		do k = rp_min, num_rps_max
			write(1,'(f17.6,a1)', advance='no') k * dt, ','
			do m = 0, max_cs - 1
				write(1, '(f15.10,a1)', advance='no') avg_csd(k, m, 2), ','
			end do
			write(1, '(f15.10)') avg_csd(k, max_cs, 2)
		end do
	close(1)

	open(unit=1, file=trim(directory) // 'cmp.csv')
		write(1,'(a)', advance='no') 'time, '
		do m = 1, num_syncytia - 1
			write(1,'(i2,a1)', advance='no') m, ','
		end do
		write(1,'(i2)') num_syncytia

		do k = rp_min, num_rps_max
			write(1,'(f17.6,a1)', advance='no') k * dt, ','
			do m = 1, num_syncytia - 1
				write(1, '(f15.10,a1)', advance='no') avg_cmp(k, m), ','
			end do
			write(1, '(f15.10)') avg_cmp(k, num_syncytia)
		end do
	close(1)

	open(unit=1,file=trim(directory) // 'size.csv')
		write(1,'(a)') 'time,avg_units,avg_cells,avg_alive,avg_progeny'
		do k = rp_min, num_rps_max
			write(1,'(f17.6,a1,f17.6,a1,f17.6,a1,f17.6,a1,f17.6)') k * dt, ',', avg_size(k, 1), ',', avg_size(k, 2), ',', avg_alive(k), ',', avg_prog_size(k)
		end do
	close(1)

end program transplantation


function rmod(n, m)
	implicit none

	integer                                               :: rmod
	integer,															   intent(in)   :: n, m

	rmod = modulo(n - 1, m) + 1
	return
end function


function unitstep(x)
	implicit none

	integer                                               :: unitstep
	integer                                               :: x

	unitstep = 0
	if (x > 0) unitstep = 1
	return
end function


subroutine random_neighbor(range, x1, x2, L1, L2, n_x1, n_x2)
	implicit none

	integer,															   intent(in)   :: range, x1, x2, L1, L2
	integer,															   intent(out)  :: n_x1, n_x2

	integer                                               :: rmod
	real(kind=8)                                          :: die1, die2

	call random_number(die1)
	call random_number(die2)
	n_x1 = rmod(x1 + (2 * range * nint(die1) - 1), L1)
	n_x2 = rmod(x2 + (2 * range * nint(die2) - 1), L2)
end subroutine


subroutine init_random_seed(m)
  implicit none
  integer, intent(in) :: m
  integer :: i, n
  integer, dimension(:), allocatable :: seed

  call random_seed(size=n)
  allocate(seed(n))

  seed = m + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put=seed)

  deallocate(seed)
end subroutine
