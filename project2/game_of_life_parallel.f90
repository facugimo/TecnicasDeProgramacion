program game_of_life
    use mpi_f08

    implicit none
    integer :: height, width
    integer :: max_gen, gen
    logical, dimension(:, :), pointer :: old_world, new_world, tmp_world

    ! New declarations:
    type(MPI_Comm) :: comm
    integer :: n_ranks, my_rank, root
    integer :: n_cols, n_rows, my_col, my_row
    integer :: my_col_beg, my_row_beg, my_col_end, my_row_end
    type(MPI_Datatype) :: a_row, a_col, a_tmp_row
    integer(kind=MPI_ADDRESS_KIND) :: low_bound, real_extent
    type(MPI_Status) :: status



    ! We start by initializing the parallel environment (and the communicator).
    call MPI_Init()
    comm = MPI_COMM_WORLD

    ! We establish the number of processes and the ID of the current process.
    call MPI_Comm_size(comm, n_ranks)
    call MPI_Comm_rank(comm, my_rank)

    ! The root process will be the one responsible for all input/output.
    root = 0

    if (my_rank == root) then
        ! First we read the first line which will give us info about the size
        ! of the board, the maximum number of generations allowed and the
        ! number of rows and columns in our torus.
        read *, height, width, max_gen, n_rows, n_cols

        ! We check that the number of processes is correct for the required
        ! geometry and abort otherwise.
        if (n_rows * n_cols /= n_ranks) then
            write (*,*) "Please, run with -np", n_rows * n_cols
            call MPI_Abort(comm, MPI_ERR_TOPOLOGY)
        end if
    end if

    ! The variables we just read need to be broadcasted to all processes frem
    ! root.
    call MPI_Bcast(height,  1, MPI_INTEGER, root, comm)
    call MPI_Bcast(width,   1, MPI_INTEGER, root, comm)
    call MPI_Bcast(max_gen, 1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_rows,  1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_cols,  1, MPI_INTEGER, root, comm)

    ! And we need to know which row and column of the torus we are in. The
    ! first process (root) takes the cell at (0, 0) and from there we advance
    ! assigning processes to the next rows of the first column until it is
    ! finished and then we go to the next column. Repeat for the whole torus.
    my_row = modulo(my_rank, n_rows)
    my_col = (my_rank - my_row) / n_rows

    ! The height and width of the board will be divided between the columns and
    ! rows of the torus. We need to know at which height and width each torus
    ! cell (process) begins and ends.

    call get_beg_end(my_row, n_rows, height, my_row_beg, my_row_end)
    call get_beg_end(my_col, n_cols, width,  my_col_beg, my_col_end)

    ! We now create a new datatype called a_col for each col in the current
    ! process (re - rb + 1). We do the same for the rows. Cols are contiguous
    ! in memory, so we can use MPI_Type_contiguous
    call MPI_Type_contiguous(my_row_end-my_row_beg+1, MPI_LOGICAL, a_col)
    call MPI_Type_commit(a_col)

    ! But rows are not, so for this datatype we need to use MPI_Type_vector.
    call MPI_Type_get_extent(MPI_LOGICAL, low_bound, real_extent)
    call MPI_Type_vector(my_col_end - my_col_beg + 1, 1, &
                         my_row_end - my_row_beg + 3, MPI_LOGICAL, a_tmp_row)
    call MPI_Type_create_resized(a_tmp_row, low_bound, real_extent, a_row)
    call MPI_Type_commit(a_row)

    ! Finally we allocate the boards old_world and new_world. Each process
    ! allocates two extra rows and columns which will be updated from the
    ! neighbouring ranks of the torus.
    allocate(old_world(my_row_beg-1 : my_row_end+1, my_col_beg-1 : my_col_end+1))
    allocate(new_world(my_row_beg-1 : my_row_end+1, my_col_beg-1 : my_col_end+1))

    call read_map( old_world, width )
    call update_borders( old_world )

    do gen = 1, max_gen
        if (my_rank == 0) then 
            print "(a, i0)", "Generation ", gen
        end if
        call print_map( old_world, width )
        call next_gen( old_world, new_world)
        call update_borders( new_world )
        ! call wait_cls( 100 )
        if (world_is_still( old_world, new_world )) exit
        ! Swap maps
        tmp_world => old_world;  old_world => new_world;  new_world => tmp_world

    end do

    ! The two derived data types have to be released.
    call MPI_Type_free(a_col)
    call MPI_Type_free(a_row)

    if (associated( old_world )) deallocate(old_world)
    if (associated( new_world )) deallocate(new_world)

    ! And we finalize MPI to end our program.
    call MPI_Finalize()

contains

    logical function world_is_still( old_map, new_map )
        logical, dimension(:, :), pointer, intent(in) :: old_map, new_map
        logical :: still

        still = all( old_map .eqv. new_map )
        call MPI_Allreduce(still, world_is_still, 1, MPI_LOGICAL, MPI_LAND, comm)
    end function world_is_still

    subroutine update_borders( map)
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer :: upper, lower, left, right
        integer :: upper_left, upper_right, lower_left, lower_right
        
        ! For the rows we need to know the rank of the cell above and below.
        upper = get_rank(my_row - 1, my_col)
        lower = get_rank(my_row + 1, my_col)
        
        ! Inner rows
        ! Send first row to upper cell.
        call MPI_Sendrecv(map(my_row_beg    , my_col_beg), 1, a_row, upper, 1, &
                          map(my_row_end + 1, my_col_beg), 1, a_row, lower, 1, comm, status)
        ! Send last row to lower cell.
        call MPI_Sendrecv(map(my_row_end    , my_col_beg), 1, a_row, lower, 2, &
                          map(my_row_beg - 1, my_col_beg), 1, a_row, upper, 2, comm, status)


        ! inner columns
        left  = get_rank(my_row, my_col - 1)
        right = get_rank(my_row, my_col + 1)
             
        ! Send first column to left cell.
        call MPI_Sendrecv(map(my_row_beg, my_col_beg    ), 1, a_col, left,  3, &
                          map(my_row_beg, my_col_end + 1), 1, a_col, right, 3, comm, status)
        ! Send last column to right cell.
        call MPI_Sendrecv(map(my_row_beg, my_col_end    ), 1, a_col, right, 4, &
                          map(my_row_beg, my_col_beg - 1), 1, a_col, left , 4, comm, status)

        ! Corners
        upper_left  = get_rank(my_row - 1, my_col - 1)
        upper_right = get_rank(my_row - 1, my_col + 1)
        lower_left  = get_rank(my_row + 1, my_col - 1)
        lower_right = get_rank(my_row + 1, my_col + 1)

        call MPI_Sendrecv(map(my_row_beg    , my_col_beg    ), 1, MPI_LOGICAL, upper_left , 4, &
                          map(my_row_end + 1, my_col_end + 1), 1, MPI_LOGICAL, lower_right , 4, comm, status)

        call MPI_Sendrecv(map(my_row_beg    , my_col_end    ), 1, MPI_LOGICAL, upper_right, 4, &
                          map(my_row_end + 1, my_col_beg - 1), 1, MPI_LOGICAL, lower_left , 4, comm, status)

        call MPI_Sendrecv(map(my_row_end    , my_col_beg    ), 1, MPI_LOGICAL, lower_left , 4, &
                          map(my_row_beg - 1, my_col_end + 1), 1, MPI_LOGICAL, upper_right , 4, comm, status)
        
        call MPI_Sendrecv(map(my_row_end    , my_col_end    ), 1, MPI_LOGICAL, lower_right, 4, &
                          map(my_row_beg - 1, my_col_beg - 1), 1, MPI_LOGICAL, upper_left , 4, comm, status)

    end subroutine update_borders

    subroutine read_map( map, w )
        logical, dimension(:, :), pointer, intent(inout) :: map
        integer, intent(in) :: w
        integer :: line
        
        
        ! Only root will actually read the map and then send the data to each
        ! process.
        if (my_rank == root) then
            block
                integer :: i
                integer :: col, row, col_beg, row_beg, col_end, row_end
                integer :: destiny
                character(len=:),      allocatable :: line_char
                logical, dimension(:), allocatable :: line_log

                allocate(character(len=w) :: line_char)
                allocate(line_log( 1:w))

                do row = 0, n_rows - 1
                    ! First we need to know how many lines the current row has.
                    call get_beg_end(row, n_rows, height, row_beg, row_end)
                    
                    ! We read each line and convert it to logical.
                    do line = row_beg, row_end
                        read *, line_char
                        do i = 1, w
                            line_log(i:i) = line_char(i:i) == "X"
                        end do

                        ! Now we cycle through the columns to send the data to the
                        ! corresponding process.
                        do col = 0, n_cols - 1
                            ! Now we need to know where each column begins and ends in the
                            ! line so we can send that part of the line to the correct rank.
                            call get_beg_end(col, n_cols, width, col_beg, col_end)
                            destiny = get_rank(row, col)
                            if (destiny == root) then
                                map(line, col_beg:col_end) = line_log(col_beg:col_end)
                            else
                                call MPI_Send(line_log(col_beg), col_end-col_beg+1, &
                                              MPI_LOGICAL, destiny, 0, comm)
                            end if
                        end do
                    end do
                end do

                ! Free memory
                if (allocated(line_char)) deallocate(line_char)
                if (allocated(line_log))  deallocate(line_log)
            end block

        ! Other processes receive the data
        else
            do line = my_row_beg, my_row_end
                call MPI_Recv(map(line, my_col_beg), 1, a_row, root, 0, comm, status)
            end do
        end if
    end subroutine read_map

    subroutine print_map( map, w )
        logical, dimension(:, :), pointer, intent(in) :: map
        integer, intent(in) :: w
        integer :: line
        
        ! Only root will actually print the map so first each process needs to
        ! send its portion of the map to root.
        ! We will send one line at the time.
        if (my_rank /= root) then
            do line = my_row_beg, my_row_end
                call MPI_Send(map(line, my_col_beg), 1, a_row, &
                root, 0, comm)
            end do
        else ! root needs to receive each portion of the map and print it.
            block
                integer :: i
                integer :: row, col, row_beg, col_beg, row_end, col_end
                integer :: source
                logical, dimension(:), allocatable :: line_log
                character(len=:),      allocatable :: line_char
                
                allocate(character(len=w) :: line_char)
                allocate(line_log(1:w))

                ! Again, we go line by line to receive the data from the
                ! different processes. we will go through all lines each row
                ! has before passing to the next row
                do row = 0, n_rows - 1
                    ! We need to know how many lines the current row has.
                    call get_beg_end(row, n_rows, height, row_beg, row_end)
                    ! Now we can cycle through all lines in the row.
                    do line = row_beg, row_end
                        ! We cycle through all columns of the torus.
                        do col = 0, n_cols - 1
                            ! We need to know which rank is sending the data
                            ! from the current coordinates (row, col).
                            source = get_rank(row, col)
                            ! And look where in the line each column begins
                            ! and ends.
                            call get_beg_end(col, n_cols, width, col_beg, col_end)
                            
                            ! We don't receive from root as it's the current
                            ! process.
                            if (source == root) then
                                line_log(col_beg:col_end) = map(line, col_beg:col_end)
                            else ! We do receive from other processes.
                                call MPI_Recv(line_log(col_beg), col_end-col_beg+1, &
                                              MPI_LOGICAL, source, 0, comm, status)
                            end if
                        end do

                        ! We convert the line to characters using 
                        ! (true=>'X' and false=>'.') and print it.
                        do i=1, w
                            line_char(i:i) = merge('X', '.', line_log(i))
                        end do
                        print "(a)",  line_char
                    end do
                end do

                ! We print an empty line at the end.
                print *
                ! And free the memory.
                if (allocated(line_log))  deallocate(line_log)
                if (allocated(line_char)) deallocate(line_char)
                end block
            end if

    end subroutine print_map

    subroutine next_gen( old_map, new_map)
        logical, dimension(:, :), pointer, intent(inout) :: old_map, new_map
        integer :: i, j
        integer :: c ! the number of live neighbors

        do j = my_col_beg, my_col_end
            do i = my_row_beg, my_row_end
                c = count( old_map(i - 1:i + 1, j - 1:j + 1) )
                if (old_map(i, j)) then ! cell is live
                    new_map(i, j) = merge( .true., .false., 3 <= c .and. c <= 4 )
                else ! cell is dead
                    new_map(i, j) = merge( .true., .false., c == 3 )
                end if
            end do
        end do
    end subroutine next_gen

    ! Wait specified number of ms and then clear the terminal screen.
    subroutine wait_cls( ms )
        integer, intent(in) :: ms
        integer :: tick, tack
        real :: rate

        call system_clock( count=tick, count_rate=rate )
        do
            call system_clock( count=tack )
            if (real( tack - tick ) / rate >= ms * 1e-3) exit
        end do
        ! Clear the terminal screen using console escape code ^[2J.
        print "(2a)", achar( 27 ), '[2J'
    end subroutine wait_cls

    ! We obtain the begining (b) and ending (e) indices for the torus cell
    ! along the height/width.
    subroutine get_beg_end(cell, n_cells, size, begin, end)
        integer, intent(in)    :: cell, n_cells, size
        integer, intent(inout) :: begin, end
        integer :: remainder, quotient

        remainder = modulo(size, n_cells)
        quotient  = (size - remainder) / n_cells
        begin = 1 + quotient * (cell    ) + min(remainder, cell    )
        end   =     quotient * (cell + 1) + min(remainder, cell + 1)
        return
    end subroutine get_beg_end


    ! We obtain the rank of a process from its coordinates.
    integer function get_rank(row, col)
        integer, intent(in) :: row, col
        integer :: r, c

        r = row
        c = col

        ! Boundary conditions:
        if (row < 0 .or. row == n_rows) then
            r = n_rows - abs(row)
        end if
        if (col < 0 .or. col == n_cols) then 
            c = n_cols - abs(col)
        end if 

        get_rank = r + c * n_rows
        return
    end function get_rank


end program game_of_life