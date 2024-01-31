module levialdi_parallel_class
    use mpi_f08
    implicit none

    public
    ! Each image object contains its row and col initial and final values as
    ! well as two arrays, one for the image proper and the other for the labels
    ! we will assign to each pixel.
    type, public :: an_image
        integer, private :: r0, rf, c0, cf
        logical, dimension(:,:), allocatable :: img
        integer, dimension(:,:), allocatable :: lbl
    end type an_image

    ! Each element of the linked list contains a pointer to an image object and
    ! to the next and previous elements of the list.
    type, public :: a_linked_list
        type(an_image),      pointer :: image
        integer                      :: rank
    
        type(a_linked_list), pointer :: previous
        type(a_linked_list), pointer :: next

        contains
            procedure, public :: init     => levialdi_init
            procedure, public :: read     => levialdi_read
            procedure, public :: upd      => levialdi_update_borders
            procedure, public :: shrink   => levialdi_shrink
            procedure, public :: grow     => levialdi_grow
            procedure, public :: print    => levialdi_print
            procedure, public :: destroy  => levialdi_destroy
    end type a_linked_list

    contains
        ! Initialize an element of the linked list. Sets height and
        ! width and allocates space for the image and label.
        subroutine levialdi_init(self, r0, rf, c0, cf, rank)
            class(a_linked_list), intent(in out) :: self
            integer, intent(in) :: r0, rf, c0, cf, rank
        
            allocate(self%image)
            allocate(self%image%img(r0-1:rf+1, c0-1:cf+1))
            allocate(self%image%lbl(r0-1:rf+1, c0-1:cf+1))
            self%rank = rank

            self%image%r0 = r0
            self%image%rf = rf
            self%image%c0 = c0
            self%image%cf = cf
            self%image%lbl = 0
            self%previous => null()
            self%next     => null()
        end subroutine levialdi_init

        ! Read the raw image data to image%img.
        subroutine levialdi_read(self, root, n_rows, n_cols, w, h, comm, a_col, a_row)
            class(a_linked_list), intent(in out) :: self
            integer,              intent(in)     :: root, n_rows, n_cols, w, h
            type(MPI_Datatype),   intent(in)     :: a_col, a_row
            type(MPI_Comm),       intent(in)     :: comm
            
            integer            :: j
            type(MPI_Status)   :: status
            
            
            if (self%rank == root) then
                block
                    integer :: i, destiny
                    integer :: r0, rf, row
                    integer :: c0, cf, col
                    character(len=:),      allocatable :: line_char
                    logical, dimension(:), allocatable :: line_log

                    allocate(character(w) :: line_char)
                    allocate(line_log(0:w+1))
                    line_log(0)   = .false.
                    line_log(w+1) = .false.

                    do row=0, n_rows - 1
                        ! We need to know where the current row begins and ends.
                        call get_beg_end(row, n_rows, h, r0, rf)
                        do j = r0, rf
                            read *, line_char
                            do i=1, w
                                line_log(i:i) = line_char(i:i) == "x"
                            end do

                            ! Now we cycle through the columns to send the data
                            ! to the corresponding process.
                            do col = 0, n_cols - 1
                                ! Now we need to know where each column begins
                                ! and ends in the line so we can send that part
                                ! of the line to the correct rank.
                                call get_beg_end(col, n_cols, w, c0, cf)
                                destiny = get_rank(row, col, n_rows, n_cols)
                                if (destiny == root) then
                                    self%image%img(j, c0:cf) = line_log(c0:cf)
                                else
                                    call MPI_Send(line_log(c0-1), cf-c0+3, MPI_LOGICAL, &
                                                  destiny, j, comm)
                                end if
                            end do
                        end do
                    end do
                    ! Free memory
                    if (allocated(line_char)) deallocate(line_char)
                    if (allocated(line_log))  deallocate(line_log)
                end block
            ! non-root processes receive the data.
            else
                do j = self%image%r0, self%image%rf
                    call MPI_Recv(self%image%img(j, self%image%c0-1), 1, a_row, &
                                  root, j, comm, status)
                end do
            end if

            ! And finally we assign the ghost cells.
            call self%upd(n_rows, n_cols, .true., comm, a_col, a_row)
        end subroutine levialdi_read

        ! Apply boundary conditions to the img array. In the serial
        ! version sets all borders to false.
        ! If data is true, update the pixels, otherwise update the labels.
        subroutine levialdi_update_borders(self, n_rows, n_cols, data, &
            comm, a_col, a_row)
            class(a_linked_list), intent(in out) :: self
            integer,              intent(in)     :: n_rows, n_cols
            logical,              intent(in)     :: data
            type(MPI_Datatype),   intent(in)     :: a_row, a_col
            type(MPI_Comm),       intent(in)     :: comm


            integer :: row, col, upper, lower, left, right
            type(MPI_Status)   :: status

                        
            ! Calculate the current row and column
            row = modulo(self%rank, n_rows)
            col = (self%rank - row)/n_rows

            ! We need to know which rank is to the left and to the right. For
            ! this we will make use of a toroidal geometry and then reset the
            ! outside borders to 0 each time
            upper = get_rank(row - 1, col,     n_rows, n_cols)
            lower = get_rank(row + 1, col,     n_rows, n_cols)
            left  = get_rank(row,     col - 1, n_rows, n_cols)
            right = get_rank(row,     col + 1, n_rows, n_cols)
    

            if (data) then
                ! Send first column to left cell.
                call MPI_Sendrecv(self%image%img(self%image%r0, self%image%c0), 1, &
                                a_col, left,  1, &
                                self%image%img(self%image%r0, self%image%cf + 1), 1, &
                                a_col, right, 1, comm, status)
                ! Send last column to right cell.
                call MPI_Sendrecv(self%image%img(self%image%r0, self%image%cf), 1, &
                                a_col, right, 2, &
                                self%image%img(self%image%r0, self%image%c0 - 1), 1, &
                                a_col, left , 2, comm, status)
                ! Send first row to upper cell.
                call MPI_Sendrecv(self%image%img(self%image%r0, self%image%c0-1), 1, &
                                a_row, upper, 3, &
                                self%image%img(self%image%rf + 1, self%image%c0-1), 1, &
                                a_row, lower, 3, comm, status)
                ! Send last row to lower cell.
                call MPI_Sendrecv(self%image%img(self%image%rf, self%image%c0-1), 1, &
                                a_row, lower, 4, &
                                self%image%img(self%image%r0 - 1, self%image%c0-1), 1, &
                                a_row, upper, 4, comm, status)
        
                ! All the pixels outside the image are black.
                if (col == 0)  self%image%img(:, self%image%c0 - 1) = .false.
                if (col == n_cols-1) self%image%img(:, self%image%cf + 1) = .false.
                if (row == 0)  self%image%img(self%image%r0 - 1, :) = .false.
                if (row == n_rows-1) self%image%img(self%image%rf + 1, :) = .false.

            else
                ! Send first column to left cell.
                call MPI_Sendrecv(self%image%lbl(self%image%r0, self%image%c0), 1, &
                                a_col, left,  1, &
                                self%image%lbl(self%image%r0, self%image%cf + 1), 1, &
                                a_col, right, 1, comm, status)
                ! Send last column to right cell.
                call MPI_Sendrecv(self%image%lbl(self%image%r0, self%image%cf), 1, &
                                a_col, right, 2, &
                                self%image%lbl(self%image%r0, self%image%c0 - 1), 1, &
                                a_col, left , 2, comm, status)

                ! Send first row to upper cell.
                call MPI_Sendrecv(self%image%lbl(self%image%r0, self%image%c0-1), 1, &
                                a_row, upper, 3, &
                                self%image%lbl(self%image%rf + 1, self%image%c0-1), 1, &
                                a_row, lower, 3, comm, status)
                ! Send last row to lower cell.
                call MPI_Sendrecv(self%image%lbl(self%image%rf, self%image%c0-1), 1, &
                                a_row, lower, 4, &
                                self%image%lbl(self%image%r0 - 1, self%image%c0-1), 1, &
                                a_row, upper, 4, comm, status)
        
                call MPI_Barrier(comm)
                ! All the pixels outside the image have label 0.

                if (col == 0)  self%image%lbl(:, self%image%c0 - 1) = 0
                if (col == n_cols-1) self%image%lbl(:,self%image%cf + 1) = 0
                if (row == 0)  self%image%lbl(self%image%r0 - 1,:) = 0
                if (row == n_rows-1) self%image%lbl(self%image%rf + 1,:) = 0
            end if
                      

        end subroutine levialdi_update_borders

        ! Perform the shrinking part of the levialdi algorithm. 
        subroutine levialdi_shrink(self)
            class(a_linked_list), target, intent(in out) :: self

            type(a_linked_list), pointer :: next
            integer :: i,j

            allocate(next)
            call next%init(self%image%r0, self%image%rf, &
                           self%image%c0, self%image%cf, self%rank)
            self%next => next
            next%previous => self
            next%image%img = self%image%img

            do j=self%image%c0, self%image%cf
                do i=self%image%r0, self%image%rf
                    if ((self%image%img(i, j)         .and. &
                         .not. self%image%img(i, j-1) .and. &
                         .not. self%image%img(i-1, j) .and. &
                         .not. self%image%img(i-1, j-1))    &
                    .or.                          &
                        (.not. self%image%img(i, j)   .and. &
                         self%image%img(i, j-1)       .and. &
                         self%image%img(i-1, j)))           &
                    then
                        next%image%img(i,j) = .not. self%image%img(i,j)
                    end if
                end do
            end do
        end subroutine levialdi_shrink

        subroutine levialdi_grow(self, k, comm)
            class(a_linked_list), intent(in out) :: self
            integer,              intent(in out) :: k
            type(MPI_Comm),       intent(in)     :: comm

            integer :: i, j
        
            ! Case last iteration.
            if (.not. associated(self%next)) then
                ! all pixels are black
                self%image%lbl = 0
                return
            else
                ! Call to allreduce in order to find the
                ! maximum number of k at this point.
                call MPI_Barrier(comm)
                call MPI_Allreduce(k, k, 1, MPI_INTEGER, &
                MPI_MAX, comm)
                do j=self%image%c0, self%image%cf
                    do i=self%image%r0, self%image%rf
                        !!!!!!!!!!!!!!!!!!!!!!
                        ! call MPI_Allreduce(k, k, 1, MPI_INTEGER, &
                        MPI_MAX, comm)
                                if (self%image%img(i,j)) then ! Found a new pixel
                            
                            if (self%next%image%img(i,j)) then
                                self%image%lbl(i, j) = self%next%image%lbl(i,j)
                            end if
                            if (.not. ANY(self%image%img(i-1:i+1, j-1)) .and. &
                                .not. ANY(self%image%img(i-1:i+1, j+1)) .and. &
                                .not. self%image%img(i-1, j) .and.            &
                                .not. self%image%img(i+1, j)                  &
                                .and. self%image%lbl(i, j) == 0               &
                               ) then ! Found a new isolated pixel.    
                                    k = k+1
                                    self%image%lbl(i, j) = k
                            else ! Not isolated.
                                if (self%next%image%img(i+1, j)) then
                                    self%image%lbl(i, j) = self%next%image%lbl(i+1, j)
                                else if (self%next%image%img(i, j+1)) then
                                    self%image%lbl(i, j) = self%next%image%lbl(i, j+1)
                                else if (self%next%image%img(i+1, j+1)) then
                                    self%image%lbl(i, j) = self%next%image%lbl(i+1, j+1)
                                else
                                end if
                            end if

                        else ! Black pixel
                            self%image%lbl(i,j) = 0
                        end if
                    end do
                end do
            end if

        end subroutine levialdi_grow


        subroutine levialdi_print(self, k, n_rows, n_cols, w, h, root, a_row, comm)
            class(a_linked_list), intent(in) :: self
            integer,              intent(in) :: k, n_rows, n_cols, root, w, h
            type(MPI_Comm),       intent(in) :: comm
            type(MPI_Datatype),   intent(in) :: a_row

            integer :: i
            
            ! Only root will actually print the map so first each process needs to
            ! send its portion of the map to root.
            ! We will send one line at the time.
            if (self%rank /= root) then
                do i = self%image%r0, self%image%rf
                    call MPI_Send(self%image%lbl(i, self%image%c0-1), 1, a_row, &
                    root, i, comm)
                end do
            else ! root needs to receive each portion of the map and print it.
                block
                    integer :: row, col, r0, c0, rf, cf, line
                    integer :: source
                    integer,   dimension(:), allocatable :: line_int
                    type(MPI_Status)   :: status

                    allocate(line_int (0:w+1))
                    ! Again, we go line by line to receive the data from the
                    ! different processes. we will go through all lines each row
                    ! has before passing to the next row.
                    do row = 0, n_rows - 1
                        ! We need to know how many lines the current row has.
                        call get_beg_end(row, n_rows, h, r0, rf)
                        ! Now we can cycle through all lines in the row.
                        do line = r0, rf
                            ! We cycle through all columns of the grid
                            do col = 0, n_cols - 1
                                ! We need to know which rank is sending the data
                                ! from the current coordinates (row, col).
                                source = get_rank(row, col, n_rows, n_cols)
                                ! And look where in the line each column begins
                                ! and ends.
                                call get_beg_end(col, n_cols, w, c0, cf)

                                ! We don't receive from root as it's the current
                                ! process.
                                if (source == root) then
                                    line_int(c0:cf) = self%image%lbl(line, c0:cf)
                                else ! We do receive from other processes.
                                    ! We are receiving a row, which has phantom cells at
                                    ! its limits, so we need to receive at (c0-1)
                                    call MPI_Recv(line_int(c0-1), cf-c0+3, &
                                                  MPI_INTEGER, source, line, comm, status)
                                end if
                            end do

                            if (k < 10) then
                                print "(100(1x, i1))",  line_int(1:w)
                            else if (k < 100) then
                                print "(100(1x, i2))",  line_int(1:w)
                            else if (k < 1000) then
                                print "(100(1x, i3))",  line_int(1:w)
                            else
                                print "(100(1x, i4))",  line_int(1:w)
                            end if
                        end do
                    end do

                    ! We print an empty line at the end.
                    print *
                    ! And free the memory.
                    if (allocated(line_int)) deallocate(line_int)
                end block
            end if
        end subroutine levialdi_print

        recursive subroutine levialdi_destroy(self)
            class(a_linked_list), intent(in out) :: self
            
            if (allocated(self%image%img)) deallocate(self%image%img)
            if (allocated(self%image%lbl)) deallocate(self%image%lbl)
            if (associated(self%image))    deallocate(self%image)

            if (associated(self%next)) then
                call self%next%destroy()
                deallocate(self%next)
            end if
        end subroutine levialdi_destroy

        function get_rank(row, col, n_rows, n_cols) result(rank)
            integer, intent(in) :: row, col, n_rows, n_cols
            integer :: rank
        
            rank = modulo( row, n_rows ) + n_rows * modulo( col, n_cols )
        end function get_rank

        ! We obtain the begining (b) and ending (e) indices.
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
            
end module levialdi_parallel_class
