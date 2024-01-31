program levialdi_parallel
    use mpi_f08
    use levialdi_parallel_class
    implicit none
    
    integer                      :: h, w
    type(a_linked_list), target  :: image
    type(a_linked_list), pointer :: current
    integer :: K = 0

    integer :: n_ranks, my_rank, root
    integer :: n_cols, n_rows, my_col, my_row
    integer :: c0, r0, cf, rf
    logical :: cont

    type(MPI_Comm)                 :: comm
    type(MPI_Datatype)             :: a_row, a_col, a_tmp_row
    type(MPI_Datatype)             :: a_row_lbl, a_col_lbl, a_tmp_row_lbl
    integer(kind=MPI_ADDRESS_KIND) :: low_bound, real_extent


    ! Initialize MPI
    call MPI_Init()
    comm = MPI_COMM_WORLD

    ! We establish the number of processes and the ID of the current process.                                                               
    call MPI_Comm_size(comm, n_ranks)
    call MPI_Comm_rank(comm, my_rank)

    ! The root process will be the one responsible for all input/output.                                                                 
    root = 0

    ! All input/output is done by root process.
    if (my_rank == root) then
        ! We start the program by reading the input file.
        ! First we read the first line.
        read *, h, w, n_rows, n_cols

        ! We check that the number of processes is correct for the required
        ! geometry.
        if (n_rows * n_cols /= n_ranks) then
            print *, "Please, run with -np", n_rows * n_cols
            call MPI_Abort(comm, MPI_ERR_TOPOLOGY)
        end if
    end if
    
    ! The variables we just read need to be broadcasted to all processes from
    ! root.
    call MPI_Bcast(h,  1, MPI_INTEGER, root, comm)
    call MPI_Bcast(w,   1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_rows,  1, MPI_INTEGER, root, comm)
    call MPI_Bcast(n_cols,  1, MPI_INTEGER, root, comm)

    ! And we need to know which row and column we are in. The first process
    ! (root) takes the cell at (0, 0) and from there we advance assigning
    ! processes to the next rows of the first column until it is finished and
    ! then we go to the next column.
    my_row = modulo(my_rank, n_rows)
    my_col = (my_rank - my_row) / n_rows

    ! The height and width of the image will be divided between in columns and
    ! rows. We need to know at which height and width each cell (process)
    ! begins and ends.

    call get_beg_end(my_row, n_rows, h, r0, rf)
    call get_beg_end(my_col, n_cols, w,  c0, cf)

    ! We now create a new datatype called a_col for each col in the current                                                                 
    ! process (rf - rb + 1). Cols are contiguous in memory, so we can use
    ! MPI_Type_contiguous.
    call MPI_Type_contiguous(rf-r0+1, MPI_LOGICAL, a_col)
    call MPI_Type_commit(a_col)

    ! But rows are not, so for this datatype we need to use MPI_Type_vector.
    ! And rows are bigger as they include ghost cells at the ends.
    call MPI_Type_get_extent(MPI_LOGICAL, low_bound, real_extent)
    call MPI_Type_vector(cf - c0 + 3, 1, &
                         rf - r0 + 3, MPI_LOGICAL, a_tmp_row)
    call MPI_Type_create_resized(a_tmp_row, low_bound, real_extent, a_row)
    call MPI_Type_commit(a_row)
    
    ! Now we do the same but for the label rows and columns. The only
    ! difference is these are made of integers instead of logicals
    call MPI_Type_contiguous(rf-r0+1, MPI_INTEGER, a_col_lbl)
    call MPI_Type_commit(a_col_lbl)

    call MPI_Type_get_extent(MPI_INTEGER, low_bound, real_extent)
    call MPI_Type_vector(cf - c0 + 3, 1, &
                         rf - r0 + 3, MPI_INTEGER, a_tmp_row_lbl)
    call MPI_Type_create_resized(a_tmp_row_lbl, low_bound, real_extent, a_row_lbl)
    call MPI_Type_commit(a_row_lbl)

    ! We can now start by initializing the image matrices and read the data.
    call image%init(r0, rf, c0, cf, my_rank)
    call image%read(root, n_rows, n_cols, w, h, comm, a_col, a_row)
    
    ! With this data we can allocate the matrix for the image and read it.
    current => image

    ! Call to all reduce in order to continue if at least one of the processes
    ! has not finished shrinking the image.
    cont = ANY(current%image%img)
    call MPI_Allreduce(cont, cont, 1, MPI_LOGICAL, MPI_LOR, comm)

    do while (cont)
        call current%shrink()
        current => current%next
        ! Now we update the ghost cells
        call current%upd(n_rows, n_cols, .true., comm, a_col, a_row)
    
        ! Call to all reduce in order to continue if at least one of the processes
        ! has not finished shrinking the image.
        cont = ANY(current%image%img)
        call MPI_Allreduce(cont, cont, 1, MPI_LOGICAL, MPI_LOR, comm)    
    end do


    do while(associated(current))
        call current%grow(K, comm)
        call MPI_Barrier(comm)
        ! Now we update the ghost cells
        call current%upd(n_rows, n_cols, .false., comm, a_col_lbl, a_row_lbl)
 
        current => current%previous
    end do

    call MPI_Allreduce(k, k, 1, MPI_INTEGER, MPI_MAX, comm)

    if (my_rank == root) then
        print "(I0)", K
    end if

    call image%print(k, n_rows, n_cols, w, h, root, a_row_lbl, comm)

    current => image

    ! Destroy recursively all links of the linked list.
    call current%destroy()

    call MPI_Type_free(a_row)
    call MPI_Type_free(a_col)
    call MPI_Type_free(a_tmp_row)
    call MPI_Type_free(a_row_lbl)
    call MPI_Type_free(a_col_lbl)
    call MPI_Type_free(a_tmp_row_lbl)

    call MPI_Finalize()


    
end program levialdi_parallel
