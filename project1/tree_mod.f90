module tree_mod
    use node_mod
    implicit none

contains
    ! logical function is_legal(transaction)
    ! subroutine tree_destroy()
    ! subroutine tree_init()
    ! subroutine tree_print()

    logical function is_legal(transaction)
        character(len=*), intent(in) :: transaction

        if (transaction == "borrowed" .or. transaction == "lent") then
            is_legal = .true.
        else
            is_legal = .false.
        end if
        return
    end function is_legal

    subroutine tree_destroy()
        if(associated(root)) call node_destroy(root)        
    end subroutine tree_destroy

    subroutine tree_init()
        character(len=20)     :: deity1, deity2, transaction, preposition
        real(kind=8)          :: amount
        integer               :: iostatus
        type(a_node), pointer :: d1_ptr, d2_ptr

        do
            !read and parse each line
            read(*, *, iostat=iostatus) deity1, transaction, amount, preposition, deity2
            if (iostatus /= 0) exit

            ! if transaction is illegal print warning and skip to next line
            if (is_legal(transaction) .eqv. .false.) then
                print "(3a)", "Wrong transaction '", trim(transaction), "'."
                cycle
            end if

            call node_insert(root, trim(deity1), d1_ptr)
            call node_insert(root, trim(deity2), d2_ptr)

            ! if "borrowed" then deity1%credit => deity2 and deity2%debit => deity1
            ! if "lent" it's the other way around
            if (transaction == "borrowed") then
                call node_add_transaction(d1_ptr, d2_ptr, amount)
            else
                call node_add_transaction(d2_ptr, d1_ptr, amount)
            end if
        end do

    end subroutine tree_init

    subroutine tree_print()
        real(kind=8) :: debit, credit
        print "(a)", ""
        if (associated(root)) then
            call node_print(root)
            call node_total_amounts(root, debit, credit)
        end if
        print "(a)", ""
        print "(a, f0.2)", "Net debit:  ", debit
        print "(a, f0.2)", "Net credit: ", debit
    end subroutine
end module tree_mod