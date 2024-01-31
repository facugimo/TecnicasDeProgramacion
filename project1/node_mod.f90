module node_mod
    implicit none
    
    type :: a_node
        character(len=:), allocatable :: name
        type(a_list_item), pointer    :: debit
        type(a_list_item), pointer    :: credit
        type(a_node), pointer         :: left
        type(a_node), pointer         :: right
    end type a_node

    type :: a_list_item
        type(a_node), pointer      :: deity
        real(kind=8)               :: amount
        type(a_list_item), pointer :: next
    end type a_list_item


    type(a_node), pointer :: root => null()

contains
    ! First node_* subroutines and functions
    ! then list_* subroutines and functions


    ! subroutine node_add_transaction(deity1, deity2, amount)
    ! recursive subroutine node_destroy(current)
    ! recursive subroutine node_insert(current,  deity, deity_ptr)
    ! recursive subroutine node_print(current)
    ! recursive subroutine node_total_amounts(current, debit, credit)

    subroutine node_add_transaction(deity1, deity2, amount)
        type(a_node), pointer, intent(in out) :: deity1
        type(a_node), pointer, intent(in out) :: deity2
        real(kind=8),          intent(in)     :: amount
        
        ! the transaction consists in adding deity1%credit => deity2 and deity2%debit => deity1
        call list_insert(deity2, amount, deity1%credit)
        call list_insert(deity1, amount, deity2%debit)
    end subroutine node_add_transaction

    recursive subroutine node_destroy(current)
        type(a_node), pointer, intent(in out) :: current

        ! first we destroy nodes to the left and right
        if (associated(current%left)) call node_destroy(current%left)
        if (associated(current%right)) call node_destroy(current%right)
        ! then we destroy the debit/credit lists
        if (associated(current%debit)) call list_destroy(current%debit)
        if (associated(current%credit)) call list_destroy(current%credit)

        ! and finally we destroy the node itself
        deallocate(current)
    end subroutine node_destroy

    recursive subroutine node_insert(current,  deity, deity_ptr)
        type(a_node), pointer, intent(in out)     :: current
        character(len=*), intent(in) :: deity
        type(a_node), pointer, intent(in out)     :: deity_ptr

        if (.not. associated(current)) then
            ! if the node is empty we create it now
            allocate(current)
            current = a_node(name=deity, debit=null(), credit=null(), left=null(), right=null())
            print "(2a)", deity, " has been added."
            deity_ptr => current
        ! if the node is a deity we check in which side we need to insert
        else if (current%name < deity) then
            call node_insert(current%right, deity, deity_ptr)
        else if (current%name > deity) then
            call node_insert(current%left, deity, deity_ptr)
        else
            ! if the deity is in the tree we just point the pointer
            deity_ptr => current
        end if
    end subroutine node_insert

    recursive subroutine node_print(current)
        type(a_node), intent(in out) :: current
        ! we print first all nodes to the left
        if (associated(current%left)) call node_print(current%left)
        ! then the current node
        print "(2a)", current%name, ":"
        if (associated(current%debit)) then
            print "(4x, a)", "debit"
            call list_print(current%debit)
        end if
        if (associated(current%credit)) then
            print "(4x, a)", "credit"
            call list_print(current%credit)
        end if
        ! and finally nodes to the right
        if (associated(current%right)) call node_print(current%right)
    end subroutine node_print

    recursive subroutine node_total_amounts(current, debit, credit)
        type(a_node), pointer, intent(in)     :: current
        real(kind=8),          intent(in out) :: debit
        real(kind=8),          intent(in out) :: credit

        if (associated(current%debit)) then
            debit = debit + list_amount(current%debit)
        end if
        if (associated(current%credit)) then
            credit = credit + list_amount(current%credit)
        end if
        if (associated(current%left)) then
            call node_total_amounts(current%left, debit, credit)
        end if
        if (associated(current%right)) then
            call node_total_amounts(current%right, debit, credit)
        end if
    end subroutine node_total_amounts

    ! recursive function list_amount(list_element) result(amount)
    ! recursive subroutine list_destroy(list_element)
    ! recursive subroutine list_insert(deity, amount,  current)
    ! recursive subroutine list_print(element)

    recursive function list_amount(list_element) result(amount)
        type(a_list_item), pointer, intent(in) :: list_element
        real(kind=8)                                   :: amount

        amount = list_element%amount
        if (associated(list_element%next)) then
            amount = amount + list_amount(list_element%next)
        end if
        return
    end function list_amount
    
    recursive subroutine list_destroy(list_element)
        type(a_list_item), pointer, intent(in out) :: list_element

        ! first we destroy the next element
        if (associated(list_element%next)) call list_destroy(list_element%next)
        ! and then we destroy the current element
        deallocate(list_element)
    end subroutine list_destroy

    recursive subroutine list_insert(deity, amount,  current)
        type(a_node), pointer,      intent(in)     :: deity
        real(kind=8),               intent(in)     :: amount
        type(a_list_item), pointer, intent(in out) ::  current

        ! if the list is empty, we just create it with the corresponding amount
        if (.not. associated(current)) then
            allocate(current)
            current = a_list_item(deity=deity, amount=amount, next=null())
        else
            ! if the list is not empty, we have 3 cases:
            ! we need to modify the current element
            ! we have to insert our element before the current list element
            ! we have to insert (recursively) our element after the current element
            if (associated(deity, current%deity)) then
                current%amount = current%amount + amount
            else if (current%deity%name > deity%name) then
                block
                    type(a_list_item), pointer :: new_element
                    allocate(new_element)
                    new_element = a_list_item(deity=deity, amount=amount, next=current)
                    current => new_element
                end block
            else if (current%deity%name < deity%name) then
                call list_insert(deity, amount, current%next)
            end if
        end if
    end subroutine list_insert

    recursive subroutine list_print(element)
        type(a_list_item), pointer, intent(out) :: element

        print "(8x, a, 1x, f0.2)", element%deity%name, element%amount
        if (associated(element%next)) call list_print(element%next)
    end subroutine list_print


end module node_mod