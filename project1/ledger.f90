program ledger
    use tree_mod
    implicit none
    
    ! read data from input, parse lines and costruct the tree
    call tree_init()
    ! print the tree in order, as well as the total debit and credit
    call tree_print()
    ! destroy the tree so that there are no memory leaks
    call tree_destroy()

end program ledger