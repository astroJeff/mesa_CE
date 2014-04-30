            
            integer function how_many_extra_binary_history_columns(b)
               use binary_def, only: binary_info
               type (binary_info), pointer :: b
            end function how_many_extra_binary_history_columns
            
            subroutine data_for_extra_binary_history_columns(b, n, names, vals, ierr)
               use const_def, only: dp
               use binary_def, only: maxlen_binary_history_column_name, binary_info
               type (binary_info), pointer :: b
               integer, intent(in) :: n
               character (len=maxlen_binary_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_binary_history_columns
