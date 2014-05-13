! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton and Pablo Marchant
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 

      module binary_lib
      
      
      implicit none


      contains


      subroutine run1_binary(tst, &
            ! star extras
            extras_controls, &
            extras_startup, &
            extras_check_model, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns, &
            extras_finish_step, &
            extras_after_evolve, &
            ! binary extras
            extras_binary_controls, &
            how_many_extra_binary_history_columns, &
            data_for_extra_binary_history_columns, &

            ierr)
            
         use run_binary_support, only: do_run1_binary
         use binary_def, only: init_binary_data
         use star_def, only: star_info
         
         logical, intent(in) :: tst
         
         interface

            subroutine extras_controls(s, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(out) :: ierr
            end subroutine extras_controls      
     
            integer function extras_startup(s, id, restart, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id
               logical, intent(in) :: restart
               integer, intent(out) :: ierr
            end function extras_startup
      
            integer function extras_check_model(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function extras_check_model

            integer function how_many_extra_history_columns(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function how_many_extra_history_columns
            
            subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
               use const_def, only: dp
               use star_def, only: star_info, maxlen_history_column_name
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra, n
               character (len=maxlen_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_history_columns
      
            integer function how_many_extra_profile_columns(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function how_many_extra_profile_columns      
      
            subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
               use const_def, only: dp
               use star_def, only: star_info, maxlen_profile_column_name
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra, n, nz
               character (len=maxlen_profile_column_name) :: names(n)
               real(dp) :: vals(nz,n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_profile_columns
      
            integer function extras_finish_step(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
               integer :: ierr
            end function extras_finish_step     
      
            subroutine extras_after_evolve(s, id, id_extra, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
               integer, intent(out) :: ierr
            end subroutine extras_after_evolve

            subroutine extras_binary_controls(ierr)
               use star_def, only: star_info
               integer, intent(out) :: ierr
            end subroutine extras_binary_controls      

            integer function how_many_extra_binary_history_columns(b)
               use binary_def, only: binary_info
               type (binary_info), pointer :: b
            end function how_many_extra_binary_history_columns
            
            subroutine data_for_extra_binary_history_columns(b, s, n, names, vals, ierr)
!            subroutine data_for_extra_binary_history_columns(b, n, names, vals, ierr)
               use const_def, only: dp
               use binary_def, only: maxlen_binary_history_column_name, binary_info
               use star_def, only: star_info
               type (binary_info), pointer :: b
               type (star_info), pointer :: s
               integer, intent(in) :: n
               character (len=maxlen_binary_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_binary_history_columns

         end interface
         
         integer, intent(out) :: ierr
         
         call init_binary_data
         
         call do_run1_binary(tst, &
            ! star extras
            extras_controls, &
            extras_startup, &
            extras_check_model, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns, &
            extras_finish_step, &
            extras_after_evolve, &
            ! binary extras
            extras_binary_controls, &
            how_many_extra_binary_history_columns, &
            data_for_extra_binary_history_columns, &

            ierr)
      
      end subroutine run1_binary
      

      end module binary_lib

