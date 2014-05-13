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
 
      module run_binary_support

      use star_lib
      use star_def
      use const_def
      use utils_lib
      use binary_def
      use binary_private_def
      use binary_ctrls_io, only: do_one_binary_setup

      
      implicit none

      character (len=256), dimension(2) :: inlist_names
      integer :: which_for_pgstar
      logical :: evolve_both_stars
      logical :: trace_binary_rlo
      real(dp) :: stopping_age

      ! model saving info
      integer, dimension(2) :: save_model_numbers
      logical, dimension(2) :: save_model_when_terminates
      character (len=256), dimension(2) :: save_model_filenames

      logical :: read_extra_binary_job_inlist1
      character (len=256) :: extra_binary_job_inlist1_name 
   
      logical :: read_extra_binary_job_inlist2
      character (len=256) :: extra_binary_job_inlist2_name 
   
      logical :: read_extra_binary_job_inlist3
      character (len=256) :: extra_binary_job_inlist3_name 
   
      logical :: read_extra_binary_job_inlist4
      character (len=256) :: extra_binary_job_inlist4_name 
   
      logical :: read_extra_binary_job_inlist5
      character (len=256) :: extra_binary_job_inlist5_name 

      namelist /binary_job/ &
         inlist_names, which_for_pgstar, &
         evolve_both_stars, stopping_age, trace_binary_rlo, &
         save_model_numbers, save_model_when_terminates, save_model_filenames, &
      ! extra files (Maybe overkill with so few inlist parameters)
         read_extra_binary_job_inlist1, extra_binary_job_inlist1_name, &
         read_extra_binary_job_inlist2, extra_binary_job_inlist2_name, &
         read_extra_binary_job_inlist3, extra_binary_job_inlist3_name, &
         read_extra_binary_job_inlist4, extra_binary_job_inlist4_name, &
         read_extra_binary_job_inlist5, extra_binary_job_inlist5_name

      contains

      subroutine do_run1_binary(tst, &
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

         use binary_mdot, only: eval_xfer_fraction
         use se_support, only: se_startup, se_finish_step, se_after_evolve
         use binary_evolve
         use mod_other_binary_jdot
         use binary_timestep
         use binary_history
         use binary_history_specs
         use run_star_support
         
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
            
!            subroutine data_for_extra_binary_history_columns(b, n, names, vals, ierr)
            subroutine data_for_extra_binary_history_columns(b, s, n, names, vals, ierr)
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
         
         integer :: id, id_extra, i, j, i_prev, result, result_reason, model_number, iounit
         type (star_info), pointer :: s
         character (len=256) :: restart_filename, photo_filename
         integer :: time0, time1, clock_rate
         logical :: doing_restart, first_try, continue_evolve_loop, just_did_backup
         real(dp) :: sum_times
         real(dp) :: dt
         real(dp) :: timestep_factor
         type (binary_info), pointer :: b
         
         include 'formats.inc'

         ierr = 0
         id_extra = 0
         call system_clock(time0,clock_rate)

         ! Find out if this is a restart
         iounit=alloc_iounit(ierr)
         open(unit=iounit, file='.restart', status='old', action='read',iostat=ierr)
         doing_restart = (ierr == 0)
         if (doing_restart) then
             read(iounit,'(a)', iostat=ierr) photo_filename ! same for both stars
             if (ierr /= 0) then
                 stop "Problem while reading restart info"
             end if
         else
             ierr = 0
         end if
         call free_iounit(iounit)

         call binary_ptr(b, ierr)

         call do_read_binary_job('inlist', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_binary_job'
            return
         end if

         if (.not. evolve_both_stars) num_stars=1
         
         call read_restart_controls('restart_inlist',ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_restart_controls'
            return
         end if
         
         write(*,*)
         write(*,*)

         result_reason = 0 ! BP: compiler warned may be used uninitialized

         ! Setup null hooks
         b% other_jdot_mb => null_other_jdot_mb
         b% other_jdot_gr => null_other_jdot_gr
         b% other_jdot_ml => null_other_jdot_ml
         b% other_jdot_tide => null_other_jdot_tide
         b% other_extra_jdot => null_other_extra_jdot
         b% other_jdot_ls => null_other_jdot_ls

         call do_one_binary_setup(b, 'inlist', ierr)
         ! extras_binary_controls is defined in run_binary_extras.f and hooks can
         ! be specified there
         call extras_binary_controls(ierr)
         
         do i = 1, num_stars
         
            call do_read_star_job(inlist_names(i), ierr)
            if (failed('do_read_star_job',ierr)) return

            id = id_from_read_star_job ! star allocated by do_read_star_job
            id_from_read_star_job = 0
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return

            save_model_numbers(i) = s% job% save_model_number
            save_model_when_terminates(i) = s% job% save_model_when_terminate
            save_model_filenames(i) = s% job% save_model_filename
            
            if (i == 1) then
               b% donor_id = id
               b% s1 => s
            else ! i > 1
               b% accretor_id = id
               b% s2 => s
            end if
         
            call starlib_init(s, ierr) ! okay to do extra calls on this
            if (failed('star_init',ierr)) return
         
            call star_set_kap_and_eos_handles(id, ierr)
            if (failed('set_star_kap_and_eos_handles',ierr)) return
            
            call star_setup(id, inlist_names(i), ierr)
            if (failed('star_setup',ierr)) return

            ! fix photo output for both stars to the one defined by binary
            s% photostep = b% photostep
            s% photo_digits = b% photo_digits

            restart_filename  = trim(s% photo_directory) // '/' // trim(photo_filename)
            
            b% star_ids(i) = id
            id_extra = b% star_extra_ids(i)
            call binary_extras_controls(s, extras_controls, ierr)
            if (failed('binary_extras_controls',ierr)) return
            
            call do_star_job_controls_before(id, s, doing_restart, ierr)
            if (failed('do_star_job_controls_before',ierr)) return    
                 
            call do_load1_star(id, s, doing_restart, restart_filename, ierr)
            if (failed('do_load1_star',ierr)) return         
            
            call do_star_job_controls_after(id, s, doing_restart, ierr)
            if (failed('do_star_job_controls_after',ierr)) return
            
            if (.not. doing_restart) then
               call before_evolve(id, ierr)
               if (failed('before_evolve',ierr)) return
            end if
            
            if (i == which_for_pgstar .or. which_for_pgstar < 0) then
               if (.not. doing_restart) then
                  call start_new_run_for_pgstar(s, ierr)
                  if (failed('start_new_run_for_pgstar',ierr)) return
               else
                  call show_terminal_header(id, ierr)
                  if (failed('show_terminal_header',ierr)) return
                  call restart_run_for_pgstar(s, ierr)
                  if (failed('restart_run_for_pgstar',ierr)) return
               end if
            end if
            
            b% star_extra_ids(i) = binary_extras_startup( &
               s, id, doing_restart, extras_startup, ierr)
            if (failed('binary_extras_startup',ierr)) return
            
            call se_startup(s, id, doing_restart, s% job% use_se_output, ierr)
            if (failed('se_startup',ierr)) return
         
            if (s% job% profile_starting_model) then
               write(*, '(a, i12)') 'save profile for model number ', s% model_number
               call save_profile(id, id_extra, &
                  how_many_extra_profile_columns, data_for_extra_profile_columns, &
                  3, ierr)
            end if
            
            s% doing_timing = .false.
            
            write(*,*)
            write(*,*)

         end do

         s% job% save_model_number = -111
         s% job% save_model_when_terminate = .false.
         s% job% save_model_filename = 'undefined'

         ! reread binary_job so model saving variables in it are preferred
         call do_read_binary_job('inlist', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_binary_job'
            return
         end if

         ! binary data must be initiated after stars
         if (.not. doing_restart) &
             call binarydata_init(b)
         b% evolve_both_stars = evolve_both_stars
         b% trace_binary_rlo = trace_binary_rlo
         call binary_history_column_names_init(ierr)
         call set_binary_history_columns(b, '', ierr)

         continue_evolve_loop = .true.
         s% doing_timing = .false.
         i_prev = 0

         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            do i = 1, num_stars

               id = b% star_ids(i)
               id_extra = b% star_extra_ids(i)
               call star_ptr(id, s, ierr)
               if (failed('star_ptr',ierr)) return
               
               if (s% model_number == s% job% first_model_for_timing) then
                  s% doing_timing = .true.
                  write(*,*) 'start timing'
                  write(*,*)
                  call system_clock(time0,clock_rate)
               end if
               
               if (s% job% auto_extend_net) then
                  call extend_net(s, ierr)
                  if (failed('extend_net',ierr)) return
               end if
               
               if (s% center_ye <= s% job% center_ye_limit_for_v_flag .and. .not. s% v_flag) then
                  write(*,1) 'have reached center ye limit', &
                     s% center_ye, s% job% center_ye_limit_for_v_flag
                  write(*,1) 'set v_flag true'
                  call star_set_v_flag(id, .true., ierr)
                  if (failed('star_set_v_flag',ierr)) return
                  if (ierr /= 0) return
               end if
             
               if (s% job% report_mass_not_fe56) call do_report_mass_not_fe56(s)
               if (s% job% report_cell_for_xm > 0) call do_report_cell_for_xm(s)
               
               model_number = get_model_number(id, ierr)
               if (failed('get_model_number',ierr)) return

            end do
             
            first_try = .true.
            just_did_backup = .false.

            call set_donor_star(b)
            step_loop: do ! may need to repeat this loop

               call set_star_timesteps()
               result = keep_going
               
               !get transfer fraction
               if (evolve_both_stars) then
                  b% companion_ratio = (b% s_accretor% photosphere_r*Rsun)/b% rl(2)
               else
                  b% companion_ratio = 0
               end if
               b% xfer_fraction = eval_xfer_fraction(s, b% mtransfer_rate, b% companion_ratio)

               do i = 1, num_stars

                  id = b% star_ids(i)
                  id_extra = b% star_extra_ids(i)
                  call star_ptr(id, s, ierr)
                  if (failed('star_ptr',ierr)) return

                  if (stop_now(s, i, id, id_extra)) then
                     result = terminate
                     result_reason = 0
                     exit step_loop
                  end if

                  result = worst_result(result, star_evolve_step(id, first_try, just_did_backup))

               end do

               call binary_evolve_step(b)

               do i = 1, num_stars
                  if (result == keep_going) then
                     id = b% star_ids(i)
                     id_extra = b% star_extra_ids(i)
                     call star_ptr(id, s, ierr)
                     result = worst_result(result, check_model(s, id, id_extra, extras_check_model))
                  end if
               end do
               if (result == keep_going) then
                  result = worst_result(result, binary_check_model(b))
               end if
               ! solve first binary timestep limit because star_pick_next_timestep needs it
               result = worst_result(result, binary_pick_next_timestep(b))
               if (result == keep_going) then
                  do i = 1, num_stars
                     id = b% star_ids(i)
                     call star_ptr(id, s, ierr)
                     if (failed('star_ptr',ierr)) return
                     result = worst_result(result, star_pick_next_timestep(id))
                  end do
               end if

               if (result == keep_going) then
                  call write_binary_history_info(b, &
                      how_many_extra_binary_history_columns, &
                      data_for_extra_binary_history_columns, ierr)
                  b% doing_first_model_of_run = .false.
                  exit step_loop
               end if

               do i = 1, num_stars

                  id = b% star_ids(i)
                  model_number = get_model_number(id, ierr)
                  if (failed('get_model_number',ierr)) return
                  
                  result_reason = get_result_reason(id, ierr)
                  if (result == retry) then
                     if (failed('get_result_reason',ierr)) return
                     if (s% job% report_retries) &
                        write(*,'(i6,3x,a,/)') model_number, &
                           'retry reason ' // trim(result_reason_str(result_reason))
                  else if (result == backup) then
                     if (failed('get_result_reason',ierr)) return
                     if (s% job% report_backups) &
                        write(*,'(i6,3x,a,/)') model_number, &
                           'backup reason ' // trim(result_reason_str(result_reason))
                  end if

               end do

               if (result == redo) then
                  do i = 1, num_stars
                     id = b% star_ids(i)
                     result = worst_result(result, star_prepare_to_redo(id))
                  end do
                  result = worst_result(result, binary_prepare_to_redo(b))
               end if
               if (result == retry) then
                  do i = 1, num_stars
                     id = b% star_ids(i)
                     result = worst_result(result, star_prepare_to_retry(id))
                  end do
                  result = worst_result(result, binary_prepare_to_retry(b))
               end if
               if (result == backup) then
                  do i = 1, num_stars
                     id = b% star_ids(i)
                     result = worst_result(result, star_do1_backup(id))
                  end do
                  result = worst_result(result, binary_do1_backup(b))
                  just_did_backup = .true.
               else
                  just_did_backup = .false.
               end if
               if (result == terminate) then
                  if (result_reason == result_reason_normal) then
                     do i = 1, num_stars
                        id = b% star_ids(i)
                        id_extra = b% star_extra_ids(i)
                        write(*, '(a, i12)') 'save profile for model number ', s% model_number
                        call save_profile(id, id_extra, &
                           how_many_extra_profile_columns, data_for_extra_profile_columns, &
                           3, ierr)
                     end do
                  end if
                  continue_evolve_loop = .false.
                  exit step_loop
               end if
               first_try = .false.
               
            end do step_loop

            if (result == keep_going) then
               ! finish binary step first or redos become inconsistent
               if(result == keep_going) result = binary_finish_step(b)
               do i = 1, num_stars
                  id = b% star_ids(i)
                  id_extra = b% star_extra_ids(i)
                  call star_ptr(id, s, ierr)
                  if (failed('star_ptr',ierr)) return
                  if (s% job% pgstar_flag .and. &
                        (i == which_for_pgstar) .or. (which_for_pgstar < 0)) &
                     call read_pgstar_inlist(s, inlist_names(i),ierr)
                  if (failed('read_pgstar_inlist',ierr)) return
                  result = extras_finish_step(s, id, id_extra)
                  if (result /= keep_going) exit evolve_loop
                  result = se_finish_step(s, id, s% job% use_se_output, &
                     how_many_extra_history_columns, data_for_extra_history_columns, &
                     how_many_extra_profile_columns, data_for_extra_profile_columns)
                  if (result /= keep_going) exit evolve_loop
                  result = star_finish_step(id, id_extra, .false., &
                        how_many_extra_profile_columns, data_for_extra_profile_columns, &
                        how_many_extra_history_columns, data_for_extra_history_columns, ierr)
                  if (failed('star_finish_step',ierr)) return
                  if (result /= keep_going) exit evolve_loop
                  if (s% job% pgstar_flag .and. &
                        (i == which_for_pgstar) .or. (which_for_pgstar < 0)) &
                     call update_pgstar_plots( &
                        s, .false., id_extra, &
                        how_many_extra_history_columns, &
                        data_for_extra_history_columns, &
                        how_many_extra_profile_columns, &
                        data_for_extra_profile_columns, &
                        ierr)
                  if (failed('update_pgstar_plots',ierr)) return
               end do
            else if (result == terminate) then
               if (result_reason == result_reason_normal) then
                  do i = 1, num_stars
                     id = b% star_ids(i)
                     id_extra = b% star_extra_ids(i)
                     call star_ptr(id, s, ierr)
                     if (failed('star_ptr',ierr)) return
                     result = star_finish_step( &
                        id, id_extra, s% job% save_photo_when_terminate, &
                        how_many_extra_profile_columns, data_for_extra_profile_columns, &
                        how_many_extra_history_columns, data_for_extra_history_columns, ierr)
                     if (failed('star_finish_step',ierr)) return
         
                     s% job% save_model_number = save_model_numbers(i)
                     s% job% save_model_filename = save_model_filenames(i)
                     if (save_model_when_terminates(i)) &
                        s% job% save_model_number = s% model_number
                     call do_saves( &
                        id, id_extra, s, &
                        how_many_extra_history_columns, &
                        data_for_extra_history_columns, &
                        how_many_extra_profile_columns, &
                        data_for_extra_profile_columns)
                  end do
                  call do_saves_for_binary_rlo
               end if
               exit evolve_loop
            end if
            
            do i = 1, num_stars
               id = b% star_ids(i)
               id_extra = b% star_extra_ids(i)
               call star_ptr(id, s, ierr)
         
               s% job% save_model_number = save_model_numbers(i)
               s% job% save_model_filename = save_model_filenames(i)
               call do_saves( &
                  id, id_extra, s, &
                  how_many_extra_history_columns, &
                  data_for_extra_history_columns, &
                  how_many_extra_profile_columns, &
                  data_for_extra_profile_columns)
            end do
            call do_saves_for_binary_rlo
            
         end do evolve_loop

         do i = 1, num_stars
         
            id = b% star_ids(i)
            id_extra = b% star_extra_ids(i)
            
            call star_ptr(id, s, ierr)
            if (failed('star_ptr',ierr)) return

            if (s% doing_timing) call show_times(id,s)
         
            result_reason = get_result_reason(id, ierr)
            if (result_reason /= result_reason_normal) then
               write(*, *) 
               write(*, *) 'terminated evolution because ' // trim(result_reason_str(result_reason))
               write(*, *)
            end if

            call extras_after_evolve(s, id, id_extra, ierr)
            if (failed('after_evolve_extras',ierr)) return

            call se_after_evolve(s, id, ierr)
            if (failed('se_after_evolve',ierr)) return

            if (s% job% pgstar_flag .and. &
                  (i == which_for_pgstar .or. which_for_pgstar < 0)) &
               call update_pgstar_plots( &
                  s, s% job% save_pgstar_files_when_terminate, id_extra, &
                  how_many_extra_history_columns, &
                  data_for_extra_history_columns, &
                  how_many_extra_profile_columns, &
                  data_for_extra_profile_columns, &
                  ierr)
            if (failed('update_pgstar_plots',ierr)) return

            call star_after_evolve(id, ierr)
            if (failed('star_after_evolve',ierr)) return

            call write_terminal_summary(id, ierr)
            if (failed('write_terminal_summary',ierr)) return
         
            call free_star(id, ierr)
            if (failed('free_star',ierr)) return
            
         end do
         
         call starlib_shutdown


         contains
                  
         
         subroutine do_saves_for_binary_rlo

            integer :: io, i
            character (len=strlen) :: str1, str2, filename, s1_photo_name, s2_photo_name

            i = index(b% s1% most_recent_photo_name,"/",.true.)
            if (i /= 0) then
               str1 = b% s1% most_recent_photo_name(i:)
            else
               str1 = b% s1% most_recent_photo_name
            end if
                  
            if (b% last_photo_filename /= str1) then
         
               b% last_photo_filename = str1
               if (evolve_both_stars) then
                  i = index(b% s2% most_recent_photo_name,"/",.true.)
                  if (i /= 0) then
                     str2 = b% s2% most_recent_photo_name(i:)
                  else
                     str2 = b% s2% most_recent_photo_name
                  end if
                  if (str1 /= str2) write(*,*) "WARNING: photos off sync"
               end if
               io = alloc_iounit(ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in alloc_iounit'
                  return
               end if
               filename = '.restart'
               open(unit=io, file=trim(filename), action='write', iostat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to open ' // trim(filename)
                  return
               end if
               write(io,'(a)') str1
               close(io)
               call free_iounit(io)
               
            end if
            
         end subroutine do_saves_for_binary_rlo


      end subroutine do_run1_binary   
      
      subroutine do_read_binary_job(filename, ierr)
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr  
         
         ierr = 0   
         
         include "binary_job.defaults"
         
         ierr = 0
         call read_inlist(filename, 1, ierr)
         
         if (ierr /= 0) then
            write(*,*) 'ierr from read_inlist ' // trim(filename)
            return
         end if
         
         !if (save_star_job_namelist) &
         !   call write_controls(star_job_namelist_name, ierr)

      end subroutine do_read_binary_job

      recursive subroutine read_inlist(filename, level, ierr)
         use utils_lib
         character(*), intent(in) :: filename
         integer, intent(in) :: level  
         integer, intent(out) :: ierr  
         
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=256) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra star_job inlist files'
            ierr = -1
            return
         end if
         
         ierr = 0
         unit=alloc_iounit(ierr)
         if (ierr /= 0) return
         
         open(unit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=binary_job, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(unit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=binary_job)
               close(unit)
            end if  
         end if
         call free_iounit(unit)
         if (ierr /= 0) return
         
         ! recursive calls to read other inlists
         
         read_extra1 = read_extra_binary_job_inlist1
         read_extra_binary_job_inlist1 = .false.
         extra1 = extra_binary_job_inlist1_name
         extra_binary_job_inlist1_name = 'undefined'
         
         read_extra2 = read_extra_binary_job_inlist2
         read_extra_binary_job_inlist2 = .false.
         extra2 = extra_binary_job_inlist2_name
         extra_binary_job_inlist2_name = 'undefined'
         
         read_extra3 = read_extra_binary_job_inlist3
         read_extra_binary_job_inlist3 = .false.
         extra3 = extra_binary_job_inlist3_name
         extra_binary_job_inlist3_name = 'undefined'
         
         read_extra4 = read_extra_binary_job_inlist4
         read_extra_binary_job_inlist4 = .false.
         extra4 = extra_binary_job_inlist4_name
         extra_binary_job_inlist4_name = 'undefined'
         
         read_extra5 = read_extra_binary_job_inlist5
         read_extra_binary_job_inlist5 = .false.
         extra5 = extra_binary_job_inlist5_name
         extra_binary_job_inlist5_name = 'undefined'
         
         if (read_extra1) then
            !write(*,*) 'read extra star_job inlist1 from ' // trim(extra1)
            call read_inlist(extra1, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra2) then
            !write(*,*) 'read extra star_job inlist2 from ' // trim(extra2)
            call read_inlist(extra2, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra3) then
            !write(*,*) 'read extra star_job inlist3 from ' // trim(extra3)
            call read_inlist(extra3, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra4) then
            !write(*,*) 'read extra star_job inlist4 from ' // trim(extra4)
            call read_inlist(extra4, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra5) then
            write(*,*) 'read extra star_job inlist5 from ' // trim(extra5)
            call read_inlist(extra5, level+1, ierr)
            if (ierr /= 0) return
         end if
         
      end subroutine read_inlist

      subroutine read_restart_controls(filename,ierr)
         use utils_lib
         character (len=*) :: filename
         integer, intent(out) :: ierr         
         integer :: unit         
         unit=alloc_iounit(ierr)
         if (ierr /= 0) return
         open(unit=unit, file=trim(filename), &
            action='read', status='old', delim='quote', iostat=ierr)
         if (ierr == 0) then
            read(unit, nml=binary_job, iostat=ierr)  
            close(unit)
         else
            ierr = 0 ! okay if no file
         end if
         call free_iounit(unit)
      end subroutine read_restart_controls

      integer function worst_result(result1, result2)
         integer, intent(in) :: result1, result2
         
         if(result1 == terminate .or. result2 == terminate) then
            worst_result = terminate
            return
         end if

         if(result1 == backup .or. result2 == backup) then
            worst_result = backup
            return
         end if
         
         if(result1 == retry .or. result2 == retry) then
            worst_result = retry
            return
         end if
         
         if(result1 == redo .or. result2 == redo) then
            worst_result = redo
            return
         end if

         worst_result = keep_going
         return
                              
      end function worst_result
      
      
      integer function check_model(s, id, id_extra, extras_check_model)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         interface
            integer function extras_check_model(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function extras_check_model
         end interface
         
         check_model = binary_extras_check_model(s, id, id_extra, extras_check_model)
         if (check_model /= keep_going) return

         check_model = star_check_model(id)
         if (check_model /= keep_going) then
            return
         end if
                              
      end function check_model
      
      
      logical function stop_now(s, i, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: i, id, id_extra
         
         stop_now = (s% star_age > stopping_age)
         if (stop_now) write(*,'(a)') 'all stars have reached stopping age'
      
      end function stop_now
      
      
! <<<<<<< "extras" interface routines
      
      
      subroutine binary_extras_controls(s, extras_controls, ierr)
         type (star_info), pointer :: s
         interface
            subroutine extras_controls(s, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(out) :: ierr
            end subroutine extras_controls      
         end interface
         integer, intent(out) :: ierr
         ierr = 0
         if (is_donor(s)) then
            call donor_controls(s, ierr)
         else
            call accretor_controls(s, ierr)
         end if
         if (ierr /= 0) return
         call extras_controls(s, ierr)
      end subroutine binary_extras_controls
      
      
      integer function binary_extras_startup(s, id, restart, extras_startup, ierr)
         use binary_tides, only : synch_spin_orbit_torque
         type (star_info), pointer :: s
         integer, intent(in) :: id
         logical, intent(in) :: restart
         interface
            integer function extras_startup(s, id, restart, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id
               logical, intent(in) :: restart
               integer, intent(out) :: ierr
            end function extras_startup
         end interface
         integer, intent(out) :: ierr

         s% use_other_adjust_mdot = .true.
         if (binary% do_j_accretion) then
            s% use_accreted_material_j = .true.
         end if
         s% accrete_given_mass_fractions = .true.
         s% accrete_same_as_surface = .false.
         if (binary% do_rotation .and. &
             (binary% do_tidal_synch .or. binary% do_initial_orbit_synch)) then
            s% other_torque => synch_spin_orbit_torque
         end if

         binary_extras_startup = extras_startup(s, id, restart, ierr)
      end function binary_extras_startup
      

      integer function binary_extras_check_model(s, id, id_extra, extras_check_model)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         interface
            integer function extras_check_model(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function extras_check_model
         end interface
         include 'formats.inc'
         binary_extras_check_model = extras_check_model(s, id, id_extra)
      end function binary_extras_check_model

      
! <<<<<<< accretor routines
      
      subroutine accretor_controls(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(b, ierr)
         s% initial_mass = b% m2
      end subroutine accretor_controls

! <<<<<<< donor routines
      
      subroutine donor_controls(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(b, ierr)
         s% initial_mass = b% m1
         s% other_photo_read => binary_photo_read
         s% other_photo_write => binary_photo_write
      end subroutine donor_controls

      subroutine binary_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         b => binary
         ierr = 0
         read(iounit, iostat=ierr) &
             b% mtransfer_rate, b% mtransfer_rate_old, b% mtransfer_rate_older, &
             b% angular_momentum_j, b% angular_momentum_j_old, b% angular_momentum_j_older, & 
             b% separation, b% separation_old, b% separation_older, &
             b% rl_relative_gap(1), b% rl_relative_gap_old(1), b% rl_relative_gap_older(1), &
             b% rl_relative_gap(2), b% rl_relative_gap_old(2), b% rl_relative_gap_older(2), &
             b% r(1), b% r_old(1), b% r_older(1), &
             b% r(2), b% r_old(2), b% r_older(2), &
             b% rl(1), b% rl_old(1), b% rl_older(1), &
             b% rl(2), b% rl_old(2), b% rl_older(2), &
             b% m(1), b% m_old(1), b% m_older(1), &
             b% m(2), b% m_old(2), b% m_older(2), &
             b% sum_div_qloc(1), b% sum_div_qloc_old(1), b% sum_div_qloc_older(1), &
             b% sum_div_qloc(2), b% sum_div_qloc_old(2), b% sum_div_qloc_older(2), &
             b% dt, b% dt_old, b% dt_older, &
             b% env, b% env_old, b% env_older, &
             b% xfer_fraction, b% xfer_fraction_old, b% xfer_fraction_older, &
             b% period, b% period_old, b% period_older, & 
             b% max_timestep, b% max_timestep_old, b% max_timestep_older, &
             b% change_factor, b% change_factor_old, b% change_factor_older, &
             b% min_binary_separation, &
             b% have_radiative_core(1), b% have_radiative_core_old(1), b% have_radiative_core_older(1), &
             b% have_radiative_core(2), b% have_radiative_core_old(2), b% have_radiative_core_older(2)
         if (ierr /= 0) stop "error in binary_photo_read"
      end subroutine binary_photo_read

      subroutine binary_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         type(binary_info), pointer :: b
         b => binary
         write(iounit) &
             b% mtransfer_rate, b% mtransfer_rate_old, b% mtransfer_rate_older, &
             b% angular_momentum_j, b% angular_momentum_j_old, b% angular_momentum_j_older, & 
             b% separation, b% separation_old, b% separation_older, &
             b% rl_relative_gap(1), b% rl_relative_gap_old(1), b% rl_relative_gap_older(1), &
             b% rl_relative_gap(2), b% rl_relative_gap_old(2), b% rl_relative_gap_older(2), &
             b% r(1), b% r_old(1), b% r_older(1), &
             b% r(2), b% r_old(2), b% r_older(2), &
             b% rl(1), b% rl_old(1), b% rl_older(1), &
             b% rl(2), b% rl_old(2), b% rl_older(2), &
             b% m(1), b% m_old(1), b% m_older(1), &
             b% m(2), b% m_old(2), b% m_older(2), &
             b% sum_div_qloc(1), b% sum_div_qloc_old(1), b% sum_div_qloc_older(1), &
             b% sum_div_qloc(2), b% sum_div_qloc_old(2), b% sum_div_qloc_older(2), &
             b% dt, b% dt_old, b% dt_older, &
             b% env, b% env_old, b% env_older, &
             b% xfer_fraction, b% xfer_fraction_old, b% xfer_fraction_older, &
             b% period, b% period_old, b% period_older, & 
             b% max_timestep, b% max_timestep_old, b% max_timestep_older, &
             b% change_factor, b% change_factor_old, b% change_factor_older, &
             b% min_binary_separation, &
             b% have_radiative_core(1), b% have_radiative_core_old(1), b% have_radiative_core_older(1), &
             b% have_radiative_core(2), b% have_radiative_core_old(2), b% have_radiative_core_older(2)

      end subroutine binary_photo_write





      end module run_binary_support
