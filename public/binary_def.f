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
 

      module binary_def

      use star_lib
      use star_def
      use const_def
      
      implicit none

      logical, parameter :: rlo_dbg = .false.
      logical, parameter :: implicit_rlo_dbg = .false.
      
      integer, parameter :: rlo_info_alloc = 1
      integer, parameter :: rlo_info_get = 2
      integer, parameter :: rlo_info_put = 3

      real(dp) :: initial_binary_period ! (seconds)
      real(dp) :: min_binary_period ! (seconds)
      
      real(dp) :: initial_mass(2) ! (msun)

      logical :: just_evolve
      
      integer :: num_stars = 2

      integer, parameter :: maxlen_binary_history_column_name = 80
         
      !interfaces for procedure pointers
      abstract interface

         subroutine other_jdot_interface(ierr)
            integer, intent(out) :: ierr
         end subroutine other_jdot_interface

      end interface

      type binary_info
         include 'binary_data.inc'
         include 'binary_controls.inc'
      end type

      !type (binary_info), target :: binary_data
      type (binary_info), target, save :: binary
      !type (binary_info), pointer :: binary => binary_data
      
      
      contains

      subroutine binary_ptr(b, ierr)
         type (binary_info), pointer, intent(out) :: b
         integer, intent(out) :: ierr
         call get_binary_ptr(b, ierr)
      end subroutine binary_ptr


      subroutine get_binary_ptr(b,ierr)
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr         
         b => binary
         ierr = 0
      end subroutine get_binary_ptr
      

      logical function is_donor(s)
         type (star_info), pointer :: s
         is_donor = (s% id == binary% star_ids(1))
      end function is_donor
      
      !Maybe add this someday...
      !subroutine result_reason_init         
      !   result_reason_str(result_reason_normal) = 'normal'
      !   result_reason_str(dt_is_zero) = 'dt_is_zero'
      !   result_reason_str(nonzero_ierr) = 'nonzero_ierr'
      !   result_reason_str(hydro_failed_to_converge) = 'hydro_failed_to_converge'
      !   result_reason_str(do_burn_failed) = 'do_burn_failed'
      !   result_reason_str(diffusion_failed) = 'element_diffusion_failed'
      !   result_reason_str(too_many_steps_for_burn) = 'too_many_steps_for_burn'
      !   result_reason_str(too_many_steps_for_diffusion) = 'too_many_steps_for_diffusion'
      !   result_reason_str(too_many_steps_for_hydro) = 'too_many_steps_for_hydro'
      !   result_reason_str(adjust_mesh_failed) = 'adjust_mesh_failed'
      !   result_reason_str(adjust_mass_failed) = 'adjust_mass_failed'
      !   result_reason_str(core_dump_model_number) = 'core_dump_model_number'
      !   result_reason_str(timestep_limits) = 'convergence problems'
      !   result_reason_str(variable_change_limits) = 'variable_change_limits'
      !end subroutine result_reason_init
      
      subroutine init_binary_data
      
      end subroutine init_binary_data

      end module binary_def
