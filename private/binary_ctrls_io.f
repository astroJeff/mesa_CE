! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton and Pablo Marchant
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module binary_ctrls_io
      
      use const_def
      use binary_def

      implicit none
      
      include "binary_controls.inc"      
      
      logical :: read_extra_binary_controls_inlist1
      character (len=256) :: extra_binary_controls_inlist1_name 
   
      logical :: read_extra_binary_controls_inlist2
      character (len=256) :: extra_binary_controls_inlist2_name 
   
      logical :: read_extra_binary_controls_inlist3
      character (len=256) :: extra_binary_controls_inlist3_name 
   
      logical :: read_extra_binary_controls_inlist4
      character (len=256) :: extra_binary_controls_inlist4_name 
   
      logical :: read_extra_binary_controls_inlist5
      character (len=256) :: extra_binary_controls_inlist5_name 
      
      namelist /binary_controls/ &
         m1, &
         m2, &
         initial_period_in_days, &
         initial_separation_in_Rsuns, &

! FOR COMMON ENVELOPE EVOLUTION
         alpha_CE, &
         do_CE, &

         fm, &
         fm_limit, &
         fa, &
         fr, &
         fr_limit, &
         fj, &
         dt_softening_factor, &
         
         cur_mdot_frac, &
         smallest_mass, &
         max_abs_mdot, &
         factor_for_contact_terminate, &
         
         max_yrs_dt, &
         
         transfer_fraction, &
         jdot_multiplier, &
         alpha, &
         tidal_Q, &
         R_companion, &

         upper_limit_on_period_in_hours, &
         lower_limit_on_period_in_hours, &

         col_depth_for_eps_extra, &
         irrad_flux_at_std_distance, &
         std_distance_for_irradiation, &

         include_accretor_mb, &
         limit_retention_by_mdot_edd, &
         companion_ratio_limit, &
         do_companion_ratio_limit, &
         accretion_powered_irradiation, &
         accretor_radius_for_irrad, &
         max_F_irr, &
         use_this_for_mdot_edd, &
         
         do_jdot_mb, &
         do_jdot_gr, &
         do_jdot_ml, &
         do_jdot_tide, &
         do_jdot_ls, &
         magnetic_braking_gamma, &
         
         mdot_scheme, &
         max_tries_to_achieve, &
         rl_rel_overlap_tolerance, &
         initial_change_factor, &
         change_factor_fraction, &
         implicit_lambda, &
         max_change_factor, &
         min_change_factor, &
         starting_mdot, &
         implicit_min_mdot, &
         do_rotation, &
         do_initial_orbit_synch, &
         do_tidal_synch, &
         do_j_accretion, &
         use_other_jdot_mb, &
         use_other_jdot_gr, &
         use_other_jdot_ml, &
         use_other_jdot_tide, &
         use_other_extra_jdot, &
         use_other_jdot_ls, &

         keep_donor_fixed, &
         
         history_name, &
         log_directory, &
         history_dbl_format, &
         history_int_format, &
         history_txt_format, &

         photostep, &
         photo_digits, &
         
         
      ! extra files
         read_extra_binary_controls_inlist1, extra_binary_controls_inlist1_name, &
         read_extra_binary_controls_inlist2, extra_binary_controls_inlist2_name, &
         read_extra_binary_controls_inlist3, extra_binary_controls_inlist3_name, &
         read_extra_binary_controls_inlist4, extra_binary_controls_inlist4_name, &
         read_extra_binary_controls_inlist5, extra_binary_controls_inlist5_name



      contains
      
      
      subroutine do_one_binary_setup(b, inlist, ierr)
         use utils_lib
         type (binary_info), pointer :: b
         character (len=*), intent(in) :: inlist
         integer, intent(out) :: ierr

         include 'formats'

         call set_default_binary_controls
         call read_binary_controls(b, inlist, ierr)

      end subroutine do_one_binary_setup


      subroutine read_binary_controls(b, filename, ierr)
         use utils_lib
         type (binary_info), pointer :: b
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr
         
         call read_binary_controls_file(b, filename, 1, ierr)
         
      end subroutine read_binary_controls
         
         
      recursive subroutine read_binary_controls_file(b, filename, level, ierr)
         use utils_lib
         character(*), intent(in) :: filename
         type (binary_info), pointer :: b
         integer, intent(in) :: level  
         integer, intent(out) :: ierr
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=256) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit 
         
         ierr = 0        
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra binary controls inlist files'
            ierr = -1
            return
         end if

         if (len_trim(filename) > 0) then
            unit=alloc_iounit(ierr); if (ierr /= 0) return
            open(unit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'Failed to open binary control namelist file ', trim(filename)
               return
            end if
            read(unit, nml=binary_controls, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 
               write(*, *) 
               write(*, *) 
               write(*, *) 
               write(*, '(a)') &
                  'Failed while trying to read binary control namelist file: ' // trim(filename)
               write(*, '(a)') &
                  'Perhaps the following runtime error message will help you find the problem.'
               write(*, *) 
               open(unit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=binary_controls)
               close(unit)
               call free_iounit(unit)
               return
            end if
            call free_iounit(unit)
         end if
         
         call store_binary_controls(b, ierr)
         
         ! recursive calls to read other inlists
         
         read_extra1 = read_extra_binary_controls_inlist1
         read_extra_binary_controls_inlist1 = .false.
         extra1 = extra_binary_controls_inlist1_name
         extra_binary_controls_inlist1_name = 'undefined'
         
         read_extra2 = read_extra_binary_controls_inlist2
         read_extra_binary_controls_inlist2 = .false.
         extra2 = extra_binary_controls_inlist2_name
         extra_binary_controls_inlist2_name = 'undefined'
         
         read_extra3 = read_extra_binary_controls_inlist3
         read_extra_binary_controls_inlist3 = .false.
         extra3 = extra_binary_controls_inlist3_name
         extra_binary_controls_inlist3_name = 'undefined'
         
         read_extra4 = read_extra_binary_controls_inlist4
         read_extra_binary_controls_inlist4 = .false.
         extra4 = extra_binary_controls_inlist4_name
         extra_binary_controls_inlist4_name = 'undefined'
         
         read_extra5 = read_extra_binary_controls_inlist5
         read_extra_binary_controls_inlist5 = .false.
         extra5 = extra_binary_controls_inlist5_name
         extra_binary_controls_inlist5_name = 'undefined'
         
         if (read_extra1) then
            write(*,*) 'read ' // trim(extra1)
            call read_binary_controls_file(b, extra1, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra2) then
            write(*,*) 'read ' // trim(extra2)
            call read_binary_controls_file(b, extra2, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra3) then
            write(*,*) 'read ' // trim(extra3)
            call read_binary_controls_file(b, extra3, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra4) then
            write(*,*) 'read ' // trim(extra4)
            call read_binary_controls_file(b, extra4, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra5) then
            write(*,*) 'read ' // trim(extra5)
            call read_binary_controls_file(b, extra5, level+1, ierr)
            if (ierr /= 0) return
         end if
         
      end subroutine read_binary_controls_file


      subroutine set_default_binary_controls
         include 'binary_controls.defaults'
      end subroutine set_default_binary_controls


      subroutine store_binary_controls(b, ierr)
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr
         
         ierr = 0
         
         b% m1 = m1
         b% m2 = m2
         b% initial_period_in_days = initial_period_in_days
         b% initial_separation_in_Rsuns = initial_separation_in_Rsuns

! FOR COMMON ENVELOPE EVOLUTION
         b% alpha_CE = alpha_CE
         b% do_CE = do_CE

         b% fm = fm
         b% fm_limit = fm_limit
         b% fa = fa
         b% fr = fr
         b% fr_limit = fr_limit
         b% fj = fj
         b% dt_softening_factor = dt_softening_factor
         
         b% cur_mdot_frac = cur_mdot_frac
         b% smallest_mass = smallest_mass
         b% max_abs_mdot = max_abs_mdot
         b% factor_for_contact_terminate = factor_for_contact_terminate
         
         b% max_yrs_dt = max_yrs_dt
         
         b% transfer_fraction = transfer_fraction
         b% jdot_multiplier = jdot_multiplier
         b% alpha = alpha
         b% tidal_Q = tidal_Q
         b% R_companion = R_companion

         b% upper_limit_on_period_in_hours = upper_limit_on_period_in_hours
         b% lower_limit_on_period_in_hours = lower_limit_on_period_in_hours

         b% col_depth_for_eps_extra = col_depth_for_eps_extra
         b% irrad_flux_at_std_distance = irrad_flux_at_std_distance
         b% std_distance_for_irradiation = std_distance_for_irradiation

         b% include_accretor_mb = include_accretor_mb
         b% limit_retention_by_mdot_edd = limit_retention_by_mdot_edd
         b% companion_ratio_limit = companion_ratio_limit
         b% do_companion_ratio_limit = do_companion_ratio_limit
         b% accretion_powered_irradiation = accretion_powered_irradiation
         b% accretor_radius_for_irrad = accretor_radius_for_irrad
         b% max_F_irr = max_F_irr
         b% use_this_for_mdot_edd = use_this_for_mdot_edd
         
         b% do_jdot_mb = do_jdot_mb
         b% do_jdot_gr = do_jdot_gr
         b% do_jdot_ml = do_jdot_ml
         b% do_jdot_tide = do_jdot_tide
         b% do_jdot_ls = do_jdot_ls
         b% magnetic_braking_gamma = magnetic_braking_gamma
         
         b% mdot_scheme = mdot_scheme
         b% max_tries_to_achieve = max_tries_to_achieve
         b% rl_rel_overlap_tolerance = rl_rel_overlap_tolerance
         b% initial_change_factor = initial_change_factor
         b% change_factor_fraction = change_factor_fraction
         b% implicit_lambda = implicit_lambda
         b% max_change_factor = max_change_factor
         b% min_change_factor = min_change_factor
         b% starting_mdot = starting_mdot
         b% implicit_min_mdot = implicit_min_mdot
         b% do_rotation = do_rotation
         b% do_initial_orbit_synch = do_initial_orbit_synch
         b% do_tidal_synch = do_tidal_synch
         b% do_j_accretion = do_j_accretion
         b% use_other_jdot_mb = use_other_jdot_mb
         b% use_other_jdot_gr = use_other_jdot_gr
         b% use_other_jdot_ml = use_other_jdot_ml
         b% use_other_jdot_tide = use_other_jdot_tide
         b% use_other_extra_jdot = use_other_extra_jdot
         b% use_other_jdot_ls = use_other_jdot_ls

         b% keep_donor_fixed = keep_donor_fixed

         b% history_name = history_name
         b% log_directory = log_directory

         b% history_dbl_format = history_dbl_format
         b% history_int_format = history_int_format
         b% history_txt_format = history_txt_format

         b% photostep = photostep
         b% photo_digits = photo_digits
         
      end subroutine store_binary_controls


      subroutine set_binary_controls_for_writing(b, ierr)
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr
         
         ierr = 0
         
         m1 = b% m1
         m2 = b% m2
         initial_period_in_days = b% initial_period_in_days
         initial_separation_in_Rsuns = b% initial_separation_in_Rsuns

! FOR COMMON ENVELOPE EVOLUTION
         alpha_CE = b% alpha_CE
         do_CE = b% do_CE

         fm = b% fm
         fm_limit = b% fm_limit
         fa = b% fa
         fr = b% fr
         fr_limit = b% fr_limit
         fj = b% fj
         dt_softening_factor = b% dt_softening_factor
         
         cur_mdot_frac = b% cur_mdot_frac
         smallest_mass = b% smallest_mass
         max_abs_mdot = b% max_abs_mdot
         factor_for_contact_terminate = b% factor_for_contact_terminate
         
         max_yrs_dt = b% max_yrs_dt
         
         transfer_fraction = b% transfer_fraction
         jdot_multiplier = b% jdot_multiplier
         alpha = b% alpha
         tidal_Q = b% tidal_Q
         R_companion = b% R_companion

         upper_limit_on_period_in_hours = b% upper_limit_on_period_in_hours
         lower_limit_on_period_in_hours = b% lower_limit_on_period_in_hours

         col_depth_for_eps_extra = b% col_depth_for_eps_extra
         irrad_flux_at_std_distance = b% irrad_flux_at_std_distance
         std_distance_for_irradiation = b% std_distance_for_irradiation

         include_accretor_mb = b% include_accretor_mb
         limit_retention_by_mdot_edd = b% limit_retention_by_mdot_edd
         companion_ratio_limit = b% companion_ratio_limit
         do_companion_ratio_limit = b% do_companion_ratio_limit
         accretion_powered_irradiation = b% accretion_powered_irradiation
         accretor_radius_for_irrad = b% accretor_radius_for_irrad
         max_F_irr = b% max_F_irr
         use_this_for_mdot_edd = b% use_this_for_mdot_edd
         
         do_jdot_mb = b% do_jdot_mb
         do_jdot_gr = b% do_jdot_gr
         do_jdot_ml = b% do_jdot_ml
         do_jdot_tide = b% do_jdot_tide
         do_jdot_ls = b% do_jdot_ls
         magnetic_braking_gamma = b% magnetic_braking_gamma
         
         mdot_scheme = b% mdot_scheme
         max_tries_to_achieve = b% max_tries_to_achieve
         rl_rel_overlap_tolerance = b% rl_rel_overlap_tolerance
         initial_change_factor = b% initial_change_factor
         change_factor_fraction = b% change_factor_fraction
         implicit_lambda = b% implicit_lambda
         max_change_factor = b% max_change_factor
         min_change_factor = b% min_change_factor
         starting_mdot = b% starting_mdot
         implicit_min_mdot = b% implicit_min_mdot
         do_rotation = b% do_rotation
         do_initial_orbit_synch = b% do_initial_orbit_synch
         do_tidal_synch = b% do_tidal_synch
         do_j_accretion = b% do_j_accretion
         use_other_jdot_mb = b% use_other_jdot_mb
         use_other_jdot_gr = b% use_other_jdot_gr
         use_other_jdot_ml = b% use_other_jdot_ml
         use_other_jdot_tide = b% use_other_jdot_tide
         use_other_extra_jdot = b% use_other_extra_jdot
         use_other_jdot_ls = b% use_other_jdot_ls

         keep_donor_fixed = b% keep_donor_fixed

         history_name = b% history_name
         log_directory = b% log_directory

         history_dbl_format = b% history_dbl_format
         history_int_format = b% history_int_format
         history_txt_format = b% history_txt_format

         photostep = b% photostep
         photo_digits = b% photo_digits
         
      end subroutine set_binary_controls_for_writing
      
      subroutine write_binary_controls(io,ierr)
         integer, intent(in) :: io
         integer, intent(out) :: ierr
         write(io, nml=binary_controls, iostat=ierr)  
      end subroutine write_binary_controls


      end module binary_ctrls_io

