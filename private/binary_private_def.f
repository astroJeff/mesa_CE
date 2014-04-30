! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module binary_private_def
      
      use binary_def

      implicit none

   ! history column options

      integer, parameter :: bh_model_number = 1
      integer, parameter :: bh_age = bh_model_number + 1
      integer, parameter :: bh_donor_index = bh_age + 1
      integer, parameter :: bh_period_days = bh_donor_index + 1
      integer, parameter :: bh_period_hr = bh_period_days + 1
      integer, parameter :: bh_period_minutes = bh_period_hr + 1
      integer, parameter :: bh_lg_separation = bh_period_minutes + 1
      integer, parameter :: bh_binary_separation = bh_lg_separation + 1
      integer, parameter :: bh_rl_1 = bh_binary_separation + 1
      integer, parameter :: bh_rl_2 = bh_rl_1 + 1
      integer, parameter :: bh_rl_overflow_1 = bh_rl_2 + 1
      integer, parameter :: bh_rl_overflow_2 = bh_rl_overflow_1 + 1
      integer, parameter :: bh_rl_relative_overflow_1 = bh_rl_overflow_2 + 1
      integer, parameter :: bh_rl_relative_overflow_2 = bh_rl_relative_overflow_1 + 1
      integer, parameter :: bh_P_rot_div_P_orb_1 = bh_rl_relative_overflow_2 + 1
      integer, parameter :: bh_P_rot_div_P_orb_2 = bh_P_rot_div_P_orb_1 + 1
      integer, parameter :: bh_star_1_mass = bh_P_rot_div_P_orb_2 + 1
      integer, parameter :: bh_lg_star_1_mass = bh_star_1_mass + 1
      integer, parameter :: bh_star_2_mass = bh_lg_star_1_mass + 1
      integer, parameter :: bh_lg_star_2_mass = bh_star_2_mass + 1
      integer, parameter :: bh_sum_of_masses = bh_lg_star_2_mass + 1
      integer, parameter :: bh_lg_mstar_dot_1 = bh_sum_of_masses + 1
      integer, parameter :: bh_lg_mstar_dot_2 = bh_lg_mstar_dot_1 + 1
      integer, parameter :: bh_lg_system_mdot_1 = bh_lg_mstar_dot_2 + 1
      integer, parameter :: bh_lg_system_mdot_2 = bh_lg_system_mdot_1 + 1
      integer, parameter :: bh_lg_wind_mdot_1 = bh_lg_system_mdot_2 + 1
      integer, parameter :: bh_lg_wind_mdot_2 = bh_lg_wind_mdot_1 + 1
      integer, parameter :: bh_star_1_div_star_2_mass = bh_lg_wind_mdot_2 + 1
      integer, parameter :: bh_delta_star_1_mass = bh_star_1_div_star_2_mass + 1
      integer, parameter :: bh_delta_star_2_mass = bh_delta_star_1_mass + 1
      integer, parameter :: bh_lg_F_irr = bh_delta_star_2_mass + 1
      integer, parameter :: bh_xfer_fraction = bh_lg_F_irr + 1
      integer, parameter :: bh_log_J_orb = bh_xfer_fraction + 1
      integer, parameter :: bh_log_J_spin_1 = bh_log_J_orb + 1
      integer, parameter :: bh_log_J_spin_2 = bh_log_J_spin_1 + 1
      integer, parameter :: bh_log_J_total = bh_log_J_spin_2 + 1
      integer, parameter :: bh_log_abs_Jdot = bh_log_J_total + 1
      integer, parameter :: bh_log_abs_jdot_mb = bh_log_abs_Jdot + 1
      integer, parameter :: bh_log_abs_jdot_gr = bh_log_abs_jdot_mb + 1
      integer, parameter :: bh_log_abs_jdot_ml = bh_log_abs_jdot_gr + 1
      integer, parameter :: bh_log_abs_jdot_tide = bh_log_abs_jdot_ml + 1
      integer, parameter :: bh_log_abs_jdot_ls = bh_log_abs_jdot_tide + 1
      integer, parameter :: bh_log_abs_extra_jdot = bh_log_abs_jdot_ls + 1
      integer, parameter :: bh_accretion_mode = bh_log_abs_extra_jdot + 1
      integer, parameter :: bh_acc_am_div_kep_am = bh_accretion_mode + 1
      
      integer, parameter :: bh_col_id_max = bh_acc_am_div_kep_am
      
      character (len=maxlen_binary_history_column_name) :: binary_history_column_name(bh_col_id_max)
      
      contains
      
      
      subroutine binary_history_column_names_init(ierr)
         integer, intent(out) :: ierr
         
         integer :: i, cnt
         ierr = 0
         cnt = 0
         binary_history_column_name(:) = ''

         binary_history_column_name(bh_model_number) = 'model_number'
         binary_history_column_name(bh_age) = 'age'
         binary_history_column_name(bh_donor_index) = 'donor_index'
         binary_history_column_name(bh_period_days) = 'period_days'
         binary_history_column_name(bh_period_hr) = 'period_hr'
         binary_history_column_name(bh_period_minutes) = 'period_minutes'
         binary_history_column_name(bh_lg_separation) = 'lg_separation'
         binary_history_column_name(bh_binary_separation) = 'binary_separation'
         binary_history_column_name(bh_rl_1) = 'rl_1'
         binary_history_column_name(bh_rl_2) = 'rl_2'
         binary_history_column_name(bh_rl_overflow_1) = 'rl_overflow_1'
         binary_history_column_name(bh_rl_overflow_2) = 'rl_overflow_2'
         binary_history_column_name(bh_rl_relative_overflow_1) = 'rl_relative_overflow_1'
         binary_history_column_name(bh_rl_relative_overflow_2) = 'rl_relative_overflow_2'
         binary_history_column_name(bh_P_rot_div_P_orb_1) = 'P_rot_div_P_orb_1'
         binary_history_column_name(bh_P_rot_div_P_orb_2) = 'P_rot_div_P_orb_2'
         binary_history_column_name(bh_star_1_mass) = 'star_1_mass'
         binary_history_column_name(bh_lg_star_1_mass) = 'lg_star_1_mass'
         binary_history_column_name(bh_star_2_mass) = 'star_2_mass'
         binary_history_column_name(bh_lg_star_2_mass) = 'lg_star_2_mass'
         binary_history_column_name(bh_sum_of_masses) = 'sum_of_masses'
         binary_history_column_name(bh_lg_mstar_dot_1) = 'lg_mstar_dot_1'
         binary_history_column_name(bh_lg_mstar_dot_2) = 'lg_mstar_dot_2'
         binary_history_column_name(bh_lg_system_mdot_1) = 'lg_system_mdot_1'
         binary_history_column_name(bh_lg_system_mdot_2) = 'lg_system_mdot_2'
         binary_history_column_name(bh_lg_wind_mdot_1) = 'lg_wind_mdot_1'
         binary_history_column_name(bh_lg_wind_mdot_2) = 'lg_wind_mdot_2'
         binary_history_column_name(bh_star_1_div_star_2_mass) = 'star_1_div_star_2_mass'
         binary_history_column_name(bh_delta_star_1_mass) = 'delta_star_1_mass'
         binary_history_column_name(bh_delta_star_2_mass) = 'delta_star_2_mass'
         binary_history_column_name(bh_lg_F_irr) = 'lg_F_irr'
         binary_history_column_name(bh_xfer_fraction) = 'xfer_fraction'
         binary_history_column_name(bh_log_J_orb) = 'log_J_orb'
         binary_history_column_name(bh_log_J_spin_1) = 'log_J_spin_1'
         binary_history_column_name(bh_log_J_spin_2) = 'log_J_spin_2'
         binary_history_column_name(bh_log_J_total) = 'log_J_total'
         binary_history_column_name(bh_log_abs_Jdot) = 'log_abs_Jdot'
         binary_history_column_name(bh_log_abs_jdot_mb) = 'log_abs_jdot_mb'
         binary_history_column_name(bh_log_abs_jdot_gr) = 'log_abs_jdot_gr'
         binary_history_column_name(bh_log_abs_jdot_ml) = 'log_abs_jdot_ml'
         binary_history_column_name(bh_log_abs_jdot_tide) = 'log_abs_jdot_tide'
         binary_history_column_name(bh_log_abs_jdot_ls) = 'log_abs_jdot_ls'
         binary_history_column_name(bh_log_abs_extra_jdot) = 'log_abs_extra_jdot'
         binary_history_column_name(bh_accretion_mode) = 'accretion_mode'
         binary_history_column_name(bh_acc_am_div_kep_am) = 'acc_am_div_kep_am'
                  
         cnt = 0
         do i=1,bh_col_id_max
            if (len_trim(binary_history_column_name(i)) == 0) then
               write(*,*) 'missing name for log column id', i
               if (i > 1) write(*,*) 'following ' // trim(binary_history_column_name(i-1))
               write(*,*) 
               cnt = cnt+1
            end if
         end do

         if (cnt > 0) then
            ierr = -1
            return
         end if

      end subroutine binary_history_column_names_init         


      end module binary_private_def

