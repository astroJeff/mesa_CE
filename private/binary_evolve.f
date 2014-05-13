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


      module binary_evolve

      use const_def
      use crlibm_lib
      use star_lib
      use star_def
      use crlibm_lib
      use binary_def

      implicit none

      contains

      subroutine binarydata_init(b)
         use utils_lib, only: is_bad_num
         type (binary_info), pointer :: b
         logical :: evolve_both_stars
         logical :: trace_binary_rlo
         integer :: finish_step_result
         include 'formats.inc'

         b% doing_first_model_of_run = .true.

         b% max_timestep = 1d99
         b% change_factor = b% initial_change_factor

         !TODO: use masses from stars to deal properly with reloads
         initial_mass(1) = b% m1
         initial_mass(2) = b% m2
         b% m(1) = initial_mass(1)*Msun
         b% m(2) = initial_mass(2)*Msun
         b% r(1) = Rsun*b% s1% photosphere_r
         if (b% evolve_both_stars) then
            b% r(2) = Rsun*b% s2% photosphere_r
         else
            b% r(2) = 0
         end if
         if (b% initial_period_in_days <= 0) then ! calculate from initial_separation_in_Rsuns
            b% separation = b% initial_separation_in_Rsuns*Rsun
            b% period = &
               (2*pi)*sqrt(b% separation*b% separation*b% separation/&
                     (standard_cgrav*(b% m(1)+b% m(2))))
            b% initial_period_in_days = b% period / (24d0*60d0*60d0)
         else
            b% period = b% initial_period_in_days*(24d0*60d0*60d0)
            b% separation = &
               pow_cr((b% s1% cgrav(1)*(b% m(1)+b% m(2)))*(b% period/(2*pi))**2,1d0/3d0)
         end if
         b% angular_momentum_j = &
            b% m(1) * b% m(2) * sqrt( b% s1% cgrav(1) * b% separation / (b% m(1) + b% m(2)) )

         b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl(2) = eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1))/b% rl(1)
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2))/b% rl(2)
         if (is_bad_num(b% rl_relative_gap(1))) stop 'binarydata_init'
         if (is_bad_num(b% rl_relative_gap(2))) stop 'binarydata_init'
         ! these will be adjusted properly by check_radiative_core
         b% have_radiative_core(1) = .true.
         b% have_radiative_core(2) = .true.

         write(*,*)
         write(*,1) 'm2', b% m2
         write(*,1) 'm1', b% m1
         write(*,1) 'initial_period_in_days', b% initial_period_in_days
         write(*,1) 'initial_separation_in_Rsun', b% separation/Rsun
         write(*,1) 'jdot_multiplier', b% jdot_multiplier
         write(*,1) 'fr', b% fr
         write(*,*)

         b% lower_limit_on_period_in_hours = -1d99
         b% have_radiative_core = .false.
         just_evolve = .false.
         b% s1% mesh_delta_coeff_pre_ms = 1
         min_binary_period = b% period
         b% min_binary_separation = b% separation

         b% num_tries = 0

         finish_step_result = binary_finish_step(b)
          
      end subroutine

      subroutine set_donor_star(b)
         use binary_mdot, only : donor_adjust_mdot, accretor_adjust_mdot
         use binary_tides, only : synch_spin_orbit_torque
         type (binary_info), pointer :: b
         include 'formats.inc'
          
         if((b% rl_relative_gap(1) > b% rl_relative_gap(2)) .or. b% keep_donor_fixed) then
            b% s_donor => b% s1
            b% s_accretor => b% s2
            b% d_i = 1
            b% a_i = 2
         else
            b% s_donor => b% s2
            b% s_accretor => b% s1
            b% d_i = 2
            b% a_i = 1
         end if
         b% s_donor% other_adjust_mdot => donor_adjust_mdot
         if (b% evolve_both_stars) b% s_accretor% other_adjust_mdot => accretor_adjust_mdot
      end subroutine

      subroutine binary_evolve_step(b)
         use utils_lib, only: is_bad_num
         use binary_jdot, only: get_jdot
         use binary_separation
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         
         integer :: ierr
         
         include 'formats.inc'

         !!check things before evolve
         !write(*,*) "check things before evolve: dt: ", b% s1% dt/secyer
         !write(*,*) "period:", b% period/(24*3600), b% period_old/(24*3600), b% period_older/(24*3600)
         !write(*,*) "sep:", b% separation/(rsun), b% separation_old/(rsun), b% separation_older/(rsun)
         !write(*,*) "am:", b% angular_momentum_j, b% angular_momentum_j_old, b% angular_momentum_j_older
         !write(*,*) "change_factor:", b% change_factor, b% change_factor_old, b% change_factor_older
         !write(*,*) "max timestep:", b% max_timestep, b% max_timestep_old, b% max_timestep_older
         !write(*,*) "m1:", b% m(1)/(msun)
         !write(*,*) "m2:", b% m(2)/(msun)
         !write(*,*) "sum of masses:", (b% m(2)+b% m(1))/(msun)
         !write(*,*) "r1:", b% r(1)/(rsun)
         !write(*,*) "r2:", b% r(2)/(rsun)
         !write(*,*) "rl1:", b% rl(1)/(rsun)
         !write(*,*) "rl2:", b% rl(2)/(rsun)

         b% m(1) = b% s1% mstar
         if (b% evolve_both_stars) then
            b% m(2) = b% s2% mstar
         else
            b% m(2) = b% m(2) &
               - b% xfer_fraction*b% mtransfer_rate*b% s1% dt
         end if
         
         b% r(1) = b% s1% photosphere_r*Rsun ! radius at photosphere in cm
         if (b% evolve_both_stars) b% r(2) = b% s2% photosphere_r*Rsun ! radius at photosphere in cm



! ------------------ New Functions to Calculate Ang Mom. Change ------------------- !

!         s => b% s_donor

         if (b% do_CE) b% check_CE = check_CE(b)
         
         if (b% do_CE .and. (dabs(b% mtransfer_rate) .gt. 1.0e-50)) then

            write(*,*) "MASS TRANSFER HAS STARTED"

            if (.not. b% started_rlof) then
!                call initial_CE_setup(b)
                b% started_rlof = .true.
                b% started_CE = .false.
                b% max_mdot_reached = 0.0
            end if

             if (b% mtransfer_rate .gt. b% max_mdot_reached) then
                b% max_mdot_reached = b% mtransfer_rate
             end if
             
         end if


         if (b% check_CE) then
            call new_separation_CE(b)
         else
            call new_separation_jdot(b)
         endif

         write(*,*) b% check_CE

         if (1.0 > 2.0) then
!         if (b% do_CE .and. b% started_rlof) then
            deallocate(b% CE_rho_old)
            deallocate(b% CE_P_old)
            deallocate(b% CE_vel_old)
            deallocate(b% CE_lnE_old)
            allocate(b% CE_rho_old(size(b% s1% rho)), stat=ierr)
            allocate(b% CE_P_old(size(b% s1% P)), stat=ierr)
            allocate(b% CE_vel_old(size(b% s1% v)), stat=ierr)
            allocate(b% CE_lnE_old(size(b% s1% lnE)), stat=ierr)
            b% CE_rho_old = b% s1% rho
            b% CE_P_old = b% s1% P
            b% CE_vel_old = b% s1% v
            b% CE_lnE_old = b% s1% lnE    
         end if
         
! ------------------ New Functions to Calculate Ang Mom. Change ------------------- !

        
        
         write(*,'(A,1pe16.9)') "New Separation = ", b% separation
 
         b% period = 2*pi*sqrt(b% separation*b% separation*b% separation/&
               (b% s1% cgrav(1)*(b% m(1)+b% m(2)))) 
         if (b% period < min_binary_period) min_binary_period = b% period
         
         ! use the new separation to calculate the new roche lobe radius
         
         b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
         b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
         b% rl_relative_gap(1) = (b% r(1) - b% rl(1))/b% rl(1) ! gap < 0 means out of contact
         b% rl_relative_gap(2) = (b% r(2) - b% rl(2))/b% rl(2) ! gap < 0 means out of contact

         if (is_bad_num(b% rl_relative_gap(1)) .or. is_bad_num(b% rl_relative_gap(2))) then
            stop 'error solving rl_rel_gap'
         end if
         
      end subroutine

      integer function binary_check_model(b)
         use binary_mdot, only: rlo_mdot, check_implicit_rlo
         use binary_irradiation
         use binary_separation
         type (binary_info), pointer :: b

         integer :: i, j, ierr
         logical :: implicit_rlo
         real(dp) :: new_mdot


         include 'formats.inc'

         binary_check_model = retry
         ierr = 0
         
         implicit_rlo = (b% max_tries_to_achieve > 0 .and. b% rl_rel_overlap_tolerance > 0d0)
         
         binary_check_model = keep_going
                  
         if (.not. just_evolve) then
            if (ierr /= 0) then
               binary_check_model = retry
               return
            end if
            if (implicit_rlo) then ! check agreement between new r and new rl
               b% s_donor% min_abs_mdot_for_change_limits = 1d99
               binary_check_model = check_implicit_rlo(new_mdot)
               if (binary_check_model == keep_going) then
                  b% num_tries = 0
               end if
            else
               new_mdot = rlo_mdot(b) ! grams per second
               ! smooth out the changes in mdot
               new_mdot = b% cur_mdot_frac*b% mtransfer_rate + (1-b% cur_mdot_frac)*new_mdot
               if (-new_mdot/(Msun/secyer) > b% max_abs_mdot) new_mdot = -b% max_abs_mdot*Msun/secyer 
            end if
            b% mtransfer_rate = new_mdot
            call adjust_irradiation(b% s_donor, b% mtransfer_rate, b% xfer_fraction)
         end if
         
         if (b% period/(60d0*60d0) < b% lower_limit_on_period_in_hours) then
            binary_check_model = terminate
            write(*,*) 'terminate because binary period < lower_limit_on_period_in_hours'
         end if
      
         if (b% period/(60d0*60d0) > b% upper_limit_on_period_in_hours) then
            binary_check_model = terminate
            write(*,*) 'terminate because binary period > upper_limit_on_period_in_hours'
         end if
         
         if (check_merger(b)) then
            deallocate(b% CE_rho_old)
            deallocate(b% CE_P_old)
            deallocate(b% CE_vel_old)
            deallocate(b% CE_lnE_old)  
            binary_check_model = terminate      
            write(*,*) 'terminate because system merged'
         end if

         if ((.not. b% started_CE) .and. b% r(b% d_i) .lt. 1.0d9) then
           write(*,*) "terminate because system avoided common envelope"
           binary_check_model = terminate               
         endif

         

         if (b% evolve_both_stars .and. b% s_accretor% photosphere_r*Rsun >= b% rl(b% a_i)) then
            if (b% s_accretor% photosphere_r*Rsun >= b% factor_for_contact_terminate * b% rl(b% a_i)) &
               binary_check_model = terminate
            write(*,'(a)') 'accretor photosphere has reached its Roche Lobe'
            write(*,2) 'accretor r/rl', b% s_accretor% model_number, b% s_accretor% photosphere_r*Rsun/b% rl(b% a_i)
            write(*,1) 'accretor photosphere_r', b% s_accretor% photosphere_r
            write(*,1) 'accretor rl/Rsun', b% rl(b% a_i)/Rsun
            write(*,1) 'donor photosphere_r', b% s_donor% photosphere_r
            write(*,1) 'sum photosphere_rs', b% s_accretor% photosphere_r + b% s_donor% photosphere_r
            write(*,1) 'sum rls', (b% rl(1) + b% rl(2))/Rsun
            write(*,1) 'binary_separation/Rsun', b% separation/Rsun
         end if
         !write(*,*) "Accretor J, accreted J, delta J: ", &
         !    s% total_angular_momentum, s% accreted_material_j * b% mtransfer_rate*s% dt, &
         !    s% total_angular_momentum - s% total_angular_momentum_old

      end function binary_check_model

      integer function binary_finish_step(b)
         type (binary_info), pointer :: b
         real(dp) :: spin_period

         if (b% do_tidal_synch .and. b% do_rotation) then
            if (mod(b% s1% model_number, b% s1% terminal_interval) == 0) then
               spin_period = 2*pi/b% s1% omega_avg_surf
               write(*,*) 'star_1_spin_period/binary_period', &
                  spin_period/b% period, spin_period/(60*60*24), b% period/(60*60*24)
               if (b% evolve_both_stars) then
                  spin_period = 2*pi/b% s2% omega_avg_surf
                  write(*,*) 'star_2_spin_period/binary_period', &
                     spin_period/b% period, spin_period/(60*60*24), b% period/(60*60*24)
               end if
            end if
         end if

         binary_finish_step = keep_going
         ! update change factor in case mtransfer_rate has changed
         if(b% mtransfer_rate_old /= b% mtransfer_rate .and. &
             b% mtransfer_rate /= 0 .and. b% mtransfer_rate_old/=0.0) then
            if(b% mtransfer_rate < b% mtransfer_rate_old) then
               b% change_factor = b% change_factor*(1.0-b% implicit_lambda) + b% implicit_lambda* &
                  (1+b% change_factor_fraction*(b% mtransfer_rate/b% mtransfer_rate_old-1))
            else
               b% change_factor = b% change_factor*(1.0-b% implicit_lambda) + b% implicit_lambda* &
                   (1+b% change_factor_fraction*(b% mtransfer_rate_old/b% mtransfer_rate-1))
            end if
            if(b% change_factor > b% max_change_factor) b% change_factor = b% max_change_factor
            if(b% change_factor < b% min_change_factor) b% change_factor = b% min_change_factor
         end if

         ! store all variables into "old" and "older"
         b% mtransfer_rate_older = b% mtransfer_rate_old
         b% angular_momentum_j_older = b% angular_momentum_j_old
         b% separation_older = b% separation_old
         b% dt_older = b% dt_old
         b% env_older = b% env_old
         b% xfer_fraction_older = b% xfer_fraction_old
         b% sum_div_qloc_older(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc_older(2) = b% sum_div_qloc_old(2)
         b% period_older = b% period_old
         b% rl_relative_gap_older(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap_older(2) = b% rl_relative_gap_old(2)
         b% r_older(1) = b% r_old(1)
         b% r_older(2) = b% r_old(2)
         b% rl_older(1) = b% rl_old(1)
         b% rl_older(2) = b% rl_old(2)
         b% m_older(1) = b% m_old(1)
         b% m_older(2) = b% m_old(2)
         b% have_radiative_core_older = b% have_radiative_core_old
         b% max_timestep_older = b% max_timestep_old
         b% change_factor_older = b% change_factor_old

         b% mtransfer_rate_old = b% mtransfer_rate
         b% angular_momentum_j_old = b% angular_momentum_j
         b% separation_old = b% separation
         b% dt_old = b% dt
         b% env_old = b% env
         b% xfer_fraction_old = b% xfer_fraction
         b% sum_div_qloc_old(1) = b% sum_div_qloc(1)
         b% sum_div_qloc_old(2) = b% sum_div_qloc(2)
         b% period_old = b% period
         b% rl_relative_gap_old(1) = b% rl_relative_gap(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap(2)
         b% r_old(1) = b% r(1)
         b% r_old(2) = b% r(2)
         b% rl_old(1) = b% rl(1)
         b% rl_old(2) = b% rl(2)
         b% m_old(1) = b% m(1)
         b% m_old(2) = b% m(2)
         b% have_radiative_core_old = b% have_radiative_core
         b% max_timestep_old = b% max_timestep
         b% change_factor_old = b% change_factor

      end function binary_finish_step

      integer function binary_prepare_to_redo(b)
         type (binary_info), pointer :: b

         binary_prepare_to_redo = redo
         ! restore variables
         ! do not restore mtransfer_rate during implicit rlo
         if (b% num_tries == 0) b% mtransfer_rate = b% mtransfer_rate_old
         b% angular_momentum_j = b% angular_momentum_j_old
         b% separation = b% separation_old
         b% dt = b% dt_old
         b% env = b% env_old
         b% xfer_fraction = b% xfer_fraction_old
         b% sum_div_qloc(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc(2) = b% sum_div_qloc_old(2)
         b% period = b% period_old
         b% rl_relative_gap(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap(2) = b% rl_relative_gap_old(2)
         b% r(1) = b% r_old(1)
         b% r(2) = b% r_old(2)
         b% rl(1) = b% rl_old(1)
         b% rl(2) = b% rl_old(2)
         b% m(1) = b% m_old(1)
         b% m(2) = b% m_old(2)
         b% have_radiative_core = b% have_radiative_core_old
         b% max_timestep = b% max_timestep_old
         b% change_factor = b% change_factor_old

      end function binary_prepare_to_redo

      integer function binary_prepare_to_retry(b)
         type (binary_info), pointer :: b

         binary_prepare_to_retry = retry
         ! restore variables
         b% mtransfer_rate = b% mtransfer_rate_old
         b% angular_momentum_j = b% angular_momentum_j_old
         b% separation = b% separation_old
         b% dt = b% dt_old
         b% env = b% env_old
         b% xfer_fraction = b% xfer_fraction_old
         b% sum_div_qloc(1) = b% sum_div_qloc_old(1)
         b% sum_div_qloc(2) = b% sum_div_qloc_old(2)
         b% period = b% period_old
         b% rl_relative_gap(1) = b% rl_relative_gap_old(1)
         b% rl_relative_gap(2) = b% rl_relative_gap_old(2)
         b% r(1) = b% r_old(1)
         b% r(2) = b% r_old(2)
         b% rl(1) = b% rl_old(1)
         b% rl(2) = b% rl_old(2)
         b% m(1) = b% m_old(1)
         b% m(2) = b% m_old(2)
         b% have_radiative_core = b% have_radiative_core_old
         b% max_timestep = b% max_timestep_old
         b% change_factor = b% change_factor_old

         b% num_tries = 0

      end function binary_prepare_to_retry

      integer function binary_do1_backup(b)
         type (binary_info), pointer :: b

         binary_do1_backup = retry
         ! restore variables
         b% mtransfer_rate = b% mtransfer_rate_older
         b% angular_momentum_j = b% angular_momentum_j_older
         b% separation = b% separation_older
         b% dt = b% dt_older
         b% env = b% env_older
         b% xfer_fraction = b% xfer_fraction_older
         b% sum_div_qloc(1) = b% sum_div_qloc_older(1)
         b% sum_div_qloc(2) = b% sum_div_qloc_older(2)
         b% period = b% period_older
         b% rl_relative_gap(1) = b% rl_relative_gap_older(1)
         b% rl_relative_gap(2) = b% rl_relative_gap_older(2)
         b% r(1) = b% r_older(1)
         b% r(2) = b% r_older(2)
         b% rl(1) = b% rl_older(1)
         b% rl(2) = b% rl_older(2)
         b% m(1) = b% m_older(1)
         b% m(2) = b% m_older(2)
         b% have_radiative_core = b% have_radiative_core_older
         b% max_timestep = b% max_timestep_older
         b% change_factor = b% change_factor_older

         b% mtransfer_rate_old = b% mtransfer_rate_older
         b% angular_momentum_j_old = b% angular_momentum_j_older
         b% separation_old = b% separation_older
         b% dt_old = b% dt_older
         b% env_old = b% env_older
         b% xfer_fraction_old = b% xfer_fraction_older
         b% sum_div_qloc_old(1) = b% sum_div_qloc_older(1)
         b% sum_div_qloc_old(2) = b% sum_div_qloc_older(2)
         b% period_old = b% period_older
         b% rl_relative_gap_old(1) = b% rl_relative_gap_older(1)
         b% rl_relative_gap_old(2) = b% rl_relative_gap_older(2)
         b% r_old(1) = b% r_older(1)
         b% r_old(2) = b% r_older(2)
         b% rl_old(1) = b% rl_older(1)
         b% rl_old(2) = b% rl_older(2)
         b% m_old(1) = b% m_older(1)
         b% m_old(2) = b% m_older(2)
         b% have_radiative_core_old = b% have_radiative_core_older
         b% max_timestep_old = b% max_timestep_older
         b% change_factor_old = b% change_factor_older

         b% num_tries = 0

      end function binary_do1_backup

      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow_cr(m1/m2,one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
      end function eval_rlobe
      
      subroutine initial_CE_setup(b)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer :: ierr
         
         s => b% s_donor

         s% dxdt_nuc_factor = 0d0  ! To stop any nuclear burning to help out the solvers
         s% mix_factor = 0d0  ! To stop mixing processes while in a CE
         b% xfer_fraction = 0d0  ! So that no amount of mass is accreted by the accretor

         allocate(b% CE_rho_old(size(s% rho)), stat=ierr)
         allocate(b% CE_P_old(size(s% P)), stat=ierr)
         allocate(b% CE_vel_old(size(s% v)), stat=ierr)
         allocate(b% CE_lnE_old(size(s% lnE)), stat=ierr)
         b% CE_rho_old = s% rho
         b% CE_P_old = s% P
         b% CE_vel_old = s% v
         b% CE_lnE_old = s% lnE    


      end subroutine initial_CE_setup

      end module binary_evolve
