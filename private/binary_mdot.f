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


      module binary_mdot

      use const_def
      use crlibm_lib
      use star_lib
      use star_def
      use binary_def

      implicit none

      contains

      integer function check_implicit_rlo(new_mdot)
         real(dp), intent(out) :: new_mdot
         
         type (star_info), pointer :: s
         real(dp), parameter :: min_abs_mdot = 0 ! may want non-zero someday.
         real(dp) :: rl_rel_overlap
         integer :: ierr
         
         type (binary_info), pointer :: b
         include 'formats.inc'
         ierr = 0
         call binary_ptr(b, ierr)
         s => b% s_donor
         
         ! NOTE: keep in mind that for mass loss, mdot is negative
         
         check_implicit_rlo = keep_going
         new_mdot = b% mtransfer_rate

         if (b% mdot_scheme == "roche_lobe") then
            rl_rel_overlap = b% rl_relative_gap(b% d_i)

            if (rl_rel_overlap < 0 .and. abs(b% mtransfer_rate) == 0) then
               if (implicit_rlo_dbg) &
                  write(*,1) 'implicit_rlo keep_going (no contact)', rl_rel_overlap
               return
            end if
         else if (b% mdot_scheme == "Ritter") then
            call get_info_for_ritter(b)
            rl_rel_overlap = (b% r(b% d_i)-b% rl(b% d_i))/b% rl(b% d_i) - &
                b% ritter_h*safe_log_cr(abs(b% mtransfer_rate / b% mdot_thin0))/b% rl(b% d_i)
            if (abs(b% mdot_thin) <= b% implicit_min_mdot*Msun/secyer) then
               new_mdot = b% mdot_thin
               if (implicit_rlo_dbg) &
                  write(*,1) 'implicit_rlo keep_going (mdot too low, using explicit scheme)', rl_rel_overlap
               return
            end if
         else if (b% mdot_scheme == "Kolb") then
            write(*,*) "mdot_scheme = Kolb not applicable for implicit scheme"
            write(*,*) "not transfering mass"
            new_mdot = 0
            return
         else
            write(*,*) "mdot_scheme = " , b% mdot_scheme , " not recognized"
            write(*,*) "not transfering mass"
            new_mdot = 0
            return
         end if
         
         if (0 <= rl_rel_overlap .and. &
             rl_rel_overlap <= b% rl_rel_overlap_tolerance) then
            if (implicit_rlo_dbg) &
               write(*,1) 'implicit_rlo keep_going (in tol)', rl_rel_overlap
            return
         end if
            
         if (b% num_tries == 0) then
            b% have_mdot_lo = .false.
            b% have_mdot_hi = .false.
            b% mdot_lo = 0
            b% mdot_hi = 0
            b% rel_overlap_mdot_hi = 0
            b% rel_overlap_mdot_lo = 0
         end if
         
         b% num_tries = b% num_tries + 1
         if (b% num_tries > b% max_tries_to_achieve) then
            check_implicit_rlo = retry
            if (implicit_rlo_dbg) &
               write(*,2) 'implicit_rlo retry (>max tries)', b% num_tries
            return
         end if
         
         ierr = 0
         !if (b% mdot_scheme == "roche_lobe" .or. abs(b% mdot_thin) >= b% implicit_min_mdot*Msun/secyer) then
            new_mdot = pick_mdot_for_implicit_rlo( &
               rl_rel_overlap, b% rl_rel_overlap_tolerance/2, b% mtransfer_rate, ierr)
         !end if
         if (ierr /= 0) then
            check_implicit_rlo = retry
            if (implicit_rlo_dbg) &
               write(*,*) 'check_implicit_rlo = retry (ierr)'
            return
         end if
         
         if (-new_mdot < b% implicit_min_mdot*Msun/secyer .and. rl_rel_overlap < 0 .and. &
             b% mdot_scheme == "roche_lobe") then
            if (implicit_rlo_dbg) &
               write(*,*) 'implicit_rlo = keep_going (no mdot)'
            new_mdot = -min_abs_mdot
            check_implicit_rlo = keep_going
            return
         end if
         if (-b% mdot_thin < b% implicit_min_mdot*Msun/secyer .and. &
             b% mdot_scheme == "Ritter") then
            if (implicit_rlo_dbg) &
               write(*,*) 'implicit_rlo = explicit scheme for Ritter (too low mdot)'
            check_implicit_rlo = keep_going
            return
         end if
         
         if (abs(new_mdot - b% mtransfer_rate) < 1d-6*abs(b% mtransfer_rate)) then
            if (implicit_rlo_dbg) &
               write(*,*) 'implicit_rlo = keep_going (tiny change in mdot)'
            new_mdot = b% mtransfer_rate
            check_implicit_rlo = keep_going
            return
         end if
         
         check_implicit_rlo = redo
         if (implicit_rlo_dbg) then
            write(*,2) 'check_implicit_rlo = redo', &
               b% num_tries, rl_rel_overlap, b% mtransfer_rate, b% mtransfer_rate_old
            write(*,*)
         end if

      end function check_implicit_rlo

      real(dp) function pick_mdot_for_implicit_rlo( &
            new_rel_overlap, rel_overlap_target, mdot_current, ierr) result(mdot_next)
         use num_lib, only: find0_quadratic, find0
         real(dp), intent(in) :: new_rel_overlap, rel_overlap_target, mdot_current
         integer, intent(out) :: ierr
         
         real(dp) :: starting_mdot, current_change_factor
         type(binary_info), pointer :: b
         include 'formats.inc'
         
         call binary_ptr(b,ierr)
         
         ! NOTE: keep in mind that for mass loss, mdot is negative
         
         
         starting_mdot = -b% starting_mdot*Msun/secyer
         current_change_factor = b% change_factor**(b% num_tries+1)

         if (implicit_rlo_dbg) then
            if (b% have_mdot_lo) then
               write(*,2) 'lo mdot, lo rel_overlap', &
                  b% num_tries, -b% mdot_lo/(Msun/secyer), b% rel_overlap_mdot_lo
            end if
            if (b% have_mdot_hi) then
               write(*,2) 'hi mdot, hi rel_overlap', &
                  b% num_tries, -b% mdot_hi/(Msun/secyer), b% rel_overlap_mdot_hi
            end if
         end if
         
         if (b% have_mdot_lo .and. b% have_mdot_hi) then ! quadratic interpolate
            ierr = 0
            mdot_next = find0_quadratic( &
               b% mdot_lo, b% rel_overlap_mdot_lo - rel_overlap_target, &
               mdot_current, new_rel_overlap - rel_overlap_target, &
               b% mdot_hi, b% rel_overlap_mdot_hi - rel_overlap_target, ierr)
            if (new_rel_overlap >= 0) then
               if (implicit_rlo_dbg) write(*,*) 'new mdot_lo'
               b% mdot_lo = mdot_current
               b% rel_overlap_mdot_lo = new_rel_overlap
            else
               if (implicit_rlo_dbg) write(*,*) 'new mdot_hi'
               b% mdot_hi = mdot_current
               b% rel_overlap_mdot_hi = new_rel_overlap
            end if
            if (ierr /= 0) then
               ierr = 0
               mdot_next = find0( &
                  b% mdot_lo, b% rel_overlap_mdot_lo - rel_overlap_target, &
                  b% mdot_hi, b% rel_overlap_mdot_hi - rel_overlap_target)
            end if
         else if (b% have_mdot_lo) then ! don't have mdot_hi
            if (new_rel_overlap < 0) then
               if (implicit_rlo_dbg) write(*,*) 'new mdot_hi'
               b% mdot_hi = mdot_current
               b% rel_overlap_mdot_hi = new_rel_overlap
               b% have_mdot_hi = .true.
               ! linear interpolate for mdot to match rel_overlap_target
               !mdot_next = find0( &
               !   mdot_lo, rel_overlap_mdot_lo - rel_overlap_target, &
               !   mdot_hi, rel_overlap_mdot_hi - rel_overlap_target)
               mdot_next = (b% mdot_hi+b% mdot_lo)/2.0
            else ! still too low
               if (implicit_rlo_dbg) write(*,*) 'new mdot_lo'
               b% mdot_lo = mdot_current
               mdot_next = b% mdot_lo*b% change_factor
               b% rel_overlap_mdot_lo = new_rel_overlap
            end if
         else if (b% have_mdot_hi) then ! don't have mdot_lo
            if (new_rel_overlap >= 0) then
               if (implicit_rlo_dbg) write(*,*) 'new mdot_lo'
               b% mdot_lo = mdot_current
               b% rel_overlap_mdot_lo = new_rel_overlap
               b% have_mdot_lo = .true.
               ! linear interpolate for mdot to match rel_overlap_target
               !mdot_next = find0( &
               !   mdot_lo, rel_overlap_mdot_lo - rel_overlap_target, &
               !   mdot_hi, rel_overlap_mdot_hi - rel_overlap_target)
               mdot_next = (b% mdot_hi+b% mdot_lo)/2.0
            else ! mdot still too high
               if (implicit_rlo_dbg) write(*,*) 'new mdot_hi'
               b% mdot_hi = mdot_current
               b% rel_overlap_mdot_hi = new_rel_overlap
               mdot_next = b% mdot_hi/b% change_factor
            end if
         else ! don't have either
            if (mdot_current > starting_mdot .and. (.not. abs(mdot_current) > 0)) then ! recall that both are negative
               mdot_next = starting_mdot
            else if (new_rel_overlap >= 0) then
               if (implicit_rlo_dbg) write(*,*) 'new mdot_lo'
               b% mdot_lo = mdot_current
               b% rel_overlap_mdot_lo = new_rel_overlap
               b% have_mdot_lo = .true.
               mdot_next = b% mdot_lo*b% change_factor
            else
               if (implicit_rlo_dbg) write(*,*) 'new mdot_hi'
               b% mdot_hi = mdot_current
               b% rel_overlap_mdot_hi = new_rel_overlap
               b% have_mdot_hi = .true.
               mdot_next = b% mdot_hi/b% change_factor
            end if
         end if
         
         if (.not. implicit_rlo_dbg) return
         
         write(*,2) 'next mdot, cur rel_overlap, prev mdot', &
            b% num_tries, -mdot_next/(Msun/secyer), new_rel_overlap, -mdot_current/(Msun/secyer), &
            b% r(b% d_i)/Rsun, b% rl(b% d_i)/Rsun
         !if (b% num_tries > 19) stop 'pick_mdot_for_implicit_rlo'
         
      end function pick_mdot_for_implicit_rlo


      real(dp) function eval_xfer_fraction(s, new_mdot, companion_ratio) result(xfer_fraction)
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_mdot, companion_ratio
         real(dp) :: mdot_edd, frac
         integer :: ierr
         type(binary_info), pointer :: b
         include 'formats.inc'
         
         call binary_ptr(b,ierr)
         xfer_fraction = b% transfer_fraction
         if (b% limit_retention_by_mdot_edd .and. new_mdot /= 0) then
            if (b% use_this_for_mdot_edd > 0) then
               mdot_edd = b% use_this_for_mdot_edd*(Msun/secyer)
            else
               mdot_edd = 4*pi*s% cgrav(1)*b% m(b% a_i)/(clight*s% opacity(1))
            end if
            xfer_fraction = min(b% transfer_fraction, mdot_edd/abs(new_mdot))
         end if
         if (b% evolve_both_stars .and. b% do_companion_ratio_limit)then
            ! xfer_fraction goes to 0 when companion fills its roche lobe
            if (companion_ratio >= 1) then
               xfer_fraction = 0
               if (b% trace_binary_rlo) write(*,1) '0 xfer_fraction', xfer_fraction
            else if (companion_ratio > b% companion_ratio_limit) then
               frac = (1 - companion_ratio)/(1 - b% companion_ratio_limit)
               xfer_fraction = min(xfer_fraction, 0.5*(1 - cospi_cr(frac)))
               if (b% trace_binary_rlo) write(*,1) 'limit xfer_fraction', xfer_fraction
            end if
         end if
         !write(*,2) 'xfer_fraction', s% model_number, xfer_fraction
      end function eval_xfer_fraction

      real(dp) function rlo_mdot(b) ! Adapted from a routine kindly provided by Anastasios Fragkos
         type(binary_info), pointer :: b
         include 'formats.inc'

         if (b% mdot_scheme == "roche_lobe") then
            write(*,*) "mdot_scheme = roche_lobe not applicable for explicit scheme"
            write(*,*) "not transfering mass"
            rlo_mdot = 0
            return
         else if (b% mdot_scheme /= "Ritter" .and. b% mdot_scheme /= "Kolb") then
            write(*,*) "mdot_scheme = " , b% mdot_scheme , " not recognized"
            write(*,*) "not transfering mass"
            rlo_mdot = 0
            return
         end if

         call get_info_for_ritter(b)
         rlo_mdot = b% mdot_thin

         if(b% mdot_scheme == "Kolb") then
             call get_info_for_kolb(b)
             rlo_mdot = rlo_mdot + b% mdot_thick
         end if

      end function rlo_mdot

      subroutine get_info_for_ritter(b)
         type(binary_info), pointer :: b
         real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(3.0 * kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% d_i)/b% m(b% a_i)
         F1 = (1.23d0  + 0.5D0* log10_cr(q))
         rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
         b% mdot_thin0 = (2.0D0*pi/exp_cr(1.0d0)) * v_th*v_th*v_th * &
             rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1   
         if (q < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q) - 0.038D0*(log10_cr(q))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q) + 0.114D0*(log10_cr(q))**2 )
         end if

         b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            if (b% ritter_exponent > 0) b% mdot_thin = -b% mdot_thin0
            if (b% ritter_exponent > -100.0D0) b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent) ! in gm per second
         else
            b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
         end if
         
      end subroutine get_info_for_ritter

      subroutine get_info_for_kolb(b)
         type(binary_info), pointer :: b
         real(dp) :: F3, FF, G1, x_L1, q, g
         real(dp) :: mdot_thick0,  R_gas, dP, rl, s_div_rl
         integer :: i, indexR
         include 'formats.inc'

         !--------------------- Optically thick MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL
         q = b% m(b% d_i)/b% m(b% a_i)
         b% mdot_thick = 0d0
         indexR=-1
         if(b% r(b% d_i)-b% rl(b% d_i) > 0.0d0) then
            i=1
            do while (b% s_donor% r(i) > b% rl(b% d_i))
               i=i+1
            end do
            indexR = i-1
            R_gas = 8.314472d7
            x_L1 = 0.5d0-0.227d0*log10_cr(q) !King, A. R., Frank, J., & Raines, D. J. 2002, Accretion Power in Astrophysics
            g = q/(x_L1*x_L1*x_L1)+1/((1d0-x_L1)*(1d0-x_L1)*(1d0-x_L1))
            rl = b% rl(b% d_i)
            s_div_rl = b% separation/rl
            FF = q*s_div_rl*s_div_rl*s_div_rl / sqrt(g*(g-1d0-q))
            mdot_thick0 = 2.0D0*pi*FF*rl*rl*rl/(b% s_donor% cgrav(1)*b% m(b% d_i))
            do i=1,indexR-1
               G1 = b% s_donor% gamma1(i)
               F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
               b% mdot_thick = b% mdot_thick + F3*sqrt(R_gas * b% s_donor% T(i) / &
                   b% s_donor% mu(i))*(b% s_donor% P(i+1)-b% s_donor% P(i))
            end do
            G1 = b% s_donor% gamma1(i)
            F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
            dP = (b% s_donor% r(indexR) - b% rl(b% d_i)) / &
                (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% P(i+1)-b% s_donor% P(i))
            b% mdot_thick = b% mdot_thick + F3*sqrt(R_gas * b% s_donor% T(i) / b% s_donor% mu(i))*dP
             
            b% mdot_thick = -mdot_thick0*b% mdot_thick
         end if

      end subroutine get_info_for_kolb

      subroutine donor_adjust_mdot(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         logical :: implicit_rlo
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(b, ierr)

         ! mdot_system_transfer is mass lost from the system due
         ! to inefficient mass transfer
         b% mdot_system_transfer(b% d_i) = 0

         if (b% mtransfer_rate >= 0) then
            b% mdot_system_transfer(b% a_i) = 0
            return
         else 
            b% mdot_system_transfer(b% a_i) = b% mtransfer_rate*(1-b% xfer_fraction)
         end if

         b% s_donor% mstar_dot = b% s_donor% mstar_dot + b% mtransfer_rate

         ! set accretion mode data to zero, will be adjusted if it corresponds
         ! in accretor_adjust_mdot

         b% accretion_mode = 0
         b% acc_am_div_kep_am = 0.0d0

      end subroutine donor_adjust_mdot

      subroutine accretor_adjust_mdot(id, ierr)
         use chem_def, only: chem_isos
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         logical, parameter :: dbg = .false.
         integer :: j
         real(dp) :: qratio, min_r
         type (binary_info), pointer :: b
         include 'formats.inc'

         ierr = 0
         call binary_ptr(b, ierr)

         b% s_accretor% mstar_dot = b% s_accretor% mstar_dot - b% mtransfer_rate*b% xfer_fraction

         !set angular momentum accretion as described in A.3.3 of de Mink et al. 2013
         if (b% do_j_accretion) then
            qratio = b% s_accretor% mstar / b% s_donor% mstar
            min_r = 0.0425*b% separation*pow_cr(qratio+qratio**2,0.25d0)

            !TODO: MUST USE EQUATORIAL RADIUS
            if (dbg) write(*,*) "radius, impact_radius, separation: ", &
                b% s_accretor% photosphere_r, min_r/rsun, b% separation/rsun
            if (b% s_accretor% photosphere_r*Rsun < min_r) then
               b% accretion_mode = 2
               b% s_accretor% accreted_material_j = &
                  sqrt(b% s_accretor% cgrav(1) * b% s_accretor% mstar * b% s_accretor% photosphere_r*Rsun) 
            else
               b% accretion_mode = 1
               b% s_accretor% accreted_material_j = &
                  sqrt(b% s_accretor% cgrav(1) * b% s_accretor% mstar * 1.7*min_r)
            end if
            b% acc_am_div_kep_am = b% s_accretor% accreted_material_j / &
                sqrt(b% s_accretor% cgrav(1) * b% s_accretor% mstar * b% s_accretor% photosphere_r*Rsun)
         end if

         !set accreted material composition
         b% s_accretor% num_accretion_species = b% s_donor% species
         do j = 1, b% s_donor% species
            b% s_accretor% accretion_species_id(j) = chem_isos% name(b% s_donor% chem_id(j))
            b% s_accretor% accretion_species_xa(j) = b% s_donor% xa(j,1)
         end do

         if (dbg) write(*,2) 'accretor mass_transfer', b% s_accretor% model_number, &
            b% mtransfer_rate, log10_cr(b% mtransfer_rate)
         
      end subroutine accretor_adjust_mdot

      end module binary_mdot
