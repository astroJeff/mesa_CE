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


      module binary_jdot

      use const_def
      use star_lib
      use star_def
      use crlibm_lib
      use binary_def

      implicit none

      contains

      real(dp) function get_jdot(mdot, xfer_fraction)
         real(dp), intent(in) :: mdot, xfer_fraction

         integer :: ierr
         type (binary_info), pointer :: b

         call binary_ptr(b, ierr)

         ! calculate jdot from magnetic braking
         if (.not. b% do_jdot_mb) then
             b% jdot_mb = 0d0
         else if (.not. b% use_other_jdot_mb) then
             call default_jdot_mb(ierr)
         else
             call b% other_jdot_mb(ierr)
         end if
         
         ! calculate jdot from gravitational wave radiation
         if (.not. b% do_jdot_gr) then
             b% jdot_gr = 0d0
         else if (.not. b% use_other_jdot_gr) then
             call default_jdot_gr(ierr)
         else
             call b% other_jdot_gr(ierr)
         end if
            
         ! calculate jdot for mass ejected from system
         if (.not. b% do_jdot_ml) then
             b% jdot_ml = 0d0
         else if (.not. b% use_other_jdot_ml) then
             call default_jdot_ml(ierr)
         else
             call b% other_jdot_ml(ierr)
         end if
         
         ! calculate tidal jdot
         if (.not. b% do_jdot_tide) then
             b% jdot_tide = 0d0
         else if (.not. b% use_other_jdot_tide) then
             call default_jdot_tide(ierr)
         else
             call b% other_jdot_tide(ierr)
         end if
         
         ! calculate extra jdot
         if (.not. b% use_other_extra_jdot) then
             b% extra_jdot = 0
         else 
             call b% other_extra_jdot(ierr)
         end if

         ! solve jdot due to L-S coupling
         if (.not. b% do_jdot_ls) then
             b% jdot_ls = 0d0
         else if (.not. b% use_other_jdot_ls) then
             call default_jdot_ls(ierr)
         else
             call b% other_jdot_ls(ierr)
         end if
         
         get_jdot = (b% jdot_mb + b% jdot_gr + b% jdot_ml + &
             b% jdot_tide + b% extra_jdot)*b% jdot_multiplier + b% jdot_ls
         
      end function get_jdot

      subroutine check_radiative_core(b)
         use mlt_def, only: convective_mixing
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         
         real(dp) :: sum_conv, q_loc, sum_div_qloc 
         integer :: i, k, id

         include 'formats.inc'

         do i=1,2
            if (i == 1) then
               s => b% s_donor
               id = b% d_i
            else if (b% evolve_both_stars .and. b% include_accretor_mb) then
               s => b% s_accretor
               id = b% a_i
            else
               exit
            end if

            if (s% star_mass > 1.5) then
                b% have_radiative_core = .false.
                cycle
            end if

            ! calculate how much of inner region is convective
            sum_conv = 0; q_loc = 0
            do k = s% nz, 1, -1
               q_loc = s% q(k)
               if (q_loc > 0.5d0) exit 
               if (s% mixing_type(k) == convective_mixing) &
                  sum_conv = sum_conv + s% dq(k)
            end do
            
            sum_div_qloc = (b% sum_div_qloc(id) + sum_conv/q_loc)/2
            b% sum_div_qloc(id) = sum_div_qloc
            
            if (b% have_radiative_core(id)) then ! check if still have rad core
               if (sum_div_qloc > 0.75d0) then
                  b% have_radiative_core(id) = .false.
                  write(*,*)
                  write(*,*) 'turn off magnetic braking because radiative core has gone away'
                  write(*,*)
                  ! required mdot for the implicit scheme may drop drastically,
                  ! so its neccesary to increase change factor to avoid implicit 
                  ! scheme from getting stuck
                  b% change_factor = b% max_change_factor
               end if
            else if (sum_div_qloc < 0.25d0) then ! check if now have rad core
               if (.not. b% have_radiative_core(id)) then
                  write(*,*)
                  write(*,*) 'turn on magnetic braking'
                  write(*,*)
               end if
               b% have_radiative_core(id) = .true.
            end if
         end do
            
     end subroutine check_radiative_core

      subroutine default_jdot_mb(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: rsun4,two_pi_div_p3
         ierr = 0
         call binary_ptr(b, ierr)
         b% jdot_mb = 0
         rsun4 = rsun*rsun*rsun*rsun
         call check_radiative_core(b)
         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)
         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
         !if (b% have_radiative_core(b% d_i)) b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
         b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
                           pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
                           two_pi_div_p3

         if (b% evolve_both_stars .and. b% include_accretor_mb .and. &
             b% have_radiative_core(b% a_i)) then
             b% jdot_mb = b% jdot_mb - &
                           3.8d-30*b% m(b% a_i)*rsun4* &
                           pow_cr(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
                           two_pi_div_p3
         end if

      end subroutine default_jdot_mb

      subroutine default_jdot_gr(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: bs4, clight5, cgrav3
         ierr = 0
         call binary_ptr(b, ierr)
         bs4 = b% separation*b% separation*b% separation*b% separation
         clight5 = clight*clight*clight*clight*clight
         cgrav3 = b% s_donor% cgrav(1)*b% s_donor% cgrav(1)*b% s_donor% cgrav(1)
         b% jdot_gr = -32d0 * cgrav3 * b% m(b% a_i) * b% m(b% d_i) * (b% m(b% a_i) + b% m(b% d_i)) / &
             (5d0 * clight5 * bs4) * b% angular_momentum_j
      end subroutine default_jdot_gr

      subroutine default_jdot_ml(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: alfa
         ierr = 0
         call binary_ptr(b, ierr)
         if (b% alpha < 0) then
            alfa = (b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i)))**2
         else
            alfa = b% alpha
         end if
         b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             alfa*(b% separation**2)*2*pi/b% period
         if (b% alpha < 0) then
            alfa = (b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i)))**2
         else
            alfa = b% alpha
         end if
         b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*&
             alfa*(b% separation**2)*2*pi/b% period
      end subroutine default_jdot_ml

      subroutine default_jdot_tide(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: rc, rcbs5
         ierr = 0
         call binary_ptr(b, ierr)
         b% jdot_tide = 0
         if (b% tidal_Q <= 0 .or. b% period <= 0) then
            b% jdot_tide = 0
         else
            rc = b% R_companion
            if (rc <= 0) rc = pow_cr(b% m(b% a_i),0.78d0)
            rcbs5 = (rc*Rsun/b% separation)*(rc*Rsun/b% separation)* &
                  (rc*Rsun/b% separation)*(rc*Rsun/b% separation)*(rc*Rsun/b% separation)
            b% jdot_tide = -(9d0/2d0)*(2*pi/b% period)* &
               (b% m(b% d_i)/b% m(b% a_i))*(1/b% tidal_Q)*rcbs5 * b% angular_momentum_j
         end if
      end subroutine default_jdot_tide

      subroutine default_jdot_ls(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: delta_J
         ierr = 0
         call binary_ptr(b, ierr)
         b% jdot_ls = 0
         ! ignore in first step, or if not doing rotation
         if (b% doing_first_model_of_run .or. .not. b% do_rotation) &
             return
         ! bulk change in spin angular momentum takes tides into account
         delta_J = b% s_donor% total_angular_momentum_old - &
             b% s_donor% total_angular_momentum
         ! ignore angular momentum lost through winds
         if (b% s_donor% mstar_dot < 0) &
             delta_J = delta_J + b% s_donor% angular_momentum_removed * &
                 abs(b% mdot_system_wind(b% d_i) / b% s_donor% mstar_dot)
         ! Repeat for accretor
         if (b% evolve_both_stars) then
             delta_J = delta_J + b% s_accretor% total_angular_momentum_old - &
                 b% s_accretor% total_angular_momentum
             if (b% s_accretor% mstar_dot < 0) then
                 ! all AM lost from the accretor is lost from the system
                 delta_J = delta_J + b% s_accretor% angular_momentum_removed
             end if
         end if
         b% jdot_ls = delta_J / b% s_donor% dt
      end subroutine default_jdot_ls

      end module binary_jdot
