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


      module binary_timestep

      use const_def
      use crlibm_lib
      use star_lib
      use star_def
      use binary_def

      implicit none

      contains

      subroutine set_star_timesteps() ! sets the smallest next timestep for all stars
         integer :: i
         real(dp) :: dt_min, rel_overlap
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         integer :: ierr
         ierr = 0
         call binary_ptr(b, ierr)
         dt_min = 1d99
         do i = 1, num_stars
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            if (s% dt_next < dt_min) then
               dt_min = s% dt_next
            end if
         end do
         dt_min = min(dt_min, b% max_timestep)
         do i = 1, num_stars
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            s% dt_next = dt_min
         end do

         if (.not. b% trace_binary_rlo) return
         call star_ptr(b% star_ids(1), s, ierr)
         if (ierr /= 0) then
             write(*, *) trim('star_ptr') // ' ierr', ierr
             return
         end if
         if (mod(s% model_number, s% terminal_interval) == 0) then
            write(*,'(99a20)') 'star', 'model', 'age', 'mass', 'lg_mdot', '(r-rl)/rl', 'last photo'
         else if (num_stars > 1) then
            call star_ptr(b% star_ids(2), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            if (mod(s% model_number, s% terminal_interval) == 0) &
               write(*,'(99a20)') 'star', 'model', 'age', 'mass', 'lg_mdot', '(r-rl)/rl', 'last photo'
         end if
         do i = 1, num_stars
            call star_ptr(b% star_ids(i), s, ierr)
            if (ierr /= 0) then
                write(*, *) trim('star_ptr') // ' ierr', ierr
                return
            end if
            rel_overlap = b% rl_relative_gap(i)
            write(*,'(2i20,4(4x,1pe16.9),8x,a)') i, &
               s% model_number, s% star_age, s% star_mass, &
               log10_cr(max(1d-99,abs(s% star_mdot))), rel_overlap, &
               trim(s% most_recent_photo_name)
         end do
         
      end subroutine set_star_timesteps

      integer function binary_pick_next_timestep(b)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         
         real(dp) :: &
            env_change, dtm, dtj, dta, dtr, &
            j_change, sep_change, rel_gap_change, set_dt
         character (len=8) :: why_str

         include 'formats.inc'

         dtm = 1d99
         dtj = 1d99
         dta = 1d99
         dtr = 1d99

         binary_pick_next_timestep = keep_going

         s => b% s_donor

         if (b% max_timestep < 0) b% max_timestep = b% s1% dt

         b% env = s% star_mass - s% he_core_mass 
         if (b% env_old /= 0) then
            env_change = b% env - b% env_old
         else
            env_change = 0
         end if
         
         if (b% rl_relative_gap_old(b% d_i) /= 0) then
            rel_gap_change = b% rl_relative_gap_old(b% d_i) - b% rl_relative_gap(b% d_i)
         else
            rel_gap_change = 0
         end if
         
         if (b% angular_momentum_j_old /= 0) then
            j_change = b% angular_momentum_j - b% angular_momentum_j_old
         else
            j_change = 0
         end if
         
         if (b% separation_old /= 0) then
            sep_change = b% separation - b% separation_old
         else
            sep_change = 0
         end if
   
         ! get limits for dt based on relative changes
         if (b% fm > 0) then
            dtm = s% time_step/(abs(env_change/max(b% env, b% fm_limit))/b% fm+1d-99)
         end if
         
         if (b% fr > 0) then
            dtr = s% time_step/ &
                (abs(rel_gap_change/max(-b% rl_relative_gap(b% d_i), b% fr_limit))/b% fr+1d-99)
         end if
         if (dtr < 10) dtr = 10

         if (b% fj > 0) then
            dtj = s% time_step/(abs(j_change/b% angular_momentum_j)/b% fj+1d-99)
         end if

         if (b% fa > 0) then
            dta = s% time_step/(abs(sep_change/b% separation)/b% fa+1d-99)
         end if

         set_dt = min(dtm, dtr, dtj, dta)
         
         if (set_dt == dtm) then
            why_str = 'dtm'
         else if (set_dt == dtr) then
            why_str = 'dtr'
         else if (set_dt == dtj) then
            why_str = 'dtj'
         else if (set_dt == dta) then
            why_str = 'dta'
         else
            why_str = '   '
         end if

         if (set_dt < 1d-7) set_dt = 1d-7 ! there's a limit to everything

         if (rlo_dbg) &
            write(*,'(i6,3x,a,3x,e20.10,12x,a,3x,f8.2,10x,a)') s% model_number, &
               '(rl-r)/r', b% rl_relative_gap(1), b% rl_relative_gap(2), &
               'signed_lg_rel_gap_x_1e4', &
               sign(1d0,b% rl_relative_gap)*log10_cr(max(1d0,1d6*abs(b% rl_relative_gap(1)))), &
               sign(1d0,b% rl_relative_gap)*log10_cr(max(1d0,1d6*abs(b% rl_relative_gap(2)))), &
               why_str

         b% max_timestep = exp10_cr(b% dt_softening_factor*log10_cr(set_dt*secyer) + &
             (1-b% dt_softening_factor)*log10_cr(b% max_timestep))
         
      end function binary_pick_next_timestep
      

      end module binary_timestep
