!   This is a test program to attempt to externalize the routines in /mesa/binary
!   for the calculation of orbital separation evolution. This is not necessary, but 
!   done to make the code more modular, specifically to ease the inclusion of
!   common envelope evolution into MESA. Woo Woo!


      module binary_separation
            
      use star_lib
      use crlibm_lib
      use star_def
      use const_def
      use binary_def
      
      implicit none

      contains
      
      
      logical function check_CE(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: mdot, rl2, R_accretor
         include 'formats.inc'
         
         b => binary
         s => b% s_donor

         mdot = b% mtransfer_rate_old

         R_accretor = accretor_mass_radius_relation(b% m2)
         rl2 = b% rl(b% a_i)

         if( dabs(mdot) > 1.0d0*(s% mstar)*dsqrt(standard_cgrav * s% rho(1))) then
            ! Dynamically Unstable RLOF
            check_CE = .true.
         else if (b% r(b% d_i) + R_accretor > b% separation) then
            ! Contact Binary
            check_CE = .true.
         else if (b% rl(b% a_i) < b% r(b% a_i)) then
            ! Double Common Envelope
            check_CE = .true.
         else
            check_CE = .false.
         endif
     
      end function check_CE

      logical function check_merger(s)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: R_donor, R_accretor
         include 'formats.inc'
         
         b => binary
         s => b% s_donor

         R_accretor = accretor_mass_radius_relation(b% m(b% a_i))

         if (s% center_he4 .lt. 1d-3) then      ! AGB star
            R_donor = s% c_core_radius*Rsun
         else if (s% center_h1 .lt. 1d-3) then  ! RGB star
            R_donor = s% he_core_radius*Rsun
         else                                   ! MS star
            R_donor = s% r(1)
         endif

         if (R_accretor + R_donor .gt. b% separation) then
            check_merger = .true.
         else
            check_merger = .false.
         endif
      end function check_merger

      real(dp) function accretor_mass_radius_relation(mass)
         real(dp), intent(in) :: mass

         ! Basic mass-radius relation for stars less than 2 Msun
         accretor_mass_radius_relation = Rsun * (mass/Msun) ** (0.9)

      end function accretor_mass_radius_relation

      real(dp) function calc_naive_A_f(s)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: E_bind_ml, E_orb_init, af_E_orb_f
         real(dp) :: E_int_spec, E_bind_spec, E_kin_spec, E_enth_spec
         real(dp) :: Delm2
         integer :: i
         include 'formats.inc'

         b => binary
         s => b% s_donor

         ! ------------ A naive calculation of A_f --------------- !

         E_bind_ml = 0.0d0
         i = 1
         Delm2 = s% mstar - s% he_core_mass*Msun
         do while ((s% mstar * (1.0d0 - s% q(i)) .le. Delm2) .and. (i < s% nz))
            E_int_spec = dexp(s% lnE(i))
            E_bind_spec = - standard_cgrav * s% mstar / s% r(i)
            E_kin_spec = s% v(i) * s% v(i) / 2.0
            E_enth_spec = s% P(i)/s% rho(i)

            E_bind_ml = E_bind_ml - (E_int_spec + E_bind_spec + E_kin_spec + E_enth_spec) &
                        * s% mstar * s% dq(i)

            i = i +1
         enddo

         if (i > 1) then
            E_int_spec =  dexp(s% lnE(i))
            E_bind_spec = - standard_cgrav * s% mstar / s% r(i)
            E_kin_spec = s% v(i) * s% v(i) / 2.0
            E_enth_spec = s% P(i)/s% rho(i)

            E_bind_ml = E_bind_ml - (E_int_spec + E_bind_spec + E_kin_spec + E_enth_spec) &
                         * (Delm2 - s% mstar * (1.0d0 - s% q(i-1)))
         endif

         E_orb_init = b% alpha_CE * standard_cgrav * b% m(b% a_i)*Msun * s% mstar / (2.0d0 * b% separation)
         af_E_orb_f = b% alpha_CE * standard_cgrav * b% m(b% a_i)*Msun * s% he_core_mass*Msun / 2.0d0
         calc_naive_A_f = af_E_orb_f / (E_orb_init + E_bind_ml)
 

      end function calc_naive_A_f

      subroutine new_separation_CE(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         include 'formats.inc'
         
         
      end subroutine new_separation_CE


      subroutine new_separation_jdot(b)
         use binary_jdot, only: get_jdot
         type(binary_info), pointer :: b
         include 'formats.inc'



         ! solve the winds in the system for jdot calculation,
         ! these don't include mass lost due to mass_transfer_efficiency < 1.0
         b% mdot_system_wind(b% d_i) = b% s_donor% mstar_dot - b% mtransfer_rate
         if (b% evolve_both_stars) then
            b% mdot_system_wind(b% a_i) = b% s_accretor% mstar_dot + &
                b% mtransfer_rate * b% xfer_fraction
         else
            b% mdot_system_wind(b% a_i) = 0.0d0
         end if

         ! get jdot and update orbital J
         b% jdot = get_jdot(b% mtransfer_rate, b% xfer_fraction)
         b% angular_momentum_j = b% angular_momentum_j + b% jdot*b% s1% time_step*secyer

         if (b% angular_momentum_j <= 0) then
            stop 'bad angular_momentum_j'
         end if
         
         ! use the new j to calculate new separation
         
         b% separation = ((b% angular_momentum_j/(b% m(1)*b% m(2)))**2) *&
             (b% m(1)+b% m(2)) / b% s1% cgrav(1)
         if (b% separation < b% min_binary_separation) &
            b% min_binary_separation = b% separation





      end subroutine new_separation_jdot


      end module binary_separation