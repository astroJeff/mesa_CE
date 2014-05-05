          ! state saved and restored at each step (only for donor if doing binary)
          
          real(dp) :: &
             mtransfer_rate, mtransfer_rate_old, mtransfer_rate_older, &
             angular_momentum_j, angular_momentum_j_old, angular_momentum_j_older, & 
             ! separation in cms
             separation, separation_old, separation_older, &
             ! (r-rl)/rl for each star
             rl_relative_gap(2), rl_relative_gap_old(2), rl_relative_gap_older(2), &
             ! photosphere radius for each star (as computed by the evolution code) (cm)
             r(2), r_old(2), r_older(2), &
             ! roche lobe radius for each star (as derived from masses and separation) (cm)
             rl(2), rl_old(2), rl_older(2), &
             ! star masses
             m(2), m_old(2), m_older(2), &
             sum_div_qloc(2), sum_div_qloc_old(2), sum_div_qloc_older(2), &
             dt, dt_old, dt_older, &
             env, env_old, env_older, &
             xfer_fraction, xfer_fraction_old, xfer_fraction_older, &
             period, period_old, period_older, & 
             max_timestep, max_timestep_old, max_timestep_older, &
             change_factor, change_factor_old, change_factor_older, &
             min_binary_separation

          logical :: &
             have_radiative_core(2), have_radiative_core_old(2), have_radiative_core_older(2)
        
          ! INFO FOR COMMON ENVELOPE EVOLUTION
!          real(dp) :: E_kin_spec, E_enth_spec, E_grav_spec
!          real(dp) :: Delm2, E_bind_ml, E_orb_init, af_E_orb_fin
!          real(dp) :: alpha_CE

          real(dp), pointer, dimension(:) :: CE_lnE_old, CE_rho_old, CE_P_old, CE_vel_old
          real(dp) :: max_mdot_reached

          logical :: in_CE_flag, merger_flag, started_CE_flag, prev_CE_flag
          logical :: check_CE 
          logical :: started_rlof, started_CE

          ! *************** end of info that is saved and restored for restarts  ***************

          ! work variables -- not saved or restored. must be set at each step.
             
          real(dp) :: jdot, jdot_mb, jdot_gr, jdot_ml, jdot_tide, jdot_ls, extra_jdot

          real(dp) :: companion_ratio

          !rate of mass loss from the system due to inefficient mass accretion by each star
          real(dp) :: mdot_system_transfer(2)
          !rate of mass loss from the system due to winds
          real(dp) :: mdot_system_wind(2)
          
         
      integer :: star_ids(2), star_extra_ids(2)
      ! pointers that are adjusted to donor and accretor
      type (star_info), pointer :: s_donor, s_accretor
      ! pointers to fixed stars
      type (star_info), pointer :: s1, s2
      
      integer :: donor_id, accretor_id
      integer :: d_i, a_i ! array indexes for the donor and the accretor
      character (len=strlen) :: last_photo_filename

      logical :: evolve_both_stars
      logical :: trace_binary_rlo
            
      ! info for implicit rlo
      real(dp) :: mdot_lo, mdot_hi, rel_overlap_mdot_hi, rel_overlap_mdot_lo
      logical :: have_mdot_lo = .false., have_mdot_hi = .false.
      integer :: num_tries = 0

      ! info for Ritter and Kolb & Ritter prescriptions
      real(dp) :: ritter_h, ritter_exponent, mdot_thin, mdot_thin0, mdot_thick

      ! type of accretion (none=0, ballistic=1, Keplerian=2)
      integer :: accretion_mode = 0
      ! ratio of accreted specific angular momentum to that of a Keplerian orbit at R star
      real(dp) :: acc_am_div_kep_am

      ! "other" procedures
      procedure(other_jdot_interface), pointer, nopass :: &
          other_jdot_mb => null()

      procedure(other_jdot_interface), pointer, nopass :: &
          other_jdot_gr => null()

      procedure(other_jdot_interface), pointer, nopass :: &
          other_jdot_ml => null()

      procedure(other_jdot_interface), pointer, nopass :: &
          other_jdot_tide => null()

      procedure(other_jdot_interface), pointer, nopass :: &
          other_extra_jdot => null()

      procedure(other_jdot_interface), pointer, nopass :: &
          other_jdot_ls => null()

      ! info for binary_history.data
      integer, pointer :: history_column_spec(:)
      logical, pointer :: history_column_is_integer(:) ! true if integer value; false if double
      real(dp), pointer :: history_column_values(:)

      logical :: doing_first_model_of_run = .false.
