
      ! these are the basic parameters of the run
      real(dp) :: m1 ! mass of first star
      real(dp) :: m2 ! mass of second star, or point mass
      real(dp) :: initial_period_in_days
      real(dp) :: initial_separation_in_Rsuns
      ! if initial_period_in_days <= 0, then use this instead to set up
 


      ! ADDED CONTROLS FOR COMMON ENVELOPE
      real(dp) :: alpha_CE
      logical :: do_CE

        
      ! other controls
      real(dp) :: fm
      real(dp) :: fm_limit
      real(dp) :: fa
      real(dp) :: fr
      real(dp) :: fr_limit
      real(dp) :: fj
      real(dp) :: dt_softening_factor
      real(dp) :: cur_mdot_frac
      real(dp) :: smallest_mass
      real(dp) :: max_abs_mdot ! in Msun/year units. e.g., 1d-6
      real(dp) :: max_yrs_dt
      real(dp) :: jdot_multiplier
      real(dp) :: upper_limit_on_period_in_hours
      real(dp) :: lower_limit_on_period_in_hours ! stop if period (in hours) is < this
      real(dp) :: factor_for_contact_terminate ! stop if accretor radius > this * its roche lobe
      
      real(dp) :: col_depth_for_eps_extra ! g cm^-2
      real(dp) :: irrad_flux_at_std_distance ! erg s^-1 cm^-2
      real(dp) :: std_distance_for_irradiation ! cm
      
      real(dp) :: magnetic_braking_gamma ! = 3.0 for normal case
      real(dp) :: alpha ! for mass loss jdot
      real(dp) :: companion_ratio_limit ! for limiting accretion as accretor fills its rl
      real(dp) :: transfer_fraction ! note: this is Saul's "beta"
      real(dp) :: tidal_Q ! if this is <= 0, then skip tidal braking. 
         ! e.g. 1e8 for hot Jupiter
      real(dp) :: R_companion ! photosphere radius (Rsun units)
         ! need this for calculating jdot_tide
         ! jdot_tide/J = -9/2(2 pi/Porb)*(Mp/Ms)*(1/Q)*(R_companion/separation)^5
         ! if set this to <= 0, use companion_mass**0.78d0 as estimate

      logical :: do_rotation, do_initial_orbit_synch, do_tidal_synch, do_j_accretion

      logical :: do_companion_ratio_limit ! limit accretion as accretor fills its rl

      logical :: do_jdot_mb
      logical :: do_jdot_gr
      logical :: do_jdot_ml
      logical :: do_jdot_tide
      logical :: do_jdot_ls

      logical :: use_other_jdot_mb
      logical :: use_other_jdot_gr
      logical :: use_other_jdot_ml
      logical :: use_other_jdot_tide
      logical :: use_other_extra_jdot
      logical :: use_other_jdot_ls

      logical :: keep_donor_fixed
      
      logical :: accretion_powered_irradiation
      real(dp) :: accretor_radius_for_irrad ! (cm)
      real(dp) :: max_F_irr ! (erg cm^-2 s^-1)

      logical :: include_accretor_mb
      logical :: limit_retention_by_mdot_edd 
      real(dp) :: use_this_for_mdot_edd

      character (len=strlen) :: mdot_scheme

      ! controls for implicit selection of mdot to get required RL relative overlap at end of step
      real(dp) :: rl_rel_overlap_tolerance ! e.g. 1d-3
         ! if at end of step, donor R >= RL, then require (R-RL)/RL <= rl_rel_overlap_tolerance.
            ! if relative overlap is too large, then select larger mdot and redo the step.
         ! if at end of step, donor R < RL, then require mdot = 0.
            ! if mdot is not 0, then select smaller mdot, perhaps 0, and redo the step.
      integer :: max_tries_to_achieve ! <= 0 means skip this process
         ! if this is > 0 and fail to achieve desired rl_rel_overlap, retry with smaller timestep.
      real(dp) :: initial_change_factor ! e.g. 2.0
         ! amount by which the current mass transfer rate will be changed during a step in
         ! the implicit scheme. can be automatically adjusted during the run if
         ! implicit_lambda /= 0
      real(dp) :: implicit_lambda ! e.g. 0.25
      real(dp) :: change_factor_fraction ! e.g 0.5
         ! each step during mass transfer the change factor will be modified according to:
         !    if(mass_transfer_rate < mass_transfer_prev) then
         !       change_factor = change_factor*(1.0-implicit_lambda) &
         !          + implicit_lambda*(1+change_factor_fraction*(mass_transfer_rate/mass_transfer_prev-1))
         !    else
         !       change_factor = change_factor*(1.0-implicit_lambda) &
         !          + implicit_lambda*(1+change_factor_fraction*(mass_transfer_prev/mass_transfer_rate-1))
         !    end if
         ! choosing implicit_lambda = 0 will keep the change factor constant.
      real(dp) :: max_change_factor
      real(dp) :: min_change_factor

      real(dp) :: starting_mdot ! starting guess for mdot used in the implicit scheme
      real(dp) :: implicit_min_mdot ! if mtransfer_rate goes below this, assume no mass transfer
      
      ! binary_history.data controls
      character (len=strlen) :: history_name
      character (len=strlen) :: log_directory
      character (len=strlen) :: history_dbl_format
      character (len=strlen) :: history_int_format
      character (len=strlen) :: history_txt_format

      integer :: photostep
      integer :: photo_digits
         
      ! end of controls
