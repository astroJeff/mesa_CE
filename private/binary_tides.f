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
 
      module binary_tides

      use star_lib
      use star_def
      use const_def
      use utils_lib
      use crlibm_lib
      use binary_def
      
      implicit none


      contains
      
      
      subroutine synch_spin_orbit_torque(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: osep ! orbital separation (cm)
         real(dp) :: qratio ! mass_other_star/mass_this_star
         real(dp) :: rlr ! roche lobe radius (cm)
         real(dp) :: dt_next ! next timestep
         real(dp) :: fsync, spin_period
            ! efficiency of tidal synchronization. (time scale × FSYNC ). 
         integer :: isb ! synchronization flag
            ! 0: no synchronisation 
            ! 1: synchronisation on Zahn’s timescale 
            ! 2: instantaneous synchronisation (e.g. for initializing a binary) 
            ! 3: synchronization on the time scale of the orbital period,
                  ! only for the roche-lobe filling star 
            ! 4: synchronisation on Tassoul’s timescale 
            ! 5: Zahn’s synchronisation timescale for radiative envelopes
         type (binary_info), pointer :: b
         ierr = 0
         if (ierr /= 0) return
         call binary_ptr(b, ierr)

         if (.not. (b% do_tidal_synch .or. &
             (b% doing_first_model_of_run .and. b% do_initial_orbit_synch))) &
             return

         call star_ptr(id, s, ierr)

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(s)) then
            rlr = b% rl(b% d_i)
         else
            qratio = 1.0/qratio
            rlr = b% rl(b% a_i)
         end if
         dt_next = s% dt
         fsync = 1.0
         if (b% doing_first_model_of_run .and. b% do_initial_orbit_synch) then
             isb = 2
         !else if (b% mtransfer_rate /= 0) then
         !   isb = 3
         else
            isb = 1
         end if

         call synch_spin_to_orbit(s, osep, qratio, rlr, dt_next, fsync, isb, ierr)
         
      end subroutine synch_spin_orbit_torque
      
      subroutine synch_spin_to_orbit(s, osep, qratio, rl, dt_next, fsync, isb, ierr)
         ! initially based on spiba.f kindly provided by Norbert Langer and group.
         use mlt_def, only: convective_mixing
         type (star_info), pointer :: s
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
         real(dp), intent(in) :: rl ! roche lobe radius (cm)
         real(dp), intent(in) :: dt_next ! next timestep
         real(dp), intent(in) :: fsync
         integer, intent(in) :: isb ! synchronization flag
         integer, intent(out) :: ierr
      
         real(dp) :: G, m, porb, m_lim, tsyn, r_phot, xig, xc, dalg, two_53, &
            tdyn, tkh, rho_face, cv_face, T_face, csound_face, ff, omega_orb, xc8
         real(dp), dimension(s% nz) :: j_synch, delta_j, tdyn_div_tkh
         integer :: k, nz, k_rl, k_xm
      
         include 'formats'
      
         ierr = 0
         if (osep <= 0) return
         if (qratio <= 0) return
         if (rl <= 0) return
         if (dt_next <= 0) return
         if (fsync <= 0) return
         
         
         
         k_rl = 0 ! BP: compiler warned may be used uninitialized
         k_xm = 0 ! BP: compiler warned may be used uninitialized
         
      
         nz = s% nz
         m = s% m(1)
         G = s% cgrav(1)
         r_phot = s% photosphere_r*Rsun
      
         ! porb = orbital period (seconds)
         porb = sqrt(osep*osep*osep*4*pi*pi/(G*m*(1+qratio)))
         ! j_synch = synchronized specific angular momentum
         dalg = 0
         omega_orb = 2*pi/porb
         do k=1,nz
            j_synch(k) = omega_orb*s% i_rot(k)
            dalg = dalg + abs(s% j_rot(k) - j_synch(k))
            delta_j(k) = 0
         end do
      
         if (isb == 0) then ! no synchronisation
         else if (isb == 2) then ! instantaneous synchronisation
            do k=1,nz
               delta_j(k) = s% j_rot(k) - j_synch(k)
            end do
         else if (isb == 3) then ! synch on timescale of orbital period
            if (r_phot >= rl) then ! filling roche lobe
               do k=1,nz
                  if (s% r(k) <= rl) then
                     k_rl = k
                     exit
                  end if
               end do
               m_lim = (1d0 - 1d-4)*s% m(k_rl)
               do k=1,nz
                  if (s% m(k) <= m_lim) then
                     k_xm = k
                     exit
                  end if
               end do
               tsyn = porb ! synchronize on timescale of orbital period
               do k=1,k_xm
                  delta_j(k) = (1d0 - exp_cr(-dt_next/tsyn))*(s% j_rot(k) - j_synch(k))
               end do
               do k=k_xm+1,nz
                  delta_j(k) = 0
               end do
            end if
         else ! isb = 1, 4, or 5
            ! get synchronization timescale
            if (isb == 1) then ! Zahn 1977, convective envelope
               tsyn = secyer*pow6(osep/r_phot)/qratio**2
            else if (isb == 4) then ! Tassoul & Tassoul
               tsyn = 1.44*4d1/(qratio*pow_cr(1d0+qratio,3d0/8d0))
               tsyn = tsyn*pow_cr(3.83d33/s% L(1),1d0/4d0)
               tsyn = tsyn*pow_cr(Msun/m,1d0/8d0)
               tsyn = tsyn*pow_cr(r_phot/Rsun,9d0/8d0)
               tsyn = tsyn*pow_cr(osep/r_phot,3.3d1/8d0)
               tsyn = tsyn*secyer
            else if (isb == 5) then ! Zahn 1977, radiative envelope
               xig = 0
               do k=2,nz
                  xig = xig + s% r(k)**2*(s% dm(k-1) + s% dm(k))/2
               end do
               xc = 0
               do k=10,nz
                  if (s% mixing_type(k) == convective_mixing .and. &
                      s% rho(k) > 1d5*s% rho(1)) then
                     xc = s% r(k)/r_phot
                     exit
                  end if
               end do
               two_53 = 3.1748021039364d0 ! 2d0**(5d0/3d0)
               if (xc <= 0) xc = 1d-7
               tsyn = 5d0*two_53*sqrt(G*m/(r_phot*r_phot*r_phot))
               tsyn = tsyn*qratio**2*pow_cr(1d0+qratio,5d0/6d0)
               tsyn = tsyn*(m*r_phot**2/xig)
               xc8 = xc*xc*xc*xc*xc*xc*xc*xc
               tsyn = tsyn*exp10_cr(-1.37d0)*xc8
               tsyn = tsyn*pow_cr(r_phot/osep,17d0/2d0)
               tsyn = 1d0/tsyn
               if (tsyn > 1d13*secyer) tsyn = 1d13*secyer
            else
               ierr = -1
               write(*,2) 'bad isb value for synch_spin_to_orbit', isb
               return
            end if
            !do k=1,nz
            !   delta_j(k) = (1d0 - exp_cr(-dt_next/tsyn))*(s% j_rot(k) - j_synch(k))
            !end do
            ! set local timescale following Gautschy & Glatzel 1990 MNRAS 245, 597-613.
            !write(*,*) "tsyn is", log10(tsyn/secyer)
            do k=1,nz
               !tdyn_div_tkh := local dynamical time-scale / local thermal time-scale
               rho_face = star_interp_val_to_pt(s% rho,k,nz,s% dq,'binary_tides')
               cv_face = star_interp_val_to_pt(s% cv,k,nz,s% dq,"binary_tides")
               T_face = star_interp_val_to_pt(s% T,k,nz,s% dq,"binary_tides")
               csound_face = star_interp_val_to_pt(s% csound,k,nz,s% dq,"binary_tides")
               tkh = 4*pi*s% r(k)**2*rho_face*cv_face*T_face/s% L(k) ! (4.4)
               tdyn = 1/csound_face
               tdyn_div_tkh(k) = tdyn/tkh
            end do
         
            dalg = max(0d0,dalg*(1d0 - exp_cr(-dt_next/(tsyn*fsync))))
            ! find ff to give desired value of dalg
            ff = rt(1d-30, 1d30, 5d-2, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'failed in tide routine'
               return
            end if
      
            do k=1,nz
               delta_j(k) = s% j_rot(k) - j_synch(k)
               if (delta_j(k) >= 0) then
                  delta_j(k) = min(delta_j(k),delta_j(k)*tdyn_div_tkh(k)**2*ff)
               else
                  delta_j(k) = max(delta_j(k),delta_j(k)*tdyn_div_tkh(k)**2*ff)
               end if
            end do
      
         end if

         do k=1,nz
            s% extra_jdot(k) = -delta_j(k)/dt_next
            s% extra_omegadot(k) = 0
         end do
   
         contains
      
      
         real(dp) function ffsync(ff)
            real(dp), intent(in) :: ff
            integer :: k
            real(dp) :: dalgtmp, dal
            dalgtmp = 0
            do k=1,nz
               dal = s% j_rot(k) - j_synch(k)
               if (dal >= 0) then
                  dal = min(dal,dal*tdyn_div_tkh(k)**2d0*ff)
               else
                  dal = max(dal,dal*tdyn_div_tkh(k)**2d0*ff)
               end if
               dalgtmp = dalgtmp + abs(dal)
            end do
            ffsync = dalg - dalgtmp
         end function ffsync


         real(dp) function rt(x1,x2,xacc,ierr)
            real(dp), intent(in) :: x1, x2, xacc
            integer, intent(out) :: ierr
         
            real(dp) :: dx, f, fmid, xmid
            integer :: j, jmax, k
            include 'formats'
            ierr = 0
            fmid = ffsync(x2)
            if (fmid == 0) then
               rt = x2; return
            end if
            f = ffsync(x1)
            if (f == 0) then
               rt = x1; return
            end if
            if (f*fmid > 0) then
               ierr = -1
               return
            end if
            if (f < 0) then
               rt = x1
               dx = x2 - x1
            else
               rt = x2
               dx = x1 - x2
            end if
            jmax = 500
            do j=1,jmax
               dx = dx/2
               xmid = rt + dx
               fmid = ffsync(xmid)
               if (fmid <= 0) rt = xmid
               if (abs(fmid/dalg) < xacc) return
            end do
            ierr = -1
         end function rt
   
   
      
      end subroutine synch_spin_to_orbit



      end module binary_tides
