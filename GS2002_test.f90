module parameters

  integer, parameter :: k_cbm = 0 ! 0 for ISM, 2 for wind-like CBM
  real(kind=8), parameter :: n_cbm = 1.d0  ! Density if k = 0 (#/cm^3)
  real(kind=8), parameter :: A_star = 3.d33 ! Density if k = 2 (#/cm)
  
  real(kind=8), parameter :: E_iso = 1.d52   ! Impulsive energy added
  real(kind=8), parameter :: Gam_coast = 1000 ! E/m for ejecta and also peak
                                             !   Lorentz factor of fireball

  integer, parameter :: n_tobs = 51  ! Number of t_obs values to use,
                                     !   *after* tobs_start
  real(kind=8), parameter :: tobs_start = 1.d1  ! Start and end times to
  real(kind=8), parameter :: tobs_end   = 1.195d6  !   consider, in Earth frame
  
  integer, parameter :: n_ph  = 640  ! Size of photon energy array
  
    ! Luminosity distance and redshift. Only set one of these, and let
    !   subroutine cosmo_calc handle the other
  real(kind=8) :: dist_lum
  real(kind=8) :: redshift = 1.d0
  
  real(kind=8), parameter :: p_spec = 2.23d0 ! Spectral index of electron
                                            !   distribution before cooling
  
  real(kind=8), parameter :: epse = 0.1d0   ! Fraction of energy in electrons
  real(kind=8), parameter :: epse_bar = (p_spec-2.d0)/(p_spec-1.d0)*epse
                                ! And a form more convenient for manipulating
  
  real(kind=8), parameter :: epsB0 = 0.01d0 ! Fraction of energy in B fields
  
  logical, parameter :: do_ssa = .true.  ! Calculate synch self-absorption?
  
end module parameters

module constants

  real(kind=8), parameter :: third = 1.d0/3.d0
  real(kind=8), parameter :: pii = 2.d0*asin(1.d0)
  real(kind=8), parameter :: degtrd = pii/180.0  ! Degrees to radians
  real(kind=8), parameter :: radtdg = 1.0d0/degtrd ! Radians to degrees
  
  real(kind=8), parameter :: xmp = 1.6726218d-24    ! proton mass in grams
  real(kind=8), parameter :: xme = 9.10938291d-28    ! electron mass
  real(kind=8), parameter :: ccgs = 2.99792458d10   ! speed of light (cm/s)
  real(kind=8), parameter :: qcgs= 4.80320451d-10     ! proton charge in cgs units:
                                                 !  _E_lectro_S_tatic _U_nities
  real(kind=8), parameter :: sigT = 6.652458d-25 ! Thomson cross section (cm^-2)
  real(kind=8), parameter :: xkb  = 1.3806488d-16    ! Boltzmann's constant (cgs)
  real(kind=8), parameter :: xh   = 6.626070d-27     ! Planck constant (erg-sec)
  real(kind=8), parameter :: xh_bar = xh/(2.0d0*pii)
  real(kind=8), parameter :: rm_prot = xmp*(ccgs**2) ! proton rest mass en. [erg]
  real(kind=8), parameter :: rm_elec = xme*(ccgs**2) ! electron rest mass en. [erg]
  real(kind=8), parameter :: o_o_mec = 1.d0/(xme*ccgs)
  real(kind=8), parameter :: o_o_rme = 1.d0/rm_elec
  
  real(kind=8), parameter :: pc_to_cm = 3.084d18 ! conversion parsec to cm
  real(kind=8), parameter :: ergtev = 6.242d11  ! conversion ergs to eV
  real(kind=8), parameter :: etkev  = 6.242d08  ! conversion ergs to keV
  real(kind=8), parameter :: etmev  = 6.242d05  ! conversion ergs to MeV
  real(kind=8), parameter :: evterg = 1.602d-12 ! conversion eV to ergs
  real(kind=8), parameter :: xkevte = 1.602d-09 ! conversion keV to ergs
  real(kind=8), parameter :: xmevte = 1.602d-06 ! conversion MeV to ergs
  
end module constants

module data_arrays
  use parameters

  ! Data used for computing Psyn_array
  real(kind=8), dimension(0:25) :: y_array
  real(kind=8), dimension(0:40) :: chi_array
  real(kind=8), dimension(0:n_tobs) :: tobs_array
  real(kind=8), dimension(0:n_ph) :: phot_en_array, tzero_array
  real(kind=8), dimension(0:n_ph,0:n_ph,0:40) :: Psyn_array

  ! Data used in photon production
  real(kind=8) :: o_o_sone
  real(kind=8), dimension(0:n_ph) :: syn_F, syn_x
  
  ! State of the afterglow at a particular t_obs
  real(kind=8) :: Gam_LOS, R_LOS, t_LOS
  
  ! Values used in Euler method subroutine
  integer :: em_steps
  real(kind=8), dimension(n_ph) :: em_size_array
  
end module data_arrays

!****************************************************************************
!****************************************************************************
! MAIN PROGRAM
!****************************************************************************
!****************************************************************************
program GS2002_test

! Evaluates afterglow emission using three methods from Granot & Sari
!   (2002ApJ...568..820G): Eqs (5)-(9), Eq (A14), and Eq (A24).
!
! Ignores synchrotron self-absorption (i.e. alpha_nu = 0 in Eq A24).

use parameters
use constants
use data_arrays

implicit none

interface
  double precision function GPS_egg(y_input, i_params, r_params)
    real(kind=8), intent(in) :: y_input
    integer, intent(in) :: i_params(7)
    real(kind=8), intent(in) :: r_params(7)
  end function

  subroutine rad_transfer_diffeq(s, n_phot, jnu, alpha, r_params, i_params)
    integer, intent(in) :: n_phot
    real(kind=8), intent(in) :: s
    integer, intent(in) :: i_params(7)
    real(kind=8), intent(in) :: r_params(7)
    real(kind=8), intent(out) :: jnu(n_phot), alpha(n_phot)
  end subroutine
end interface  

logical :: fast_cooling
integer :: i_tm, j_ph, k_chi, l_y, k_x
real(kind=8) :: t_fac, dist_lum_Mpc, Psyn, y_integral_running, y_hi,      &!&
     y_integrand_upper, y_integrand_lower, y_in, y_lo,                    &!&
     chi_integral_running, chi_integrand_lower, chi_integrand_upper,      &!&
     chi_hi, integrand, chi_in, chi_lo, t_d, nu_m, nu_c, F_ext_c, F_ext_m,&!&
     nu_phot, F11, F9_tilde, F2, F3_tilde, R_perp_max, y_crit,            &!&
     eqA22_split, y_max, y_min, s_far, s_near, nu_ac, nu_sa, F_ext_ac,    &!&
     F_ext_sa, F7, F10_tilde, F11_tilde, F1, F2_tilde
real(kind=8), dimension(0:n_ph) :: F_nu, I_nu_in, I_nu_out, I_nu_hi,      &!&
     I_nu_lo, nuFnu
integer, dimension(7) :: i_params
real(kind=8), dimension(7) :: r_params


  ! Pre-compute F(x), which will save many, many calls to the Bessel function
  !   in the loops to follow
  call calc_syn_Fx(syn_F, syn_x)
  o_o_sone = 1.d0 / syn_x(1)
  
  
  ! Fill tobs_array
  t_fac = (tobs_end/tobs_start)**(1.d0/float(n_tobs))
  tobs_array(0) = tobs_start
  do i_tm = 1, n_tobs
    tobs_array(i_tm) = tobs_array(i_tm-1) * t_fac
  enddo
  
  ! Create y_array and chi_array
  call set_y_array(y_array)
  call set_chi_array(chi_array)
  call set_tzero_array(tzero_array)
  call set_phot_en_array(phot_en_array)
  
  
  ! Find whichever of redshift or comoving distance was not specified; must
  !   convert dist_lum to megaparsecs since that's what cosmo_calc expects
  if( ((dist_lum .eq. 0.d0) .and. (redshift .gt. 0.d0)) .or.              &!&
      ((dist_lum .gt. 0.d0) .and. (redshift .eq. 0.d0)) ) then
    dist_lum_Mpc = dist_lum / (pc_to_cm*1.d6)
    call cosmo_calc(dist_lum_Mpc, redshift)
    if( dist_lum .eq. 0.d0 ) dist_lum = dist_lum_Mpc * (1.d6*pc_to_cm)
  else
    write(*,"(2A)") 'ERROR: set one of dist_lum/redshift to zero, and ',  &!&
         'the other to a positive number.'
    write(*,"(A)") 'Stopping program now'
    stop
  endif
  
  
  ! Pre-compute Psyn_array, which can be reduced to a 3-D array:
  !   (1) nu
  !   (2) Gam0 (conditions at original shock crossing)
  !   (3) chi
!  do j_ph = 0, n_ph
!    do i_tm = 0, n_ph
!      do k_chi = 0, 40
!
!        call compute_Psyn(j_ph, i_tm, k_chi, Psyn)
!        Psyn_array(j_ph, i_tm, k_chi) = Psyn
!        
!      enddo  ! Loop over k_chi
!    enddo  ! Loop over i_tm
!    ! Progress printout
!    if( mod(j_ph, 10) .eq. 0 ) write(*,"(2X,A,I0)") 'j_ph = ',j_ph
!  enddo  ! Loop over j_ph
  
  
  ! For each t_obs and each (observer-frame) frequency, evaluate Eq (A14).
  !   This is a double integral over y and chi, which we will handle using
  !   the pre-computed y_array and chi_array values and the trapezoid rule.
  do i_tm = 0,-1!0, n_tobs
    
    ! Compute Gam_LOS and R_LOS for this value of t_obs, since they're the
    !   same at every point during the inner loops and there's no need to
    !   waste computation cycles re-computing them
    ! R_LOS from Eq (A15)
    R_LOS = (17.d0-4.d0*k_cbm) * (4.0-k_cbm) * E_iso * tobs_array(i_tm)   &!&
           / (1.d0+redshift) / (4.d0 * pii * xmp * ccgs)
    if( k_cbm .eq. 0 ) then
      R_LOS = ( R_LOS / n_cbm )**(1.d0/(4.d0-k_cbm))
    else
      R_LOS = ( R_LOS / A_star )**(1.d0/(4.d0-k_cbm))
    endif
    
    ! Gam_LOS from Eq (A15)
    Gam_LOS = (17.d0-4.d0*k_cbm) * E_iso                                  &!&
             / ( 4.d0**(5.d0-k_cbm) * (4.d0-k_cbm)**(3.d0-k_cbm) * pii    &!&
                * rm_prot * ( ccgs*tobs_array(i_tm)                       &!&
                             / (1.d0+redshift) )**(3.d0-k_cbm) )
    if( k_cbm .eq. 0 ) then
      Gam_LOS = sqrt(2.d0) * (Gam_LOS / n_cbm)**(0.5d0/(4.d0-k_cbm))
    else
      Gam_LOS = sqrt(2.d0) * (Gam_LOS / A_star)**(0.5d0/(4.d0-k_cbm))  
    endif
    
    
    ! Write the identifier for this time step to the output file
    write(800,"(A,I3,A,ES12.3E2)") 'i_tm = ', i_tm, ';  t_obs = ',        &!&
       tobs_array(i_tm)
    
    
    do j_ph = 0, n_ph
      
      ! Outer integral: over y.  We work outside-in, taking advantage of two
      !   approximations.  First, when y = 1 the integrand evaluates to 0 --
      !   emission drops to zero precisely at the rim.  Second, emissivity is
      !   almost flat at the center of the image, which means that the value
      !   at y_array(0) is almost equal to the value at y_array(1)
      y_integral_running = 0.d0
      y_hi = 1.d0
      
      y_integrand_upper = 0.d0
      y_integrand_lower = 0.d0
      do l_y = 24, 1, -1
        
        y_in = y_array(l_y)
        y_lo = y_in
        
        ! Inner integral: over chi.  Again, work outside in.  First point of
        !   evaluation is set by y_in, and we will work all the way down to
        !   chi = 1.
        chi_integral_running = 0.d0
        
        chi_integrand_upper = 0.d0
        chi_integrand_lower = 0.d0
        chi_hi = y_in**(k_cbm-4.d0)
        ! Evaluate integrand at chi_hi, then start rest of loop
        call evaluate_A14_integrand(i_tm, j_ph, 40, y_in, chi_hi,         &!&
           integrand)
        chi_integrand_upper = integrand
        
        do k_chi = count(chi_array .le. chi_hi)-1, 0, -1
          
          ! Evaluate next lower value on chi_array, and use that to update
          !   the running total using the trapezoid rule
          chi_in = chi_array(k_chi)
          chi_lo = chi_in
          call evaluate_A14_integrand(i_tm, j_ph, k_chi, y_in, chi_in,    &!&
             integrand)
          chi_integrand_lower = integrand
          
          chi_integral_running = chi_integral_running                     &!&
                                + 0.5d0 * ( chi_integrand_upper           &!&
                                           + chi_integrand_lower )        &!&
                                        * ( chi_hi - chi_lo )
          
          ! Set up for next pass through loop
          chi_hi = chi_lo
          chi_integrand_upper = chi_integrand_lower
          
        enddo  ! loop over chi values
        
        
        ! Now have y_integrand_lower (= chi_integral_running).  Use that to
        !   update y_integral_running
        y_integrand_lower = chi_integral_running
        y_integral_running = y_integral_running                           &!&
                            + 0.5d0 * ( y_integrand_upper                 &!&
                                       + y_integrand_lower )              &!&
                                    * ( y_hi - y_lo )
        
        ! Set up for next pass through loop
        y_hi = y_lo
        y_integrand_upper = y_integrand_lower
        
      enddo  ! loop over y values
      
      ! With the loop finished, take care of the final y value: 0.  Figure
      !   11 of GPS99 suggests that the integral is flat in this part of the
      !   afterglow, so just reuse the preceding y_integrand value
      y_integral_running = y_integral_running                             &!&
                          + 0.5d0 * ( y_integrand_lower                   &!&
                                     + y_integrand_lower )                &!&
                                  * ( y_hi - y_array(0) )
      
      ! Multiply by the prefactor from GS2002 Eq (A14)
      y_integral_running = y_integral_running                             &!&
                          * 2.d0*(4.d0-k_cbm) * R_LOS**3                  &!&
                          * (1.d0 + redshift) / dist_lum**2
      if( y_integral_running .lt. 1.d-55 ) y_integral_running = 1.d-99
      
      
      ! Write the result to file
      write(800,"(2i4,2es13.3e2,es13.4e2,es13.3e2)") i_tm, j_ph,          &!&
         log10(phot_en_array(j_ph)*ergtev),                               &!&
         log10(phot_en_array(j_ph)/xh),                                   &!&
         log10(y_integral_running * phot_en_array(j_ph)/xh),    &!& nu*Fnu
         log10(y_integral_running * 1.d23)                       ! Fnu (Jy)
      
    enddo  ! loop over photon energies
    
    ! Write a line to explain the columns
    write(800,"(12X,A)") 'log10(eV)    log10(Hz)  log10(nuFnu)   log10(Jy)'
    ! Skip a line before next time step
    write(800,*)
    
    ! Progress meter
    if(mod(i_tm, 5).eq.0) write(*,"(A,I0)") 'i_tm = ',i_tm
    
  enddo  ! loop over observer times


   !-------------------------------------------------------------------------
!--! For each time step, fill out the flux array nu F_nu for each entry in
   !   the array phot_en_array.
   ! We will need different computations based on whether the spectrum is
   !   fast- or slow-cooling.
   ! In the slow-cooling case, breaks (1), (2), and (3) separate spectral
   !   regions B, D, G, and H respectively.  These breaks are called nu_sa,
   !   nu_m, and nu_c below due to their physical origins: synchtroton
   !   self-absorption, the minimum electron energy, and the cooling break.
   ! In the fast-cooling case, breaks (7), (10), (11), and (9) separate
   !   spectral regions B, C, E, F, and H respectively.  The breaks are
   !   nu_ac, nu_sa, nu_c, and nu_m.  The new break, nu_ac, is synchrotron
   !   self-absorption by electrons that have fully cooled (Absorption by
   !   Cooled electrons).
   !
   ! As stated in GS2002, there are three conditions we could use to decide
   !   whether we're in fast- or slow-cooling:
   !  (1) nu_2 = nu_3  (nu_m = nu_c in the slow-cooling case)
   !  (2) nu_7 = nu_10 (nu_ac = nu_sa, implying no electron cooling)
   !  (3) nu_9 = nu_11 (nu_m = nu_c in the fast-cooling case)
   ! We'll use option (1) since, per GS2002, nu_10 and nu_11 cannot be
   !   calculated in the k = 2 scenario.  When k = 2, we will assume that
   !   nu_11 = nu_3 and nu_10 = nu_1.
   !
   ! Only worry about self-absorption, though, if instructed to do so by the
   !   do_ssa parameter.
   !-------------------------------------------------------------------------
  do i_tm = 0, n_tobs
    
    t_d = tobs_array(i_tm) / 86400.d0   ! t_obs in days
    
    
    ! Compute the locations of the nu_m and nu_c breaks using Table 2 of
    !   GS2002.
    if( k_cbm .eq. 0 ) then
      nu_m = 3.73d15 * (p_spec - 0.67d0) * sqrt(1.d0+redshift)            &!&
            * sqrt(E_iso/1.d52) * sqrt(epsB0) * epse_bar**2 * t_d**(-1.5d0)
      nu_c = 6.37d13 * (p_spec - 0.46d0) * exp(-1.16d0*p_spec)            &!&
            / sqrt(1.d0+redshift) * epsB0**(-1.5d0) / n_cbm               &!&
            / sqrt(E_iso/1.d52) / sqrt(t_d)
    else
      nu_m = 4.02d15 * (p_spec - 0.69d0) * sqrt(1.d0+redshift)            &!&
            * sqrt(E_iso/1.d52) * epse_bar**2 * sqrt(epsB0) * t_d**(-1.5d0)
      nu_c = 4.40d10 * (3.45d0 - p_spec) * exp(0.45d0*p_spec)             &!&
            * (1.d0+redshift)**(-1.5d0) * epsB0**(-1.5d0)                 &!&
            * sqrt(E_iso/1.d52) * (A_star*xmp/5.d11)**(-2.d0) * sqrt(t_d)
    endif

    
    ! Decide whether we're in fast or slow cooling.  If fast cooling, we will
    !   need to recalculate nu_m and nu_c -- unless k = 2, in which case we
    !   only recalculate nu_m (we assume nu_11 = nu_3)
    ! Note that F_ext quantities are given in units of Jy, not mJy as in
    !   GS2002
    if( nu_m .gt. nu_c ) then
      
      fast_cooling = .true.
      
      if( k_cbm .eq. 0 ) then
        
        nu_m = 3.94d15 * (p_spec - 0.74d0) * sqrt(1.d0+redshift)          &!&
              * epse_bar**2 * sqrt(epsB0) * sqrt(E_iso/1.d52) * t_d**(-1.5d0)
        nu_c = 5.86d12 / sqrt(1.d0+redshift) * epsB0**(-1.5d0) / n_cbm    &!&
              / sqrt(E_iso/1.d52) / sqrt(t_d)
        
        F_ext_c = 2.84d-2 * (1.d0+redshift) * sqrt(epsB0) * sqrt(n_cbm)   &!&
                 * (E_iso/1.d52) * (dist_lum/1.d28)**(-2)
        
        if( do_ssa ) then
          
          nu_ac = 1.12d8 * ( (3.d0*p_spec-1)/(3.d0*p_spec+2) )**1.6d0     &!&
                 * (1.d0+redshift)**(-1.3d0) * epse_bar**(-1.6d0)         &!&
                 * epsB0**(-0.4d0) * n_cbm**0.3d0                         &!&
                 * (E_iso/1.d52)**(-0.1d0) * t_d**0.3d0
          nu_sa = 1.32d10 / sqrt(1.d0+redshift) * epsB0**1.2d0            &!&
                 * n_cbm**1.1d0 * (E_iso/1.d52)**0.7d0 / sqrt(t_d)
          
          F_ext_ac = 5.27d-6 * ( (3.d0*p_spec-1)/(3.d0*p_spec+2) )**2.2d0 &!&
                    * (1.d0+redshift)**(-0.1d0) * epse_bar**(-2.2d0)      &!&
                    * epsB0**(-0.8d0) * n_cbm**0.1d0                      &!&
                    * (E_iso/1.d52)**0.3d0 * t_d**1.1d0                   &!&
                    / (dist_lum/1.d28)**2
          
        endif
        
      else
        
        nu_m = 3.52d15 * (p_spec - 0.31d0) * sqrt(1.d0 + redshift)        &!&
              * epse_bar**2 * sqrt(epsB0) * sqrt(E_iso/1.d52) * t_d**(-1.5d0)
        ! nu_c is same as for slow cooling, for lack of a better option
        
        ! F_ext_c will be same as slow cooling, for lack of better option
        F_ext_c = 8.02d2 * exp(7.02d0 * (p_spec-2.5d0))                   &!&
                 * (1.d0+redshift)**(p_spec + 0.5d0)                      &!&
                 * epse_bar**(p_spec - 1.d0)                              &!&
                 * epsB0**(p_spec - 0.5d0)                                &!&
                 * (A_star*xmp/5.d11)**p_spec                             &!&
                 * sqrt(E_iso/1.d52)                                      &!&
                 * t_d**(0.5d0 - p_spec)                                  &!&
                 * (dist_lum/1.d28)**(-2)
        
        if( do_ssa ) then
          
          nu_ac = 1.68d8 * ( (3.d0*p_spec-1)/(3.d0*p_spec+2) )**1.6d0     &!&
                 / (1.d0+redshift) * epse_bar**(-1.6d0)                   &!&
                 * epsB0**(-0.4d0) * (A_star*xmp/5.d11)**0.6d0            &!&
                 * (E_iso/1.d52)**(-0.4d0)
          ! nu_sa is same as for slow cooling, for lack of a better option
          nu_sa = 8.31d9 * ( (p_spec-1)/(3.d0*p_spec+2) )**0.6d0          &!&
                 * (1.d0+redshift)**(-0.4d0) / epse_bar * epsB0**0.2d0    &!&
                 * (A_star*xmp/5.d11)**1.2d0 * (E_iso/1.d52)**(-0.4d0)    &!&
                 * t_d**(-0.6d0)
          
          F_ext_ac = 5.27d-6 * ( (3.d0*p_spec-1)/(3.d0*p_spec+2) )**2.2d0 &!&
                    * (1.d0+redshift)**(-0.1d0) * epse_bar**(-2.2d0)      &!&
                    * epsB0**(-0.8d0) * n_cbm**0.1d0                      &!&
                    * (E_iso/1.d52)**0.3d0 * t_d**1.1d0                   &!&
                    / (dist_lum/1.d28)**2
          
        endif
        
      endif
      
    else
      
      fast_cooling = .false.
      
      if( k_cbm .eq. 0 ) then
        
        F_ext_m = 9.93d-3 * (p_spec + 0.14d0) * (1.d0+redshift)           &!&
                 * sqrt(epsB0) * sqrt(n_cbm) * (E_iso/1.d52)              &!&
                 * (dist_lum/1.d28)**(-2)
        
        if( do_ssa ) then
          
          nu_sa = 1.24d9 * ( (p_spec-1)/(3.d0*p_spec+2) )**0.6d0          &!&
                 / (1.d0+redshift) / epse_bar * epsB0**0.2d0              &!&
                 * n_cbm**0.6d0 * (E_iso/1.d52)**0.2d0
          
          F_ext_sa = 6.47d-4 * (p_spec-1)**1.2d0 / (3.d0*p_spec-1)        &!&
                    / (3.d0*p_spec+2)**0.2d0 * sqrt(1.d0+redshift)        &!&
                    / epse_bar * epsB0**0.4d0 * n_cbm**0.7d0              &!&
                    * (E_iso/1.d52)**0.9d0 * sqrt(t_d)                    &!&
                    / (dist_lum/1.d28)**2
          
        endif
        
      else
        
        F_ext_m = 7.69e-2 * (p_spec + 0.12d0) * (1.d0+redshift)**1.5d0    &!&
                 * sqrt(epsB0) * (A_star*xmp/5.d11) * sqrt(E_iso/1.d52)   &!&
                 / sqrt(t_d) * (dist_lum/1.d28)**(-2)
        
        if( do_ssa ) then
          
          nu_sa = 8.31d9 * ( (p_spec-1)/(3.d0*p_spec+2) )**0.6d0          &!&
                 * (1.d0+redshift)**(-0.4d0) / epse_bar * epsB0**0.2d0    &!&
                 * (A_star*xmp/5.d11)**1.2d0 * (E_iso/1.d52)**(-0.4d0)    &!&
                 * t_d**(-0.6d0)
          
          F_ext_sa = 9.19d-3 * (p_spec-1)**1.2d0 / (3.d0*p_spec-1)        &!&
                    / (3.d0*p_spec+2)**0.2d0 * (1.d0+redshift)**1.2d0     &!&
                    / epse_bar * epsB0**0.4d0 * (A_star*xmp/5.d11)**1.4d0 &!&
                    * (E_iso/1.d52)**0.2d0 * t_d**(-0.2d0)                &!&
                    / (dist_lum/1.d28)**2
          
        endif
        
      endif
      
    endif  ! check on fast/slow cooling
    
    
    ! Equation 5 defines F_nu as the product of several terms.  Each of these
    !   terms depends on photon frequency, so place remaining calculations
    !   in a loop over photon energy.
    F_nu(:) = 1.d-99
    do j_ph = 0, n_ph
      
      nu_phot = phot_en_array(j_ph) / xh
      
      if( fast_cooling ) then
        
        if( .not. do_ssa ) then
          
          F11 = (nu_phot/nu_c)**(-0.597d0 * third)                        &!&
               + (nu_phot/nu_c)**(0.597d0 * 0.5d0)
          F11 = F_ext_c / F11**(1.d0/0.597d0)
          
          F9_tilde = 1.d0  +  (nu_phot/nu_m)**( (3.34d0 - 0.82d0*p_spec)  &!&
                                               * (-0.5d0 + 0.5d0*p_spec) )
          F9_tilde = 1.d0 / F9_tilde**(1.d0/(3.34d0 - 0.82d0*p_spec))
          
          F_nu(j_ph) = F11 * F9_tilde
          
        else
          
          F7 = (nu_phot/nu_ac)**(-(1.99d0-4.d-2*p_spec) * 2.d0)           &!&
              + (nu_phot/nu_ac)**(-(1.99d0-4.d-2*p_spec) * 1.375d0)
          F7 = F_ext_ac / F7**(1.d0/(1.99d0-4.d-2*p_spec))
          
          F10_tilde = 1.d0  +  (nu_phot/nu_sa)**( 1.213d0                 &!&
                                                 * (1.375d0 - third) )
          F10_tilde = 1.d0 / F10_tilde**(1.d0/1.213d0)
          
          F11_tilde = 1.d0  +  (nu_phot/nu_c)**( 0.597d0 * (third + 0.5d0) )
          F11_tilde = 1.d0 / F11_tilde**(1.d0/0.597d0)
          
          F9_tilde = 1.d0  +  (nu_phot/nu_m)**( (3.34d0 - 0.82d0*p_spec)  &!&
                                               * (-0.5d0 + 0.5d0*p_spec) )
          F9_tilde = 1.d0 / F9_tilde**(1.d0/(3.34d0 - 0.82d0*p_spec))
          
          F_nu(j_ph) = F7 * F10_tilde * F11_tilde * F9_tilde
          
        endif
        
      else
        
        if( .not. do_ssa ) then
          
          F2 = (nu_phot/nu_m)**( -(1.84d0 - 0.4d0*p_spec)*third )         &!&
              + (nu_phot/nu_m)**( -(1.84d0 - 0.4d0*p_spec)                &!&
                                   * (1.d0-p_spec)*0.5d0 )
          F2 = F_ext_m / F2**( 1.d0/(1.84d0 - 0.4d0*p_spec) )
          
          F3_tilde = 1.d0  +  (nu_phot/nu_c)**( (1.15d0-0.06d0*p_spec)    &!&
                                               * 0.5d0 ) ! (1-p)/2 - (-p/2)
          F3_tilde = 1.d0 / F3_tilde**(1.d0/(1.15d0-0.06d0*p_spec))
          
          F_nu(j_ph) = F2 * F3_tilde
          
        else
          
          F1 = (nu_phot/nu_sa)**(-1.64d0 * 2.d0)                          &!&
              + (nu_phot/nu_sa)**(-1.64d0 * third)
          F1 = F_ext_sa / F1**( 1.d0/1.64d0 )
          
          F2_tilde = 1.d0  +  (nu_phot/nu_m)**( (1.84d0-0.4d0*p_spec)     &!&
                                               * (third - 0.5d0*(1.d0-p_spec)) )
          F2_tilde = 1.d0 / F2_tilde**( 1.d0/(1.84d0 - 0.4d0*p_spec) )
          
          F3_tilde = 1.d0  +  (nu_phot/nu_c)**( (1.15d0-0.06d0*p_spec)    &!&
                                               * 0.5d0 ) ! (1-p)/2 - (-p/2)
          F3_tilde = 1.d0 / F3_tilde**(1.d0/(1.15d0-0.06d0*p_spec))
          
          F_nu(j_ph) = F1 * F2_tilde * F3_tilde
          
        endif
        
      endif
      
    enddo  ! loop over photon_energies
    
    
    ! Now convert F_nu into nu*F_nu, go from Jy to cgs, and write it out
    nuFnu(:) = F_nu(:) * phot_en_array(:) / xh * 1.d-23
    where( F_nu .lt. 1.d-55 ) F_nu = 1.d-99
    ! Write the identifier for this time step to the output file
    write(810,"(A,I3,A,ES12.3E2)") 'i_tm = ', i_tm, ';  t_obs = ',        &!&
       tobs_array(i_tm)
    do j_ph = 0, n_ph
      write(810,"(2i4,2es13.3e2,es14.5e2,es13.3e2)") i_tm, j_ph,    &!&
           log10(phot_en_array(j_ph)*ergtev),                       &!&
           log10(phot_en_array(j_ph)/xh),                           &!&
           log10(nuFnu(j_ph)),                                      &!&
           log10(F_nu(j_ph))
    enddo
    ! Explanation line
    write(810,"(12X,A)") 'log10(eV)    log10(Hz)   log10(nuFnu)   log10(Jy)'
    ! Skip a line before next time step
    write(810,*)
    
  enddo  ! loop over observer times
   !-------------------------------------------------------------------------
   ! Analytical solution calculated
   !-------------------------------------------------------------------------
  

end program GS2002_test


!****************************************************************************
!****************************************************************************
subroutine set_y_array(y_array)

! Initializes the values of the y_array, which is one of the arrays we will
!   use for applying Eqs (A14) and (A24)
!
! Output argument:
!  1) y_array: array of y values at which we will calculate emission

implicit none

  ! Output argument
real(kind=8), dimension(0:25), intent(out) :: y_array

  ! Set y_array values
  y_array( 0) = 0.0d0
  y_array( 1) = 0.1d0
  y_array( 2) = 0.2d0
  y_array( 3) = 0.3d0
  y_array( 4) = 0.4d0
  y_array( 5) = 0.5d0
  y_array( 6) = 0.6d0
  y_array( 7) = 0.65d0
  y_array( 8) = 0.70d0
  y_array( 9) = 0.75d0
  y_array(10) = 0.80d0
  y_array(11) = 0.82d0
  y_array(12) = 0.84d0
  y_array(13) = 0.86d0
  y_array(14) = 0.88d0
  y_array(15) = 0.90d0
  y_array(16) = 0.91d0
  y_array(17) = 0.92d0
  y_array(18) = 0.93d0
  y_array(19) = 0.94d0
  y_array(20) = 0.95d0
  y_array(21) = 0.96d0
  y_array(22) = 0.97d0
  y_array(23) = 0.98d0
  y_array(24) = 0.99d0
  y_array(25) = 1.00d0

return
end subroutine set_y_array


!****************************************************************************
!****************************************************************************
subroutine set_chi_array(chi_array)

! Initializes the values of the chi_array, which is one of the arrays we will
!   use for applying Eqs (A14) and (A24)
!
! Output argument:
!  1) chi_array: array of chi values at which we will calculate emission

implicit none

  ! Output argument
real(kind=8), dimension(0:40), intent(out) :: chi_array

  ! Set chi_array values
  chi_array( 0) = 1.d0
  chi_array( 1) = 1.000001d0
  chi_array( 2) = 1.000002d0
  chi_array( 3) = 1.000003d0
  chi_array( 4) = 1.000005d0
  chi_array( 5) = 1.00001d0
  chi_array( 6) = 1.00002d0
  chi_array( 7) = 1.00003d0
  chi_array( 8) = 1.00005d0
  chi_array( 9) = 1.0001d0
  chi_array(10) = 1.0002d0
  chi_array(11) = 1.0003d0
  chi_array(12) = 1.0005d0
  chi_array(13) = 1.001d0
  chi_array(14) = 1.002d0
  chi_array(15) = 1.003d0
  chi_array(16) = 1.005d0
  chi_array(17) = 1.01d0
  chi_array(18) = 1.02d0
  chi_array(19) = 1.03d0
  chi_array(20) = 1.05d0
  chi_array(21) = 1.1d0
  chi_array(22) = 1.2d0
  chi_array(23) = 1.3d0
  chi_array(24) = 1.5d0
  chi_array(25) = 2.d0
  chi_array(26) = 3.d0
  chi_array(27) = 5.d0
  chi_array(28) = 10.d0
  chi_array(29) = 20.d0
  chi_array(30) = 30.d0
  chi_array(31) = 50.d0
  chi_array(32) = 100.d0
  chi_array(33) = 200.d0
  chi_array(34) = 300.d0
  chi_array(35) = 500.d0
  chi_array(36) = 1000.d0
  chi_array(37) = 2000.d0
  chi_array(38) = 3000.d0
  chi_array(39) = 5000.d0
  chi_array(40) = 10000.d0

return
end subroutine set_chi_array


!****************************************************************************
!****************************************************************************
subroutine set_tzero_array(tzero_array)

! Initializes the values of the tzero_array, which is one of the arrays we
!   will use for pre-computing Eq (A19)
!
! Output argument:
!  1) tzero_array: array of tzero values at which we will calculate emission

use parameters
use constants

implicit none

  ! Output argument
real(kind=8), dimension(0:n_ph), intent(out) :: tzero_array

  ! Local variables
integer :: j_ph
real(kind=8) :: R_decel, t_decel, Gam_end, R_end, t_end, tzero_start,     &!&
     tzero_end, t_fac

  ! Compute engine-frame times marking start and end of tzero_array.  First
  !   entry will be 0.1 times t_decel (the time at which the Blandford-McKee
  !   solution has Gamma = Gam_cost).  Final entry will be 2 times tobs_end
  !   after converting from observer frame to engine frame.
  ! tzero_start precomputation
  if( k_cbm .eq. 0 ) then
    R_decel = E_iso*(17.d0-4.d0*k_cbm)                                    &!&
             / (8.d0*pii * n_cbm*rm_prot * Gam_coast**2)
  else
    R_decel = E_iso*(17.d0-4.d0*k_cbm)                                    &!&
             / (8.d0*pii * A_star*rm_prot * Gam_coast**2)
  endif
  R_decel = R_decel**(1.d0/(3.d0-k_cbm))
  t_decel = R_decel/ccgs
  
  ! tzero_end precomputation
  Gam_end = (17.d0 - 4.d0*k_cbm) * E_iso * (ccgs*tobs_end)**(k_cbm-3.d0)  &!&
           / ( 4.d0**(5.d0-k_cbm) * (4.d0-k_cbm)**(3.d0-k_cbm) * pii )
  if( k_cbm .eq. 0 ) then
    Gam_end = Gam_end / (n_cbm * rm_prot)
  else
    Gam_end = Gam_end / (A_star * rm_prot)
  endif
  Gam_end = sqrt(2.d0) * Gam_end**(0.5d0/(4.d0-k_cbm))
  if( k_cbm .eq. 0 ) then
    R_end = E_iso*(17.d0-4.d0*k_cbm)                                      &!&
             / (8.d0*pii * n_cbm*rm_prot * Gam_end**2)
  else
    R_end = E_iso*(17.d0-4.d0*k_cbm)                                      &!&
             / (8.d0*pii * A_star*rm_prot * Gam_end**2)
  endif
  R_end = R_end**(1.d0/(3.d0-k_cbm))
  t_end = R_end/ccgs
  
  ! Now actually compute tzero_start & tzero_end
  tzero_start = 0.1d0 * t_decel
  tzero_end   = 2.d0 * t_end
  
  
  ! Set tzero_array values
  t_fac = (tzero_end/tzero_start)**(1.d0/float(n_ph))
  tzero_array(0) = tzero_start
  do j_ph = 1, n_ph
    tzero_array(j_ph) = tzero_array(j_ph-1) * t_fac
  enddo

return
end subroutine set_tzero_array


!****************************************************************************
!****************************************************************************
subroutine set_phot_en_array(phot_en_array)

! Initializes the values of phot_en_array, which is one of the arrays we will
!   use for pre-computing Eq (A19), and for computing Eqs (A14) and (A24)
!
! Output argument:
!  1) phot_en_array: array of photon energies, in cgs units

use parameters
use constants

implicit none

  ! Output argument
real(kind=8), dimension(0:n_ph), intent(out) :: phot_en_array

  ! Local variables
integer :: j_ph
real(kind=8) :: en_fac

  ! Bottom of phot_en_array will have a value of 10^-7.5 eV, to include
  !   100 MHz = 4e-7 eV.  We'll use 8 steps per decade of photon energy
  phot_en_array(0) = 10.d0**(-8.5d0) * evterg
  en_fac = 10.d0**0.025d0
  do j_ph = 1, n_ph
    phot_en_array(j_ph) = phot_en_array(j_ph-1) * en_fac
  enddo

return
end subroutine set_phot_en_array


!****************************************************************************
!****************************************************************************
subroutine calc_syn_Fx(syn_F, syn_x)

! Pre-calculates F(x) in Rybicki & Lightman's Eq. 6.31c.  The arrays syn_F
!   and syn_x will be used in computing synchrotron emission.
!
! Input arguments:
!   None
! Output arguments:
!  1) syn_F: array containing F(x) values from Eq. 6.31c
!  2) syn_x: array containing x values associated with syn_F

use parameters, only: n_ph

implicit none

  ! Output arguments
real(kind=8), dimension(0:n_ph), intent(out) :: syn_F, syn_x

  ! Local variables
integer :: j, i, n_xxx
real(kind=8) :: xnu, xxx_max, xxx, xxx_log, del_xxx, xxx_fac, sum_k, xx1, &!&
     xx, xx2, del, ri, rk, rip, rkp


  ! Set xnu for use in bessik
  xnu = 5.d0/3.d0
  
  ! Initialize syn_x and syn_F
  syn_x(:) = 0.d0
  syn_F(:) = 0.d0
  
  
!--! Fill syn_x by hand
   !-------------------------------------------------------------------------
  syn_x(0:118) = (/  0.d0,                                                &!&
     1.0d-6, 2.0d-6, 3.0d-6, 5.0d-6, 1.0d-5, 2.0d-5, 3.0d-5, 5.0d-5,      &!&
     1.0d-4, 2.0d-4, 3.0d-4, 5.0d-4, 1.0d-3, 2.0d-3, 3.0d-3, 5.0d-3,      &!&
     1.0d-2, 2.0d-2, 3.0d-2, 4.0d-2, 5.0d-2, 6.0d-2, 7.0d-2, 8.0d-2,      &!&
     9.0d-2, 1.0d-1, 1.1d-1, 1.2d-1, 1.3d-1, 1.4d-1, 1.5d-1, 1.6d-1,      &!&
     1.7d-1, 1.8d-1, 1.9d-1, 2.0d-1, 2.1d-1, 2.2d-1, 2.3d-1, 2.4d-1,      &!&
     2.5d-1, 2.6d-1, 2.7d-1, 2.8d-1, 2.9d-1, 3.0d-1, 3.1d-1, 3.2d-1,      &!&
     3.3d-1, 3.4d-1, 3.5d-1, 3.6d-1, 3.7d-1, 3.8d-1, 3.9d-1, 4.0d-1,      &!&
     4.1d-1, 4.2d-1, 4.3d-1, 4.4d-1, 4.5d-1, 4.6d-1, 4.7d-1, 4.8d-1,      &!&
     4.9d-1, 5.0d-1, 6.0d-1, 7.0d-1, 8.0d-1, 9.0d-1, 1.0d+0, 1.1d+0,      &!&
     1.2d+0, 1.3d+0, 1.4d+0, 1.5d+0, 1.6d+0, 1.7d+0, 1.8d+0, 1.9d+0,      &!&
     2.0d+0, 2.1d+0, 2.2d+0, 2.3d+0, 2.4d+0, 2.5d+0, 2.6d+0, 2.7d+0,      &!&
     2.8d+0, 2.9d+0, 3.0d+0, 4.0d+0, 5.0d+0, 6.0d+0, 7.0d+0, 8.0d+0,      &!&
     9.0d+0, 1.0d+1, 1.1d+1, 1.2d+1, 1.3d+1, 1.4d+1, 1.5d+1, 1.6d+1,      &!&
     1.7d+1, 1.8d+1, 1.9d+1, 2.0d+1, 2.1d+1, 2.2d+1, 2.3d+1, 2.4d+1,      &!&
     2.5d+1, 2.6d+1, 2.7d+1, 2.8d+1, 2.9d+1, 3.0d+1 /)
  
  
!--! Now calculate F(x) at each value of syn_x
   !-------------------------------------------------------------------------
  do j = 1, 117
    
    xxx_max = 30.d0
    xxx     = syn_x(j)
    xxx_log = log10(xxx)
    
    n_xxx   = 200
    del_xxx = (log10(xxx_max) - xxx_log)/real(n_xxx)
    xxx_fac = 10.d0**del_xxx
    
    sum_k = 0.0d00
    xx1   = xxx
    xx    = xx1 / sqrt(xxx_fac)  ! Both of these are initialized so they
    xx2   = xxx                  !   have the right values in the loop
    do i = 1, n_xxx
      xx2 = xx2 * xxx_fac  !  = 10.0**(xxx_log + (i)*del_xxx)
      xx  = xx  * xxx_fac  !  = sqrt(xx1*xx2)
      del = xx2 - xx1
      
!     ! Modified Bessel function of fractional order
      call bessik(xx, xnu, ri, rk, rip, rkp, 1.d-8)
      
!     ! Here is the factor accounting for isotropy:
      sum_k = sum_k + sqrt(1.0d00-(xxx/xx)**2)*rk*del
!     ! Use this if not accounting for isotropy:
!      sum_k = sum_k + rk*del
      xx1   = xx2
    enddo ! loop over n_xxx
    
    syn_F(j) = xxx*sum_k
  
  enddo
  
  syn_F(118) = 0.d0  ! Set this one by hand at the end
   !-------------------------------------------------------------------------
   ! F(x) calculated

return
end subroutine calc_syn_Fx


!****************************************************************************
!****************************************************************************
subroutine bessik(x,xnu,ri,rk,rip,rkp,EPS)

! Returns the modified Bessel functions ri = I_nu, rk = K_nu, and their
! derivatives rip = I'_nu, rkp = K'_nu, for positive x and xnu (= nu) .ge. 0.
!
! The relative accuracy is within one or two significant digits of the input
! argument EPS. FPMIN is a parameter set close to the machine's smallest
! floating point number. All internal arithmetic is in double precision.
!
! Uses beschb, which uses chebev.
!
! Subroutine taken almost unmodified from Numerical Recipes in Fortran 77,
!  2nd edition.

implicit none

  ! Input arguments
real(kind=8) :: x, xnu, EPS
  ! Parameters
integer :: MAXIT
real(kind=8) :: XMIN, FPMIN, PI
parameter (FPMIN=1.0d-30, MAXIT=10000, XMIN=2.0, PI=3.141592653589793d0)
  ! Output arguments
real(kind=8) :: ri, rip, rk, rkp
  ! Local variables
integer :: i,l,nl
real(kind=8) :: a, a1, b, c, d, del, del1, delh, dels, e, f, fact, fact2, &!&
     ff, gam1, gam2, gammi, gampl, h, p, pimu, q, q1, q2, qnew, ril, ril1,&!&
     rimu, rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2,&!&
     xi, xi2, xmu, xmu2

if(x .le. 0.0 .or. xnu .lt. 0.0) then
  write(*,*) 'bad arguments in bessik', x, xnu
  read(*,*)
endif 

! n is the number of downward recurrences of the I's and upward recurrences
!  of the K's. xmu lies between -1/2 and +1/2
nl  = int(xnu + 0.5d0)
xmu  = xnu-nl
xmu2 = xmu*xmu
xi  = 1.0d0/x
xi2 = 2.0d0*xi
h   = xnu*xi

! Evaluate CF1 by modified Lentz's method (section 5.2 in Numerical Recipes
!  in Fortran 77, 2nd ed.)
if(h .lt. FPMIN) h = FPMIN
b = xi2*xnu
d = 0.0d0
c = h

do i = 1, MAXIT
  b   = b + xi2
  d   = 1.0d0/(b+d)
  c   = b + 1.d0/c
  del = c*d
  h   = del*h
  if(abs(del-1.0d0) .lt. EPS) exit
enddo
if(i .gt. MAXIT) then
  write(*,*) ' x too large in bessik; try asymptotic expansion'
  read(*,*)
endif

! Initialize I_nu and I'_nu for downward recurrence, and store values for
!  later rescaling
ril  = FPMIN
ripl = h*ril
ril1 = ril
rip1 = ripl
fact = xnu*xi

do l = nl,1,-1
  ritemp = fact*ril+ripl
  fact   = fact-xi
  ripl   = fact*ritemp+ril
  ril    = ritemp
enddo

! We now have unnormalized I_mu and I'_mu.  Use the series
f = ripl/ril
if(x .lt. XMIN) then
  x2   = 0.5d00*x
  pimu = PI*xmu

  if(abs(pimu) .lt. EPS) then
    fact = 1.0d00
  else
    fact = pimu/sin(pimu)
  endif

  d = -log(x2)
  e = xmu*d

  if(abs(e) .lt. EPS) then
    fact2 = 1.0d00
  else
    fact2 = sinh(e)/e
  endif

  ! Chebyshev evaluation of Gamma_1 and Gamma_2
  call beschb(xmu, gam1, gam2, gampl, gammi)

  ff  = fact*(gam1*cosh(e) + gam2*fact2*d)  ! f0
  sum = ff
  e   = exp(e)
  p = 0.5d0*e/gampl    ! p0
  q = 0.5d0/(e*gammi)  ! q0
  c = 1.0d00
  d = x2*x2
  sum1 = p

  do i = 1, MAXIT
    ff = (i*ff+p+q)/(i*i-xmu2)
    c  = c*d/i
    p  = p/(i-xmu)
    q  = q/(i+xmu)
    del = c*ff
    sum = sum + del
    del1 = c*(p-i*ff)
    sum1 = sum1 + del1
    if(abs(del) .lt. abs(sum)*EPS) exit
  enddo
  if(i .gt. MAXIT) then
    write(*,*) 'bessk series failed to converge'
    read(*,*) 
  endif

  rkmu = sum
  rk1  = sum1*xi2

! Otherwise, evaluate CF2 by Steed's algorithm (also section 5.2 of Numerical
!  Recipes in Fortran 77, 2nd ed.).  This is OK because there can be no zero
!  denominators
else

  b = 2.0d00*(1.0d0 + x)
  d = 1.0d00/b
  delh = d
  h    = delh
  q1 = 0.0d00   ! Initializations for recurrence 6.7.35
  q2 = 1.0d00
  a1 = 0.25d00 - xmu2
  c = a1
  q = c         ! First term in equation 6.7.34
  a = -a1
  s = 1.0d00 + q*delh

  do i = 2, MAXIT
    a = a - 2*(i-1)
    c = -a*c/i
    qnew = (q1 - b*q2)/a
    q1 = q2
    q2 = qnew
    q  = q + c*qnew
    b  = b + 2.0d0
    d  = 1.0d0/(b + a*d)
    delh = (b*d - 1.0d00)*delh
    h    = h + delh
    dels = q*delh
    s    = s + dels
    if(abs(dels/s) .lt. EPS) exit ! Only need to test convergence of sum
                                  !  since CF2 itself converges more quickly
  enddo
  if(i .gt. MAXIT) then
    write(6,"(A)") ' bessik: failure to converge in cf2'
    stop
  endif

  h = a1*h
  if(x .gt. 350.0d00) then
    rk = 1.0d-150
    return
  endif
  
  rkmu = sqrt(PI/(2.d0*x))*exp(-x)/s  ! Can omit the factor exp(-x) to scale
  rk1  = rkmu*(xmu+x+.5d0-h)*xi       !  all the returned functions by exp(x)
                                      !  for all x .ge. XMIN
endif

rkmup = xmu*xi*rkmu - rk1
rimu  = xi/(f*rkmu - rkmup)  ! Get I_mu from the Wronskian
ri    = (rimu*ril1)/ril      ! Scale original I_nu and I'_nu
rip   = (rimu*rip1)/ril

! Now use the upward recurrence of K_nu
do i = 1, nl
  rktemp = (xmu + i)*xi2*rk1 + rkmu
  rkmu   = rk1
  rk1    = rktemp
enddo

rk  = rkmu
rkp = xnu*xi*rkmu - rk1

return
end subroutine bessik


!****************************************************************************
!****************************************************************************
subroutine beschb(x,gam1,gam2,gampl,gammi)

! Evaluates Gamma_1 and Gamma_2 by Chebyshev expansion for abs(x) .le. 1/2.
! Also returns 1/Gamma(1+x) and 1/Gamma(1-x).
!
! Uses chebev_subr
!
! Subroutine taken almost unmodified from Numerical Recipes in Fortran 77,
!  2nd edition.

implicit none

  ! Input arguments
real(kind=8) :: x
  ! Parameters
integer :: NUSE1, NUSE2
parameter ( NUSE1 = 7, NUSE2 = 8 )
  ! Output arguments
real(kind=8) :: gam1, gam2, gammi, gampl
  ! Local variables
real(kind=8) :: xx, c1(7), c2(8), chebev_val

c1(1) = -1.142022680371172d00
c1(2) =  6.516511267076d-03
c1(3) =  3.08709017308d-04
c1(4) = -3.470626964d-06
c1(5) =  6.943764d-09
c1(6) =  3.6780d-11
c1(7) = -1.36d-13

c2(1) =  1.843740587300906d00
c2(2) = -0.076852840844786d00
c2(3) =  1.271927136655d-03
c2(4) = -4.971736704d-06
c2(5) = -3.3126120d-08
c2(6) =  2.42310d-10
c2(7) = -1.70d-13
c2(8) = -1.0d-15

xx = 8.0d0*x*x - 1.0d0 ! Multiply x by 2 to convert range to -1 to 1, and
                       !  then apply the transformation for evaluating even
                       !  Chebyshev series
call chebev_subr(-1.0d00, 1.0d00, c1, NUSE1, xx, chebev_val)
gam1 = chebev_val
call chebev_subr(-1.0d00, 1.0d00, c2, NUSE2, xx, chebev_val)
gam2  = chebev_val

gampl = gam2 - x*gam1
gammi = gam2 + x*gam1

return
end subroutine beschb


!****************************************************************************
!****************************************************************************
subroutine chebev_subr(a, b, c, m, x, chebev_val)

! Chebyshev evaluation. All arguments except the last are input. c(1:m) is
!  an array of Chebyshev coefficients, the first m elements of c, which is
!  the output of another subroutine (chebft, which must have previously been
!  called with the same value of a and b).
!
! The Chebyshev polynomial is evaluated at a point y defined below, and the
!  result is returned in chebev_val.
!
! Subroutine taken from Numerical Recipes in Fortran 77, 2nd edition.

  ! Input arguments
integer :: m
real(kind=8) :: a, b, x, c(m)
  ! Output arguments
real(kind=8) :: chebev_val
  ! Local variables
integer :: j
real(kind=8) :: d, dd, sv, y, y2

if ((x-a)*(x-b) .gt. 0.0d00) then
  write(*,*) 'x not in range in chebev'; read(*,*)
endif 

d  = 0.0d00
dd = 0.0d00
y  = (2.0d00*x - a - b)/(b -a)
y2 = 2.0d00*y

do j = m, 2, -1
  sv = d
  d  = y2*d - dd + c(j)
  dd = sv
enddo

chebev_val = y*d - dd + 0.5d00*c(1)

return
end subroutine chebev_subr


!****************************************************************************
!****************************************************************************
subroutine compute_Psyn(j_ph, i_tm, k_chi, Psyn)

! Uses Equation (A19) of Granot & Sari (2002ApJ...568..820G) to compute the
!   synchrotron emissivity.
!
! Input arguments:
!  1) j_y: index of current photon energy
!  2) i_tm: index of current tzero
!  3) k_chi: index of current chi
! Output arguments:
!  1) Psyn: integrated emissivity

use parameters
use constants
use data_arrays

implicit none

interface
  double precision function Psyn_integrand(log_gam, r_params, i_params)
    real(kind=8), intent(in) :: log_gam
    integer, intent(in) :: i_params(7)
    real(kind=8), intent(in) :: r_params(7)
  end function
end interface

  ! Input arguments
integer, intent(in) :: i_tm, j_ph, k_chi
  ! Output argument
real(kind=8), intent(out) :: Psyn

  ! Local variables
integer :: istat
real(kind=8) :: phot_en, t_zero, chi, R_zero, n_ext, Gam_zero, n_zero,    &!&
     e_zero, Bsq_zero, gam_min_uncool, gam_max_uncool, K_zero,            &!&
     gam_max_cool, gam_min_cool, B_sq, u_bound, l_bound
integer, dimension(7) :: i_params
real(kind=8), dimension(7) :: r_params

   !-------------------------------------------------------------------------
!--! (0) Get values from data arrays
   !-------------------------------------------------------------------------
  phot_en = phot_en_array(j_ph)
  t_zero  = tzero_array(i_tm)
  chi     = chi_array(k_chi)
  
   !-------------------------------------------------------------------------
!--! (1) Compute conditions at time fluid element crossed shock
   !-------------------------------------------------------------------------
  R_zero   = ccgs * t_zero
  if( k_cbm .eq. 0 ) then
    n_ext = n_cbm
  else
    n_ext = A_star * R_zero**(-k_cbm)
  endif
  Gam_zero = (17.d0 - 4.d0*k_cbm) * E_iso * (ccgs*t_zero)**(k_cbm-3.d0)   &!&
            / ( 8.d0*pii * n_ext*rm_prot )
  Gam_zero = sqrt(Gam_zero)
  
  n_zero   = sqrt(8.d0) * Gam_zero * n_ext
  e_zero   = 2.d0 * Gam_zero**2 * n_ext*rm_prot
  Bsq_zero = 8.d0*pii * epsB0 * e_zero
   !-------------------------------------------------------------------------
   ! Zero-time conditions found
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
!--! (2) Compute gam_min_cool and gam_max_cool
   !-------------------------------------------------------------------------
  ! Compute minimum and maximum electron Lorentz factors for the integral,
  !   using Equations (A1) for uncooled gamma_min, (A12) for uncooled
  !   gamma_max, and (A11) for the cooled versions of both.  Also find
  !   K_zero, which sets the normalization.
  ! Uncooled gamma_max should be infinity, according to Granot & Sari (2002),
  !   but we use here a Lorentz factor of 6e14 as a finite but absurdly large
  !   upper bound on the electron distribution, corresponding to an electron
  !   with an energy of 3x10^20 eV.
  gam_min_uncool = epse_bar * e_zero / (n_zero * rm_elec)
  
  gam_max_uncool = 6.d14
  
  K_zero = (p_spec - 1.d0) * n_zero * gam_min_uncool**(p_spec-1.d0)
  
  
  ! Cool the electrons
  if( chi .gt. 1.d0 ) then
    gam_max_cool = sqrt(2.d0)*(19.d0-2.d0*k_cbm) * pii * xme*ccgs         &!&
                  * Gam_zero * chi**( (25.d0-2.d0*k_cbm)                  &!&
                                     /(24.d0-6.d0*k_cbm) )                &!&
                  / ( sigT * 8.d0*pii * epsB0 * e_zero * t_zero           &!&
                     * ( chi**( (19.d0-2.d0*k_cbm)                        &!&
                               /(12.d0-3.d0*k_cbm) )  -  1.d0 ) )
  else
    gam_max_cool = gam_max_uncool
  endif
  gam_max_cool = min( gam_max_cool, gam_max_uncool )

  if( chi .gt. 1.d0 ) then
    gam_min_cool = gam_min_uncool / ( chi**( (13.d0-2.d0*k_cbm)           &!&
                                            /(24.d0-6.d0*k_cbm) )         &!&
                                     + gam_min_uncool/gam_max_cool )
  else
    gam_min_cool = gam_min_uncool
  endif
  gam_min_cool = min( gam_min_cool, gam_min_uncool )
   !-------------------------------------------------------------------------
   ! gam_min_cool and gam_max_cool found
   !-------------------------------------------------------------------------
  
  
   !-------------------------------------------------------------------------
!--! (3) Calculate synchrotron production in the plasma frame
   !-------------------------------------------------------------------------
  ! Local magnetic field value
  B_sq = 8.d0*pii * epsB0 * e_zero * chi**(- (26.d0-4.d0*k_cbm)           &!&
                                            /(12.d0-3.d0*k_cbm) )
  
  ! Integrate the emissivity if electron Lorentz factors are high enough
  !   (gamma > 6).  Otherwise return 0
  if( gam_max_cool .lt. 6.d0 ) then
    Psyn = 1.d-99
  else
    ! r_params(1): product of prefactors in Eqs (A13) and (A17)
    ! r_params(2): photon energy, needed to compute syn_x
    ! r_params(3): numerical factor in nu_syn
    ! r_params(4): maximum electron Lorentz factor for Eq (A13)
    r_params(1) = K_zero * chi**( (2.d0*k_cbm-13.d0)*(p_spec+2.d0)        &!&
                                 /(24.d0-6.d0*k_cbm) )                    &!&
                 * sqrt(3.d0*B_sq) * qcgs**3.d0 / rm_elec
    r_params(2) = phot_en
    r_params(3) = 3.d0 * qcgs * sqrt(B_sq) / (4.d0*pii * xme*ccgs)
    r_params(4) = gam_max_cool
    
    u_bound = log(gam_max_cool)
    l_bound = log( max( 6.d0, gam_min_cool ) )
    
    call romb_quad( Psyn_integrand, l_bound, u_bound, r_params, i_params, &!&
       Psyn, istat )
  endif
   !-------------------------------------------------------------------------
   ! Synchrotron emissivity calculated
   !-------------------------------------------------------------------------

return
end subroutine compute_Psyn


!****************************************************************************
!****************************************************************************
double precision function Psyn_integrand(log_gam, r_params, i_params)

! Evaluates the integrand of Eq (A19), but in log space (i.e. with an extra
!   factor of gam_e to make up for the dgam_e --> dln(gam_e) conversion)
!
! Input arguments:
!   1) log_gam: current value being tried for ln(gam_e)
!   2) i_params: array (length 7) of up to 7 integer parameters to be used
!     in evaluating the function
!   3) r_params: array (length 7) of up to 7 floating point parameters to
!     be used in evaluating the function
! Output argument:
!   1) Psyn_integrand: the value of the integrand in Eq (A19)

use parameters
use constants
use data_arrays

implicit none

  ! Input arguments
real(kind=8), intent(in) :: log_gam
integer, dimension(7), intent(in) :: i_params
real(kind=8), dimension(7), intent(in) :: r_params
  ! Output argument is the function name; no need to declare it, apparently

  ! Local variables
integer :: i_x, i_x_lo
real(kind=8) :: gam_e, Psyn_prefac, phot_en, nu_c_prefac, gam_max,        &!&
     dist_term, xxx, F_x

  ! We might be integrating in log space, but the function uses gamma_e in
  !   linear space for all its evaluation.  So return it to linear space here
  gam_e = exp(log_gam)
  
  
  ! Read in the data from r_params
  Psyn_prefac = r_params(1)
  phot_en     = r_params(2)
  nu_c_prefac = r_params(3)
  gam_max     = r_params(4)
  
  
  ! Electron distribution is Equation (A13) in Granot & Sari (2002)
  ! Extra factor of gam_in is because dgam = gam*dgam/gam = gam * d(ln[gam]),
  !   and we are integrating in log space.
  if( gam_e .lt. gam_max ) then
    dist_term = gam_e**(-p_spec+1.d0)                                    &!&
               * ( 1.d0 - gam_e/gam_max )**(p_spec - 2.d0)
  else
    dist_term = 0.d0
  endif
  
  
  ! The "isotropic" form of syn_F in calc_syn_Fx is within 5% of the correct
  !   value everywhere in the domain of computation, and is typically closer
  !   than 1%.  Use that instead of repeated calls to integrate Equation
  !   (A18).
  ! The "isotropic" form of syn_F can be calculated directly.  Note the
  !   factor of the Planck constant in xxx, which converts nu_c frequency to
  !   energy to match phot_en
  !--------------------------------------------------------------------------
  xxx = phot_en / (nu_c_prefac * xh * gam_e**2)
  
  ! At large x, no synchrotron power produced
  if( xxx .ge. 30.d0 ) then
    Psyn_integrand = 1.d-99
    return
  endif
  
  ! Find the appropriate value of syn_x for this photon energy
  do i_x = 0, n_ph-1
    if( syn_x(i_x+1) .ge. xxx ) then
      i_x_lo = i_x
      exit
    endif
  enddo
  
  ! Interpolate between values of syn_F to get F(x).  If x < syn_x(1), then
  !   it falls in the power-law part of the function F(x)
  if( xxx .lt. syn_x(1) ) then
    F_x = (xxx*o_o_sone)**third * syn_F(1)
  else
    F_x = syn_F(i_x_lo)                                                   &!&
         +  (xxx - syn_x(i_x_lo)) / (syn_x(i_x_lo+1) - syn_x(i_x_lo))     &!&
           * (syn_F(i_x_lo+1) - syn_F(i_x_lo))
  endif
  !--------------------------------------------------------------------------
  ! Isotropic synchrotron power found
  
  
  ! Integrand (and so result of this function) is the product of the number
  !   of electrons and power per electron
  Psyn_integrand = Psyn_prefac * dist_term * F_x

return
end function Psyn_integrand


!****************************************************************************
!****************************************************************************
subroutine romb_quad(func, x_start, x_end, r_params, i_params, integral,  &!&
     istat)

! Integrates a function from provided starting point to ending point.  Uses
!   Romberg's method of order 2*k (where, e.g., k = 2 is Simpson's rule).
!
! Input arguments
!   1) func: name of (external) function computing the value of the
!     integrand.  Must be double precision function
!   2) x_start: initial value of x
!   3) x_end: final value of x
!   4) r_params: array (length 7) of up to 7 floating point parameters to
!     be used by func
!   5) i_params: array (length 7) of up to 7 integer parameters to be used
!     by func
! Output arguments
!   1) integral: value of the integral
!   2) istat: status variable.  1 if integral converged, -1 otherwise
!
! Subroutine taken from qromb in Numerical Recipes in Fortran 77, 2nd ed.
!   It is explained in Section 4.3.

implicit none

interface
  double precision function func(x_in, r_params, i_params)
    real(kind=8), intent(in) :: x_in
    integer, intent(in) :: i_params(7)
    real(kind=8), intent(in) :: r_params(7)
  end function
end interface

  ! Input arguments
real(kind=8), intent(in) :: x_start, x_end
integer, dimension(7), intent(in) :: i_params
real(kind=8), dimension(7), intent(in) :: r_params
  ! Output arguments
integer, intent(out) :: istat
real(kind=8), intent(out) :: integral

  ! Local variables
integer, parameter :: jmax = 20         ! Maximum number of refining
integer, parameter :: jmaxp = jmax + 1  !   iterations to perform
integer, parameter :: k_ord = 5         ! Order of Romberg's method; number
integer, parameter :: km = k_ord - 1    !   of points used in extrapolation
real(kind=8), parameter :: tol = 1.d-5  ! Accuracy desired
integer :: j
real(kind=8) :: s_in, d_int
real(kind=8), dimension(jmaxp) :: h, s
integer, dimension(7) :: i_params_int
real(kind=8), dimension(7) :: r_params_int
  
  h(1) = 1.d0
  s_in = 0.d0 ! Irrelevant since it will be ignored when j = 1 and replaced
              !   immediately thereafter
  
  ! Loop over refinements to the integral, halving the step size at
  !   each cycle
  do j = 1, jmax
    
    i_params_int(:) = i_params(:)
    i_params_int(5) = j  ! 1 = i_tm, 2 = k_ang, 3 = n_step, 4 = m_ph,
                         !   and 5 = j
    r_params_int(:) = r_params(:)
    call trapzd(func, x_start, x_end, r_params_int, i_params_int, s_in, j,&!&
       s(j))
    
    if( j .ge. k_ord ) then
      ! Warning: polint expects arrays for h and s, so what we are
      !   actually passing are sub-arrays of h and s running from element
      !   j-km to j.
      call polint(h(j-km), s(j-km), k_ord, 0.d0, integral, d_int)
      if( abs(d_int) .le. (tol * abs(integral)) ) exit
    endif
    
    s_in   = s(j)
    h(j+1) = 0.25d0 * h(j)  ! This is a key step.  The factor is 0.25 even
                            !   though we only decrease our step size by 0.5.
                            !   This makes the extrapolation a polynomial in
                            !   h^2 rather than h, as allowed by Eq. 4.2.1,
                            !   the Euler-Maclaurin Summation formula.
  enddo
  
  if( j .gt. jmax ) then
    istat = -1
  else
    istat = 1
  endif
  
return
end subroutine romb_quad


!****************************************************************************
!****************************************************************************
subroutine trapzd(func, x_start, x_end, r_params, i_params, s_in, n, s_out)

! Computes the nth stage of refinement of an extended trapezoidal rule.
!   When called with n = 1, the routine returns the crudest estimate of the
!   integral of func from x_start to x_end.  Subsequent calls with n = 2+
!   will improve the integral by adding 2^(n-2) additional interior points.
! Do not modify s betweeen sequential calls to trapzd!
!
! Input arguments
!   1) func: name of (external) function computing the value of the
!     integrand.  Must be double precision function
!   2) x_start: initial value of x
!   3) x_end: final value of x
!   4) r_params: array (length 7) of up to 7 floating point parameters to
!     be used by func
!   5) i_params: array (length 7) of up to 7 integer parameters to be used
!     by func
!   6) s_in: value of integral at previous level of refinement (irrelevant
!     if n = 1)
!   7) n: level of refinement
! Output argument
!   1) s_out: value of the integral after refinement
!
! Subroutine taken from trapzd in Numerical Recipes in Fortran 77, 2nd ed.
! It is explained in Section 4.2.

implicit none

interface
  double precision function func(x_in, r_params, i_params)
    real(kind=8), intent(in) :: x_in
    integer, intent(in) :: i_params(7)
    real(kind=8), intent(in) :: r_params(7)
  end function
end interface

  ! Input arguments
integer, intent(in) :: n
real(kind=8), intent(in) :: x_start, x_end, s_in
integer, dimension(7), intent(in) :: i_params
real(kind=8), dimension(7), intent(in) :: r_params
  ! Output argument
real(kind=8), intent(out) :: s_out

  ! Local variables
integer :: it, j
real(kind=8) :: del, x, running_tot
  
  if( n .eq. 1 ) then
    
    s_out = 0.5d0 * (x_end-x_start) * ( func(x_end, r_params, i_params)   &!&
                                       + func(x_start, r_params, i_params) )
    
  else
    
    it  = 2**(n - 2)
    del = (x_end - x_start) / (1.d0*it)
    x   = x_start + 0.5d0*del
    running_tot = 0.d0
    
    do j = 1, it
      running_tot = running_tot  +  func(x, r_params, i_params)
      x = x + del
    enddo
    
    s_out = 0.5d0 * ( s_in  +  del * running_tot )
    
  endif
  
return
end subroutine trapzd


!****************************************************************************
!****************************************************************************
subroutine polint(x_array, y_array, n, x_in, y_out, dy_out)

! Uses Neville's algorithm to interpolate (or extrapolate) a polynomial
!   passing through the points provided.
!
! Input arguments
!   1) x_array: array containing x values of points
!   2) y_array: array containing y values of points
!   3) n: number of points provided, and one more than degree of
!     interpolating polynomial
!   4) x_in: value at which we want the polynomial evaluated
! Output arguments
!   1) y_out: value of polynomial at x_in
!   2) dy_out: error estimate associated with y_out
!
! Subroutine taken from polint in Numerical Recipes in Fortran 77, 2nd ed.
! It is explained in Section 3.1.

implicit none

  ! Input arguments
integer, intent(in) :: n
real(kind=8), intent(in) :: x_in
real(kind=8), dimension(n), intent(in) :: x_array, y_array
  ! Output arguments
real(kind=8), intent(out) :: y_out, dy_out

  ! Local variables
integer :: ns, i, m
real(kind=8) :: dif, dift, ho, hp, w, den
real(kind=8), dimension(n) :: c, d
  
  ns = 1
  dif = abs( x_in - x_array(1) )
  
  ! Find the index ns of the closest table entry and initialize the tableau
  !   of c's and d's
  do i = 1, n
    dift = abs( x_in - x_array(i) )
    if( dift .lt. dif ) then
      ns  = i
      dif = dift
    endif
    
    c(i) = y_array(i)
    d(i) = y_array(i)
  enddo
  
  
  ! The initial approximation to y_out
  y_out = y_array(ns)
  
  
  ! For each column of the tableau, loop over the current c's and d's and
  !   update them
  ns = ns - 1
  do m = 1, n-1
    
    do i = 1, n-m
      ho = x_array(i)   - x_in
      hp = x_array(i+m) - x_in
      w  = c(i+1) - d(i)
      
      den = ho - hp
      ! This error will only occur if two input entries of x_array are
      !   identical to within roundoff error
      if( den .eq. 0.d0 ) then
        write(*,"(A)") 'ERROR in subroutine polint: den = 0'
        write(*,"(A)") 'Stopping program now'
        stop
      endif
      
      den = w / den
      
      ! Update the c's and d's
      d(i) = hp*den
      c(i) = ho*den
    enddo
    
    ! After each column in the tableau is computed, decide which correction
    !   (c or d) we want to add to our accumulating value of y_out.  That is,
    !   decide whith path to take through the tableau--forking up or down.
    ! Do this in such a way as to take the most "straight line" route through
    !   the tableau to its apex, updating ns accordingly to keep track of
    !   where we are.  This route keeps the partial approximations centered
    !   (so far as possible) on the target location x_in.
    ! The last dy_out added is thus the error approximation.
    if( (2*ns) .lt. (n-m) ) then
      dy_out = c(ns + 1)
    else
      dy_out = d(ns)
      ns = ns - 1
    endif
    
    y_out = y_out + dy_out
    
  enddo
  
return
end subroutine polint


!****************************************************************************
!****************************************************************************
subroutine cosmo_calc(d_CM, z)

! Calculator to shift between redshift, lookback time, and comoving distance.
! Subroutine adapted from Hogg (1999). For more detail see
! [http://adsabs.harvard.edu/abs/1999astro.ph..5116H]
!
! Note one major difference: Equation (13) in Hogg uses Omega_r where this
! program uses Omega_k, and does not include a term for what this program
! calls Omega_r. For more information, see Wright (2006):
! [http://adsabs.harvard.edu/abs/2006PASP..118.1711W]
!
! The program calculates lookback time (in years) as well, but since this
!  result is not necessary for GRB redshift calculation it is not passed
!  back to the calling subroutine.
!
! Input values:
!  1) One of (redshift, comoving distance in Mpc)
! Output values:
!  1) Other of (redshift, comoving distance in Mpc)

implicit none

  ! Input/output arguments
real(kind=8), intent(inout) :: z, d_CM

  ! Parameters; cosmology values pulled from Planck 2013 results
                     ! Hubble constant, km/s Mpc^-1
real(kind=8), parameter :: H0 = 67.7d0
                     ! Fraction of density in neutrinos?
real(kind=8), parameter :: Omega_r = 0.4165d0 / H0**2
                     ! Fraction of density in dark energy and in matter
real(kind=8), parameter :: Omega_vac = 0.689d0 - 0.5d0*Omega_r
real(kind=8), parameter :: Omega_m   = 0.311d0 - 0.5d0*Omega_r
                     ! Assume flat Universe, i.e. Omega_k = 1 - sum(Omega_*)
real(kind=8), parameter :: Omega_k = 0.d0
                     ! Speed of light, km/s
real(kind=8), parameter :: c = 2.99792458d5
                     ! Hubble distance, Mpc
real(kind=8), parameter :: d_H = c / H0
                     ! Hubble time, years
real(kind=8), parameter :: t_H = 9.778d11 / H0
                     ! Number of steps to use in integration
integer, parameter ::  num_steps = 1000

  ! Local variables
integer :: code_stat, i
real(kind=8) :: t_look, z_save, d_save, d_old, d_new, z_old, E_old, z_new,&!&
     E_new, slope

1001  format(A, ES10.4E2, A)
2001  format(2(A, ES10.4E2))


!--! Set either redshift or comoving distance to a nonzero value. Other
   !   should be set to zero
  z_save = z
  d_save = d_CM
  code_stat = -1 ! Status of code
                 !  -1: initialization value
                 !   1: value of z for provided d_CM located. Proceed to
                 !       integration for lookback time
                 !   2: d_CM less than critical value. Skip stepping through
                 !       z integral and integration to find t_look
  
  
!--! Quick error check that at least one input is physically reasonable and
   !   that only one is positive
  if( ((z .le. 0.d0) .and. (d_CM .le. 0.d0)) .or.                         &!&
      ((z .gt. 0.d0) .and. (d_CM .gt. 0.d0)) ) then
    write(*,"(2A)") 'Invalid inputs given.  Exactly one of z and d_CM ',  &!&
                    'can/must be positive.'
    write(*,2001)   'Inputted values were z = ', z_save, ' and d_CM = ', d_CM
    write(*,"(A)") 'Stoping program now.'
    stop
  endif


!--! If d_CM is given, step through integral in z (Equation (14) in Hogg
   !   1999) until next step exceeds d_CM
   !-------------------------------------------------------------------------
  if( (z .eq. 0.d0) .and. (d_CM .gt. 0.d0)                                &!&
    ! First, check to make sure d_CM is large enough to get a reasonable
    !  redshift -- otherwise return zero
    .and. (d_CM .lt. 0.443d0) ) then
    
    z      = 0.d0
    t_look = d_CM * 3.2616d6 ! convert Mpc to years
    code_stat = 2 ! d_CM less than critical value. Skip stepping through
                    !  z integral and integration to find t_look
  
  
  ! d_CM was a reasonable value
  else if( (z .eq. 0.d0) .and. (d_CM .gt. 0.d0) ) then
  
    d_old = 0.d0
    d_new = 0.d0
    z_old = 0.d0
    E_old = 1.d0 ! E(z) = sqrt[Omega_r*(1+z)^4 + Omega_m*(1+z)^3 
                 !           + Omega_k*(1+z)^2 + Omega_vac]
                 !  Equation (13) in Hogg (1999)
  
    
    do i = 1, 100
      z_new = 1.d-4 * real(i)
      E_new = sqrt( Omega_r*(1+z_new)**4                  &!& Equation (13)
                   + Omega_m*(1+z_new)**3                 &!&  in Hogg (1999)
                   + Omega_vac)
      
      ! Use Equation (14) in Hogg (1999) to update d_new
      d_new  = d_new  +  d_H*1.d-4 * 0.5d0*(1.d0/E_old + 1.d0/E_new)
      
      ! If current step in z exceeded provided value of d_CM, perform linear
      !  interpolation between the d_old and d_new
      if( (d_old .lt. d_CM) .and. (d_new .ge. d_CM) ) then
        
        slope = (d_new - d_old) / (z_new - z_old)
        z = (d_CM - d_old) / slope  +  z_old
        code_stat = 1 ! value of z for provided d_CM located. Proceed to
                      !  integration for lookback time
        exit
        
      endif
      
      ! Set variables for next pass through loop
      d_old = d_new
      E_old = E_new
      z_old = z_new
    enddo
    
    
    ! If a value of z wasn't already found stepping from z = 0 to 0.01 with
    !   steps of size 1.d-4, start stepping up from z = 0.01 with steps of
    !   size 0.01
    i = 0
    do while (code_stat .lt. 1)
      i     = i + 1
      z_new = 1.d-2 * real(i)
      E_new = sqrt( Omega_r*(1+z_new)**4                  &!& Equation (13)
                   + Omega_m*(1+z_new)**3                 &!&  in Hogg (1999)
                   + Omega_vac)
      
      ! Use Equation (14) in Hogg (1999) to update d_new
      d_new  = d_new  +  d_H*1.d-2 * 0.5d0*(1.d0/E_old + 1.d0/E_new)
      
      ! If current step in z exceeded provided value of d_CM, perform linear
      !  interpolation between the d_old and d_new
      if( (d_old .lt. d_CM) .and. (d_new .ge. d_CM) ) then
        
        slope = (d_new - d_old) / (z_new - z_old)
        z = (d_CM - d_old) / slope  +  z_old
        code_stat = 1 ! value of z for provided d_CM located. Proceed to
                      !  integration for lookback time
        exit
        
      endif
      
      ! Set variables for next pass through loop
      d_old = d_new
      E_old = E_new
      z_old = z_new
      
      
      ! Quick and dirty error check
      if( i .gt. 1000) then
        write(*,"(A,F5.2)") "z_new exceeds reasonable bounds: ",z_new
        write(*,"(A)") "Stopping program"
        stop
      endif
      
    enddo
    
  endif
   !-------------------------------------------------------------------------
   ! End calculation of z for a provided d_CM


!--! If z is given, perform integral to find d_CM using the trapezoid rule
   !   and Equation (14) in Hogg (1999).
   ! Independently of the previous, calculate t_look using the trapezoid rule
   !   and Equations (28) in Hogg (1999).
   !-------------------------------------------------------------------------
  if( ((z .gt. 0.d0).and.(d_CM .eq. 0.d0)) .or. (code_stat .eq. 1) ) then
    
    t_look = 0.d0
    z_old = 0.d0
    E_old = 1.d0 ! E(z) = sqrt[Omega_m*(1+z)^3 + Omega_r*(1+z)^2 + Omega_vac]
                 !  Equation (13) in Hogg (1999)
  
    
    do i = 1, num_steps
      z_new = z*real(i)/num_steps
      E_new = sqrt( Omega_r*(1+z_new)**4                  &!& Equation (13)
                   + Omega_m*(1+z_new)**3                 &!&  in Hogg (1999)
                   + Omega_vac)
  
      
      ! If d_CM wasn't provided, use Equation (14) in Hogg (1999) to update
      !  running total
      if( code_stat .lt. 1 )  d_CM = d_CM  +  0.5d0*(1.d0/E_old + 1.d0/E_new)
      
      ! Independently of previous calculation, use Equation (28) in Hogg
      !   (1999) to update t_look
      t_look = t_look  +  0.5d0 * (  1.d0/((1+z_old)*E_old)               &!&
                                   + 1.d0/((1+z_new)*E_new) )
      
      ! Set variables for next pass through loop
      z_old = z_new
      E_old = E_new
    enddo
    
    ! Finally, multiply by appropriate scale factors out front and width of
    !  the trapezoids used for integration
    if( code_stat .lt. 1 )  d_CM = d_H * z/real(num_steps) * d_CM
    t_look = t_H * z/real(num_steps) * t_look
  
  endif
   !-------------------------------------------------------------------------
   ! End calculation of t_look, and possibly also d_CM


!--! Debugging output lines
!comm   write(*,"(A)") "Program finished."
!comm   write(*,2001) "Inputted values were z = ", z_save,                &!&
!comm                 " and d_CM = ", d_save
!comm   write(*,*)
!comm   write(*,1001) "Calculated values: z = ", z
!comm   write(*,1001) "                d_CM = ", d_CM, " Mpc"
!comm   write(*,1001) "              t_look = ", t_look, " years"

return
end subroutine cosmo_calc


!****************************************************************************
!****************************************************************************
subroutine evaluate_A14_integrand(i_tm, j_ph, k_chi, y_in, chi_in, integrand)

! Evaluates the integrand of GS2002's (A14) at the specified time, frequency,
!   y, and chi.
!
! Input arguments
!  1) i_tm: index in tobs_array we're using
!  2) j_ph: index in phot_en_array we're using
!  3) k_chi: index in chi_array we're using
!  4) y_in: value of y we're using
!  5) chi_in: value of chi we're using
! Output argument
!  1) integrand: value of the integrand at that location/time/frequency

use parameters
use constants
use data_arrays

implicit none

  ! Input arguments
integer, intent(in) :: i_tm, j_ph, k_chi
real(kind=8), intent(in) :: y_in, chi_in
  ! Output arguments
real(kind=8), intent(out) :: integrand
  
  ! Local variables
integer :: i_tm_lo, j_ph_lo
real(kind=8) :: t_obs, phot_en_obs, R_zero, t_zero, mu, Gam_emis,         &!&
     little_gam, beta, phot_en_pf, Psyn_ll, Psyn_lh, Psyn_hl, Psyn_hh,    &!&
     tzero_denom, phot_en_denom, Psyn_l, Psyn_h, Psyn_int

  ! Convert array indices to values
  t_obs = tobs_array(i_tm)
  phot_en_obs = phot_en_array(j_ph)
  
   !-------------------------------------------------------------------------
!--! In order to evaluate Psyn, we need to find (1) initial shock conditions
   !   associated with this t_obs/y/chi, and (2) the local plasma frame
   !   photon energy associated with phot_en_array(j_ph)
   !-------------------------------------------------------------------------
  ! Compute initial time this fluid parcel crossed the shock
  !      y = R/R_LOS
  !      chi = (R/R_0)^(4-k)    <== Eq (A8)
  !      R_0 = c*t_0
  R_zero = R_LOS*y_in / chi_in**(1.d0/(4.d0-k_cbm))
  t_zero = R_zero / ccgs
  
  ! mu from Eq (A16)
  mu = 1.d0  -  (1.d0 - chi_in*y_in**(4.d0-k_cbm))                        &!&
               / ( 2.d0*(4.d0-k_cbm) * Gam_LOS**2 * y_in )
  
  ! little_gam from Eq (A9)
  Gam_emis = Gam_LOS * y_in**(0.5d0*(k_cbm-3.d0))
  little_gam = Gam_emis / sqrt( 2.d0 * chi_in )
  
  ! If we are far enough from the shock, the BM solution breaks down and
  !   allows for little_gam < 1.  Should that occur simply exit the
  !   subroutine and return zero emission
  if( little_gam .lt. 1.05d0 ) then
    integrand = 1.d-99
    return
  endif
  
  ! Compute the local photon energy of interest using the Doppler factor
  beta = sqrt( 1.d0 - little_gam**(-2) )
  phot_en_pf = phot_en_obs * little_gam * ( 1.d0 - beta*mu )              &!&
              * (1.d0 + redshift)
   !-------------------------------------------------------------------------
   ! t_zero and phot_en_pf computed
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
!--! Use t_zero and phot_en_pf to locate this point in the precomputed
   !   Psyn_array.  If the point does not fall within the bounds of
   !   Psyn_array, set emission to zero and exit the subroutine
   !-------------------------------------------------------------------------
  i_tm_lo = count( tzero_array .le. t_zero ) - 1
  if( (i_tm_lo .lt. 0) .or. (i_tm_lo .ge. n_ph) ) then
    integrand = 1.d-99
    return
  endif
  ! Check to make sure t_zero is properly bracketed
  if( (t_zero .lt. tzero_array(i_tm_lo)) .or.                             &!&
      (t_zero .ge. tzero_array(i_tm_lo+1)) ) then
    write(*,"(A)") 'ERROR in evaluate_A14_integrand: t_zero not bracketed'
    write(*,"(3x,3es18.7e2)") tzero_array(i_tm_lo), t_zero,               &!&
       tzero_array(i_tm_lo+1)
    write(*,"(A)") 'Stopping program now'
    stop
  endif
  
  j_ph_lo = count( phot_en_array .le. phot_en_pf ) - 1
  if( (j_ph_lo .lt. 0) .or. (j_ph_lo .ge. n_ph) ) then
    integrand = 1.d-99
    return
  endif
  ! Check to make sure phot_en_pf is properly bracketed
  if( (phot_en_pf .lt. phot_en_array(j_ph_lo)) .or.                       &!&
      (phot_en_pf .ge. phot_en_array(j_ph_lo+1)) ) then
    write(*,"(2A)") 'ERROR in evaluate_A14_integrand: phot_en_pf not ',   &!&
       'bracketed'
    write(*,"(3x,3es18.7e2)") phot_en_array(j_ph_lo), phot_en_pf,         &!&
       phot_en_array(j_ph_lo+1)
    write(*,"(A)") 'Stopping program now'
    stop
  endif
   !-------------------------------------------------------------------------
   ! Point located
   !-------------------------------------------------------------------------
  
  
   !-------------------------------------------------------------------------
!--! Perform the bilinear interpolation within Psyn_array.  It is not tri-
   !   linear interpolation because chi is guaranteed to be one of the index
   !   values.
   !-------------------------------------------------------------------------
  Psyn_ll = Psyn_array(j_ph_lo,   i_tm_lo,   k_chi)
  Psyn_lh = Psyn_array(j_ph_lo,   i_tm_lo+1, k_chi)
  Psyn_hl = Psyn_array(j_ph_lo+1, i_tm_lo,   k_chi)
  Psyn_hh = Psyn_array(j_ph_lo+1, i_tm_lo+1, k_chi)
  
  tzero_denom   = 1.d0 / (tzero_array(i_tm_lo+1) - tzero_array(i_tm_lo))
  phot_en_denom = 1.d0 / (phot_en_array(j_ph_lo+1) - phot_en_array(j_ph_lo))
  
  Psyn_l =  (tzero_array(i_tm_lo+1) - t_zero) * tzero_denom * Psyn_ll     &!&
          + (t_zero - tzero_array(i_tm_lo))   * tzero_denom * Psyn_lh
  Psyn_h =  (tzero_array(i_tm_lo+1) - t_zero) * tzero_denom * Psyn_hl     &!&
          + (t_zero - tzero_array(i_tm_lo))   * tzero_denom * Psyn_hh
  
  Psyn_int =  (phot_en_array(j_ph_lo+1) - phot_en_pf)                     &!&
                                               * phot_en_denom * Psyn_l   &!&
            + (phot_en_pf - phot_en_array(j_ph_lo))                       &!&
                                               * phot_en_denom * Psyn_h
   !-------------------------------------------------------------------------
   ! Interpolation complete
   !-------------------------------------------------------------------------
   
   
   !-------------------------------------------------------------------------
!--! Combine everything into the integrand of Eq (A14) and return that as the
   !   output of the subroutine
   !-------------------------------------------------------------------------
  integrand = chi_array(k_chi) * y_in**(10.d0-2.d0*k_cbm) * Psyn_int      &!&
             / ( 1.d0  +  (7.d0-2.d0*k_cbm) * chi_array(k_chi)            &!&
                         * y_in**(4.d0-k_cbm) )**2
   !-------------------------------------------------------------------------
   ! Integrand computed
   !-------------------------------------------------------------------------
   
return
end subroutine evaluate_A14_integrand


!****************************************************************************
!****************************************************************************
subroutine Bisection_method(func_name, i_params, r_params,      &!&
     x_root, l_bound, u_bound)

! Uses bisection to find the root of a function.  Function name is passed to
!   the subroutine, and interface is written below.  Note that the function
!   has a very specific list of arguments -- any additional inputs must be
!   included from modules!
!
! Input arguments
!   1) func_name: name of the function we're trying to find the root of
!   2) i_params: array (length 7) of up to 7 integer parameters to be used
!     in evaluating the function
!   3) r_params: array (length 7) of up to 7 floating point parameters to
!     be used in evaluating the function
!   4) l_bound: minimum allowed value for x
!   5) u_bound: maximum allowed value for x
! Output argument
!   1) x_root: location of the actual root

implicit none

interface
  double precision function func_name(x, i_params, r_params)
    real(kind=8), intent(in) :: x
    integer, intent(in) :: i_params(7)
    real(kind=8), intent(in) :: r_params(7)
  end function
end interface

  ! Input argument (func_name handled in interface above)
real(kind=8), intent(in) :: l_bound, u_bound
integer, dimension(7), intent(in) :: i_params
real(kind=8), dimension(7), intent(in) :: r_params
  ! Output argument
real(kind=8), intent(out) :: x_root

  ! Local variables
real(kind=8), parameter :: tol = 1.d-11
integer, parameter :: max_itrs = log(tol)/log(0.5d0)+10
integer :: i
real(kind=8) :: x_lo, y_lo, x_hi, y_hi, x_mid, y_mid

  ! Compute the function values at the two endpoints, and make sure they
  !   genuinely bracket the root
  x_lo = l_bound
  y_lo = func_name(x_lo, i_params, r_params)
  x_hi = u_bound
  y_hi = func_name(x_hi, i_params, r_params)
  if( (y_lo*y_hi) .gt. 0.d0 ) then
    write(*,"(2A)") 'Bisection_method failed to find root ',&!&
         'because l_bound & u_bound do not bracket it'
    write(*,"(A)") 'Stopping program now'
    stop
  elseif( y_lo .eq. 0.d0 ) then
    x_root = x_lo
    return
  elseif( y_hi .eq. 0.d0 ) then
    x_root = x_hi
    return
  endif
  
  do i = 1, max_itrs
    
    ! Compute the value of f at the midpoint between x_lo and x_hi
    x_mid = 0.5d0 * (x_lo + x_hi)
    y_mid = func_name(x_mid, i_params, r_params)
    
    ! If f(x_mid) is exactly zero, we've found the root.  Otherwise,
    !  If f(x_lo) has same sign as f(x_mid) then x_mid replaces x_lo.
    !  If f(x_hi) has same sign as f(x_mid) then x_mid replaces x_hi.
    if( y_mid .eq. 0.d0 ) then
      x_root = x_mid
      exit
    else
      if( (y_mid*y_hi) .gt. 0.d0 ) then
        x_hi = x_mid
      else
        x_lo = x_mid
      endif
    endif
    
  enddo
  
  x_root = x_mid

return
end subroutine Bisection_method


!****************************************************************************
!****************************************************************************
double precision function GPS_egg(y_input, i_params, r_params)

! Start with and equation for chi and Eq. (A16) of GranotSari2002:
!
!      chi = 1 + 2*(4-k)*Gam_emis^2*(1 - r/y)            [1]
!      mu  = 1 - [1-chi*y^(4-k)]/[2(4-k)*Gam_axis^2*y]   [2]
!
!   r is the radius of the emission site in units of R_LOS, y is the shock
!   radius at time of emission in units of R_LOS, and mu is the fraction of
!   y lying along the line of sight.
! Note that Eq. [1] comes from BM76's Eq. (27), with ct replaced by their
!   Eq. (26), and truncated at order Gam^-2.
!
! Equation [2] can be rewritten
!
!      0 = 1 - mu - [1 - chi*y^(4-k)]/[2*(4-k)*Gam_axis^2*y]   [3]
!
!   Multiply by (2*(4-k) Gam_axis^2 y) to eliminate a division to arrive at
!
!      0 = (2*(4-k) Gam_axis^2 y)*(1 - mu) + chi*y^(4-k) - 1   [4]
!
! Now multiply Eq. [1] by y (again to eliminate a division) and we get
!
!      chi*y = y + 2*(4-k)*Gam_emis^2*(y - r)                  [5]
!
! Assuming that we know Gam_axis and mu, Eqs. [4] and [5] are a system of
!   two equations for the two unknowns chi and y.
!
! The quantities r, y, and mu are assumed to lie between -1 and +1 (with r
!   and y being further restricted to positive values only).
!
! Input arguments:
!   1) y_input: current value being tried for y
!   2) i_params: array (length 7) of up to 7 integer parameters to be used
!     in evaluating the function
!   3) r_params: array (length 7) of up to 7 floating point parameters to
!     be used in evaluating the function
! Output argument:
!   1) GPS_egg: the difference between the right and left sides of
!     Equation [4] above -- obviously we want that to be 0

use parameters

implicit none

  ! Input arguments
real(kind=8), intent(in) :: y_input
integer, dimension(7), intent(in) :: i_params
real(kind=8), dimension(7), intent(in) :: r_params
  ! Output argument is the function name; no need to declare it, apparently

  ! Local variables
integer :: isign
real(kind=8) :: r, mu, Gam_axis, Gam_emis, chiy
  
  
  select case( i_params(1) )
    case( 1 )  ! The typical use case: inside the egg
      r        = r_params(1)
      mu       = r_params(2)
      Gam_axis = r_params(3)
      Gam_emis = sqrt( Gam_axis**2 * y_input**(k_cbm - 3.d0) )
      chiy = y_input  +  2.d0*(4.d0-k_cbm)*(Gam_emis**2) * (y_input - r)
    case( 2 )  ! An alternate use case: along the shell given mu
      mu       = r_params(1)
      Gam_axis = r_params(2)
      chiy = y_input
    case( 3 )  ! Another use case: along shell given R_perp
      isign    = i_params(2)
      r        = r_params(1)  ! R_perp as used here
      mu       = isign * sqrt( y_input**2 - r**2 ) / y_input
      Gam_axis = r_params(2)
      chiy = y_input
    case default
      write(*,"(2A,I0)") 'ERROR in GPS_egg: i_params(1) not a ',          &!&
          'valid choice: ', i_params(1)
      write(*,"(8X,A)") '1: solve for y inside shell'
      write(*,"(8X,A)") '2,3: solve for y at shell'
      write(*,"(A)") 'Stopping program now'
      stop
  end select
  
  GPS_egg = 2.d0*(4.d0-k_cbm)*(Gam_axis**2)*y_input*(1.d0 - mu)           &!&
           +  chiy*y_input**(3.d0-k_cbm)  -  1.d0

return
end function GPS_egg


