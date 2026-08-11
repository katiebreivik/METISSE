subroutine comenv_lambda(KW,M0,L,R,MENVD,LAMBDA,id,LAMBF)
    use track_support
    implicit none
    real(dp), intent(in):: M0,L,R,MENVD,LAMBDA
    integer, intent(in) :: KW
    real(dp), intent(out) :: LAMBF

    REAL(dp) :: CELAMF, RZAMSF
    EXTERNAL CELAMF, RZAMSF
    integer, intent(in), optional :: id
    real(dp):: RZAMS, BE, alpha_rec
    integer :: idd
    type(track), pointer :: t

    idd = 1
    if(present(id)) idd = id

    t => tarr(idd)

    if (t% is_he_track) then
        RZAMS = 10.d0**t% tr(i_logR,ZAMS_HE_EEP)
    elseif(kw>= He_MS .and. kw<=He_GB .and. use_sse_NHe)then
        RZAMS = RZAMSF(M0)
    else
        RZAMS = 10.d0**t% tr(i_logR,ZAMS_EEP)
    endif

    ! LAMBDA <= -0.5 selects the MESA-derived binding-energy prescription;
    ! alpha_rec (recombination-energy fraction, Yamaguchi et al. 2024 Eq.
    ! 6-8: E_bind = Egrav + alpha_th*Eth + alpha_rec*Erec, alpha_th=1
    ! always) is then a continuous linear function of LAMBDA: -0.5 -> 0 (no
    ! recombination), -1.0 -> 0, -2.0 -> 1 (full recombination), linear in
    ! between -- keeps existing lambdaf=-1/-2 usage exactly backward
    ! compatible while allowing any intermediate value, e.g. -1.5 -> 0.5.
    ! binding_energy/binding_energy_re bracket the alpha_rec=0/1 endpoints,
    ! so no third stored column is needed for the interpolation.
    if (LAMBDA <= -0.5d0 .and. (i_binding_energy > 0 .or. i_binding_energy_re > 0)) then
        if (i_binding_energy > 0 .and. i_binding_energy_re > 0) then
            alpha_rec = MIN(1.0d0, MAX(0.0d0, -1.0d0 - LAMBDA))
            BE = t%pars%binding_energy + alpha_rec*(t%pars%binding_energy_re - t%pars%binding_energy)
        else if (i_binding_energy > 0) then
            BE = t%pars%binding_energy  ! no with-recombination column available
        else
            BE = t%pars%binding_energy_re  ! no without-recombination column available
        endif
        LAMBF = - (t%pars%mass * (t%pars%mass - t%pars%core_mass)) / &
         ((BE / 3.8d48) * R) !units core_mass in solar, BE in ergs
    else
        LAMBF = CELAMF(KW,M0,L,R,RZAMS,MENVD,LAMBDA)
    endif

    ! BE>0 (envelope formally unbound via recombination, Yamaguchi et al.
    ! 2024 Sec 6.1.2) flips LAMBF negative here -- treat as trivially
    ! ejectable rather than let a negative lambda propagate downstream.
    if (LAMBF < 0.d0) LAMBF = 100.d0
    LAMBF = MIN(100.0d0, LAMBF)

    nullify(t)
end subroutine comenv_lambda
