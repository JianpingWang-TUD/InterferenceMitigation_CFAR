function sig = FMCW_sweep(t, T_sw, t_delay, f_c, sweep_slope, amp)
% This function generates the FMCW signal.
%
sig  =  amp * rectpuls(t-T_sw/2-t_delay, T_sw)...
    .*exp(1i*2*pi*(f_c + 0.5*sweep_slope*(t-T_sw/2-t_delay)).*(t-T_sw/2-t_delay));