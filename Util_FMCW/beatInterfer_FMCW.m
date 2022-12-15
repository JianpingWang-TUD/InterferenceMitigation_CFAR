function sig_Interfer = beatInterfer_FMCW(amp_intf, fc_intf, fr_intf, T_sw_intf, t_d_intf,...
                            t, fc, fr, T_sw, freqCut_LP)
% beatInterfer_FMCW generates the interference signal after dechirping and
% low-pass filtering in FMCW radar systems
%
% Parameters:
%   amp_intf --- Amplitudes of interference, scalar/vector
%   fc_intf  --- center frequency of interference, scalar/vector
%   fr_intf  --- chirp rate of interference, scalar/vector
%   T_sw_intf --- Duration of interference, scalar/vector
%   t_d_intf --- delay time of interference relative to (dechirping) 
%                 reference signal, scalar/vector
%   t        --- time samples, vector
%   fc       --- center frequency of dechirping reference signal
%   fr       --- chirp rate of dechirping ref signal
%   T_sw     --- duration of dechirping ref signal
%   freqCut_LP --- cut frequency of low-pass filter after dechirping
%
% author: Jianping Wang @  MS3, TU Delft
%

numIntf = length(amp_intf);

sig_Interfer = zeros(1,length(t));

for kk = 1:numIntf
    f_const = fc_intf(kk) - fc - fr_intf(kk) * (t_d_intf(kk) + 0.5*T_sw_intf(kk))...
               + 0.5*fr*T_sw;
    f_lin = (fr_intf(kk) - fr) * t;
    
    ph_const = -fc_intf(kk)*(T_sw_intf(kk)/2 + t_d_intf(kk)) + fc * T_sw/2 ...
                + 0.5 * fr_intf(kk) * ( 0.5*T_sw_intf(kk) + t_d_intf(kk) )^2 ...
                - 0.5 * fr * ( 0.5*T_sw )^2; 
            
    bf_I = f_const + f_lin;
    I = bf_I >= -freqCut_LP & bf_I <= freqCut_LP;
    
    sig_I_temp = amp_intf(kk) * rectpuls(t-T_sw_intf(kk)/2-t_d_intf(kk), T_sw_intf(kk))...
        .*exp(1i*2*pi* (f_const.*t + 0.5*f_lin.*t + ph_const ) );
    sig_Interfer = sig_Interfer + I.* sig_I_temp;
end