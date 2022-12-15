function sig_dechirp_LPF = LPF_AftDechirp(Range_max, T_sw, BW, f_s, sig_dechirp)

R_p = 3;                                 % R_p-passband ripple;
R_s = 20;                                % R_s-attenuation in the stopband
fb_max    = (2*BW*Range_max)/(3e8*T_sw);   % fb_max =  13.3 MHz cutoff frequency
Wp = fb_max/(f_s/2);
f_stopband = fb_max+2e6;
Ws = f_stopband/(f_s/2);                 % f_stop = 15.3 MHz
% n-lowest order; Wn-scalar of cutoff frequency;
[n,Wn] = buttord (Wp,Ws,R_p,R_s);
[b,a]  = butter(n,Wn);

sig_dechirp_LPF = filter(b,a,sig_dechirp);