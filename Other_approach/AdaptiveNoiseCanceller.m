function [sig_interfMit, sig_fft_pos,sig_fft_neg] ...
    = AdaptiveNoiseCanceller(sig, powThreshold, filterLength, marginStep, N_fft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interference mitigation with adaptive noise canceller 
%
% Reference: 
%   F. Jin and S. Cao, "Automotive Radar Interference Mitigation 
%    Using Adaptive Noise Canceller," IEEE Transactions on Vehicular 
%    Technology, vol. 68, no. 4, pp. 3747-3754, 2019.
%
% J.Wang @MS3 TU Delft
% Created: Dec 5, 2019
% email: jianpingwang87@gmail.com or J.Wang-4@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_fft = fft(sig, N_fft);
len_sig = length(sig_fft);

if mod(len_sig,2)==0
    sig_fft_pos = sig_fft(1:len_sig/2);
    sig_fft_neg = sig_fft(len_sig/2+1:end);
    
else
    sig_fft_pos = sig_fft(1:(len_sig+1)/2);
    sig_fft_neg = sig_fft((len_sig+1)/2+1 : end);
end

len_output = length(sig_fft_pos);

ref_int = conj(flip(sig_fft_neg));     %

P_ref = sum(ref_int .* conj(ref_int));

if P_ref > powThreshold
    w0 = zeros(filterLength,1);
    w0(1) = 1;
    
    fi = zeros(len_output,1);
    epsilon = zeros(len_output,1);
    delta_w = 2/(P_ref * marginStep);
    for hh = 1:len_output
        fi = [ref_int(hh); fi(1:filterLength-1)];
        fo = w0' * fi;
        epsilon(hh) = sig_fft_pos(hh) - fo;
        
        P_ref_it = sum(abs(fi).^2);
        delta_w = 2/(P_ref_it * marginStep);
        w0 = w0 + delta_w * fi * conj(epsilon(hh));
    end
    sig_interfMit = epsilon;
else
    sig_interfMit = sig_fft_pos;
end
