function y = ampCorrect(sig_mtr, det_map)
%
% Parameter:
%   sig_mtr --- t-f domain signal matrix
%   det_map --- the detection map
%
% Author: Jianping Wang
%

[num_f, ~] = size(sig_mtr);
amp_mtr = abs(sig_mtr);
ph_mtr = angle(sig_mtr); 
y = sig_mtr;
for kk = 1:num_f
    ind = det_map(kk,:);
    sig_amp = mean(amp_mtr(kk,~ind));
    y(kk,ind) = sig_amp * exp(1i*ph_mtr(kk,ind));
end