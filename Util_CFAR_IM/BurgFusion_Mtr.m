function sig_mtr_rec = BurgFusion_Mtr(sig_mtr, det_map, mod_order)
% Implement data extrapolation for each row of a matrix (e.g., t-f
% diagram).
%
% Parameter:
%   sig_mtr  ---  data matrix of signal 
%   det_map  ---  the detection map of the signal 
%                  (position indicator of (un)available data).
% mod_order  ---  model order used for data extrapolation
%
% Output:
%  sig_mtr_rec ---  the data matrix of recovered signal
%
% Author: Jianping Wang @ MS3, TUDelft
% Date:   May 19, 2021
% $$

[N1,N2] = size(sig_mtr);

sig_mtr_rec = zeros(N1,N2);
for kk = 1:N1
    sig_mtr_rec(kk,:) = BurgFusion_MultiSeg(sig_mtr(kk,:),...
                                            det_map(kk,:),...
                                            mod_order);
end