function  sig_rec = BurgFusion_MultiSeg(sig, det_map, mod_ord)
% Using the Burg method to extrapolate data in the gap, front and end.
%
% Parameters:
%   sig    ---   the data vector of signal 
%  det_map ---   the detected map (position indicators) of the 
%                 available/unavailable data
%  mod_ord ---   model order
%
% Output:
%   sig_rec ---  the data vector of recovered signal
%
% Author: Jianping Wang @ MS3, TUDelft
% Date:   May 19, 2021
% $$

[ind_ava, ~] = data_segmentation(det_map);

%==================================
%eliminate the data segments with length < model order
Num_ind_ava = ind_ava(:,2) - ind_ava(:,1) + 1;
ind_ava(Num_ind_ava<mod_ord+1,:) = [];
%===================================


d1_ind_ava = size(ind_ava,1);

sig_rec = sig;
for kk = 1:d1_ind_ava-1
    ind_L0 = ind_ava(kk,1);
    ind_L1 = ind_ava(kk,2);
    ind_H0 = ind_ava(kk+1,1);
    ind_H1 = ind_ava(kk+1,2);
    [~, sig_interp_tmp] = BurgFusion_CFAR(sig, ind_L0, ind_L1, ind_H0, ind_H1, mod_ord);
    ind_tmp = ind_ava(kk,2)+1: ind_ava(kk+1,1)-1;
    sig_rec(ind_tmp) = sig_interp_tmp;
end

% backward extrapolation
if ind_ava(1,1) > 1
    sig_temp_flip = sig( ind_ava(1,2) : -1 : ind_ava(1,1) );
    N_extrap = ind_ava(1,1)-1;
    sig_extrap = BurgExtrapolate( sig_temp_flip, mod_ord, N_extrap );
    sig_rec(1:ind_ava(1,2)) = flip(sig_extrap);
end

% forward extrapolation
if ind_ava(end,2) < length(sig)
    sig_temp = sig( ind_ava(end,1):1:ind_ava(end,2) );
    N_extrap = length(det_map) - ind_ava(end,2);
    sig_extrap = BurgExtrapolate( sig_temp, mod_ord, N_extrap);
    sig_rec(ind_ava(end,1):end) = sig_extrap;
end