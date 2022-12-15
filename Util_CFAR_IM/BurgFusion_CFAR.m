function [sigFuse, sig_interp] = BurgFusion_CFAR(sig, ind_L0, ind_L1, ind_H0, ind_H1, ord_P)
    ii = sqrt(-1);
    
    [N_r, N_c] = size(sig);
    sig_fw = sig(:, ind_L0:ind_L1);
    sig_bw = sig(:, ind_H0:ind_H1);
    
    N_fw = size(sig_fw, 2);
    N_bw = size(sig_bw, 2);
    
    N_extp = ind_H0-ind_L1;
    %foward
    sig_fw_extp = BurgExtrapolate(sig_fw, ord_P, N_extp);
    phase_e_fw = angle(sig_bw(:,1)) - angle(sig_fw_extp(:, end));
    phase_e_p_fw = phase_e_fw / (N_extp-1) * (0:N_extp-1);
    sig_fw_extp(:, N_fw+1:N_fw+N_extp) = sig_fw_extp(:, N_fw+1:N_fw+N_extp) .* exp(ii*phase_e_p_fw);
    sig_fw_extp = sig_fw_extp(:, 1:end-1);
    
    %backward
    sig_bw_extp = BurgExtrapolate( flip(sig_bw,2), ord_P, N_extp);
    phase_e_bw = angle(sig_fw(:, end)) - angle(sig_bw_extp(:, end));
    phase_e_p_bw = phase_e_bw / (N_extp-1) * (0:N_extp-1);
    sig_bw_extp(:, N_bw+1 : N_bw+N_extp) = sig_bw_extp(:, N_bw+1 : N_bw+N_extp) .* exp(ii*phase_e_p_bw);
    sig_bw_extp = sig_bw_extp(:, 1:end-1);
    sig_bw_extp = flip(sig_bw_extp,2);
    
    %interpolation
    p_b = ind_L1+1;
    p_e = ind_H0-1;
    p   = p_b : p_e;
    
    if N_fw>N_bw
        Ic = 0.5*(1+cos(pi*(1+(p-p_b)/(p_b-p_e))));
        r = log(0.5)/log( 0.5*(1 + cos(pi*(1 + N_fw/(N_fw+N_bw))) ) );
        cw = Ic.^r;
    else
        Ic = 1 - 0.5*(1+cos(pi*(1+(p-p_b)/(p_b-p_e))));
        r = log(0.5)/log( 0.5*(1 + cos(pi*(1 + N_bw/(N_fw+N_bw))) ) );
        cw = 1 - Ic.^r;
    end
%     l = round((p_b-1)/((p_b-1)+(N_c-p_e))*(length(p)+1));
%     r = log(0.5)/log(Ic(l));
%     cw = Ic.^r;
    
    cw = repmat(cw,N_r,1);
    
    sig_interp = (1-cw).*sig_fw_extp(:,N_fw+1:N_fw+N_extp-1)...
        +cw.*sig_bw_extp(:,1:N_extp-1);
    
    sigFuse = [sig(:, 1:ind_L0-1), sig_fw, sig_interp, sig_bw, sig(:, ind_H1+1:end)];
end