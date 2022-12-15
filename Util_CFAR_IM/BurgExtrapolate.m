function sig_extp = BurgExtrapolate(sig, ord_P, N_extp)
    
    [N_r, N_c] = size(sig);
    
    sig_extp = [sig, zeros(N_r, N_extp)];
    a = arburg(sig.', ord_P);
    for kk = N_c+1:N_c+N_extp
        est_n = diag( a(:,2:end) * (sig_extp(:,kk-1:-1:kk-ord_P).' ) );
        sig_extp(:,kk) = -est_n;
    end
end