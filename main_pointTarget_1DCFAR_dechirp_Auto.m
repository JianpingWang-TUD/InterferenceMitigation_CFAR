%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFAR-based interference suppression
%
% J.Wang @MS3, TU Delft, Dec 5, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear;  clc;

flag_plot = 0;   % flag to turn on/off the figure plot
Ftsz = 12;
fig_dir = '/result/';
%% FMCW radar
ii = sqrt(-1);
c      = 3e8;
f_c    = 77e9; %3e9;             % center frequency 3GHz
T_sw   = 100e-6;%500e-6;          % Sweep time 200 us
BW     = 600e6;%40e6;
Num    = 1;
sweep_slope = -BW/T_sw;
R_max   = 250;%8e3;
tau_max = 2*R_max/c;
fb_max = abs(sweep_slope*tau_max);

BW_I          = BW;
T_sw_I        = T_sw/Num;
sweep_slope_I = BW_I/T_sw_I;
sweep_slope_I = sweep_slope;

f_s           = 40e6;          % Sampling frequency
Npoint        = f_s*T_sw;
Npoint_I      = f_s*T_sw_I;

P_tx = 1;                % Transmit Power 1 Watt

W_len = 256;
Win_len = hamming(W_len,'periodic');
hop   = 4;
noverlap = W_len-hop;

mul_factor_fft = 1;
% SNR   = inf;              % dB
SNR = 5; %15
%% Transmitter
t       = (0 : 1 : Npoint-1)/f_s;
N_t = length(t);
amp     = sqrt(P_tx);
sig_Tx  =  FMCW_sweep(t, T_sw, 0, f_c, -sweep_slope, amp);

%% Receiver
d_tar = [30, 80, 150, 153];%[2e3,  3e3,  5e3 ];
N_tar = length(d_tar);
scat_coeff_tar = [1, 0.1, 0.7, 0.7] .* exp( 1i*2*pi*rand(1, N_tar) );

sig_Rx = beatSig_FMCW(scat_coeff_tar, d_tar, t, f_c, T_sw, sweep_slope, c);

sigN_rx = awgn(sig_Rx,SNR,'measured');
% sigN_rx = sig_Rx;
%% FMCW Interference
amp_intf = [10, 10, 15];
fr_intf = [1.1*sweep_slope_I, -2*sweep_slope_I, 0.9*sweep_slope_I ];
fc_intf = [f_c, f_c, 1*f_c];
T_sw_intf = [T_sw_I, T_sw_I, T_sw_I];
t_d_intf = [0e-6, -50e-6, -2e-6];%[10e-6, 100e-6,  -150e-6];

sig_I = beatInterfer_FMCW(amp_intf, fc_intf, fr_intf, T_sw_intf, t_d_intf,...
    t, f_c, sweep_slope, T_sw, fb_max);

%% full signal synthesis
I = rectpuls(t-T_sw/2-tau_max, T_sw)>0.5;    %range truncation
sig_Rx_trc = sig_Rx(I);
sigN_rx = awgn(sig_Rx_trc,SNR,'measured');
sig_full = sigN_rx+sig_I(I);
t = t(I);

SNR_0 = db( norm(sig_Rx_trc) / norm(sigN_rx-sig_Rx_trc) );
SINR_0 = db( norm(sig_Rx_trc) / norm(sig_full-sig_Rx_trc) );  % initial SINR


hf = figure;
plot(t*1e6, real(sig_full))
grid on
axis tight
xlabel('Time [\mus]','fontsize',Ftsz);
ylabel('Amplitude','fontsize',Ftsz)
title('Real part','fontsize',Ftsz)
if flag_plot == 1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_' num2str(SNR) 'dB.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_' num2str(SNR) 'dB.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_' num2str(SNR) 'dB.fig'])
end

sig_full_fft = fft(sig_full);
%====================================
sigN_ref = fft(sigN_rx);
%====================================
sig_fft_max = max([max(abs(sig_full_fft)), max(abs(sigN_ref))]);
% sig_full_fft_nor = abs(sig_full_fft)/max(abs(sig_full_fft));
sig_full_fft_nor = abs(sig_full_fft)/sig_fft_max;
sigN_ref_nor     = abs(sigN_ref)/sig_fft_max;


N_sig = length(sig_full);
r_axis = f_s/N_sig*(0:N_sig-1)/abs(sweep_slope)*c/2;

hf = figure
plot(r_axis, db(sigN_ref_nor),'r-',...
     r_axis, db(sig_full_fft_nor),'b--')
grid on
axis([0 250 -50 10]);
xlabel('Range [km]', 'fontsize', Ftsz)
ylabel('Normalized amplitude [dB]', 'fontsize', Ftsz)
title('Range profile', 'fontsize', Ftsz)
legend('Interference-free','Interfered','Location','southeast')
if flag_plot==1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3]);
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RP_' num2str(SNR) 'dB.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RP_' num2str(SNR) 'dB.png'])
    saveas(hf, [fig_dir 'pt_RP_' num2str(SNR) 'dB.fig'])
end

%% Interference mitigation
% [sig_TF,f_slice,t_frame] =  spectrogram( (sig_dechirp_LPF),W_len,noverlap,W_len,f_s, 'centered');
[sig_TF,f_slice,t_frame] =  stft(sig_full,f_s,'Window',...
                                            Win_len,'OverlapLength',...
                                            noverlap,'FFTLength',...
                                            mul_factor_fft*W_len);
sig_TF_norm = sig_TF/max(sig_TF(:));
hfig = figure;
imagesc(t_frame*1e6,f_slice/1e6,db(sig_TF_norm))
colorbar;
% colormap('jet')
set(gca,'YDir','normal');
set(get(colorbar,'title'),'string','dB');
xlabel('Time [\mus]', 'FontSize', Ftsz)
ylabel('Frequency [MHz]', 'FontSize', Ftsz)
title('t-f spectrum of beat signal', 'FontSize', Ftsz)
% xlim([40,max(t_frame)*1e6])
%ylim([0,40])
caxis([-80,0])

if flag_plot==1
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])    
print(hfig,'-dpng', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB.png']);
print(hfig,'-dpdf', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB.pdf']);
saveas(hfig, [fig_dir 'pt_TF_' num2str(SNR) 'dB.fig'])
end
% figure
% plot(db(abs(sig_TF(70,:))))

%% Spectrogram
[Num_row,Num_col] = size(sig_TF);
sig_specgram = abs(sig_TF).^2;

%% 1-D CFAR detector
cfar = phased.CFARDetector('NumTrainingCells', 150, 'NumGuardCells',50,...
                           'ThresholdFactor','auto',...
                           'ProbabilityFalseAlarm',1e-4); %450 110 1e-4
det_map = cfar(sig_specgram.', 1:Num_col);
det_map_dlted0 = maskDilate(det_map.', 12);  % 15
%---------------------------------
cfar1 = phased.CFARDetector('NumTrainingCells', 250, 'NumGuardCells',130,...
                            'ThresholdFactor','auto',...
                            'ProbabilityFalseAlarm',1e-4); %450 110 1e-4
det_map1 = cfar1(sig_specgram.', 1:Num_col);
det_map_dlted1 = maskDilate(det_map1.', 12);  % 15
%---------------------------------
cfar2 = phased.CFARDetector('NumTrainingCells', 60, 'NumGuardCells',30,...
                            'ThresholdFactor','auto',...
                            'ProbabilityFalseAlarm',1e-4); %450 110 1e-4
det_map2 = cfar2(sig_specgram.', 1:Num_col);
det_map_dlted2 = maskDilate(det_map2.', 12);  % 15
%----------------------------------------------

hf = figure
imagesc(t_frame*1e6,f_slice/1e6, det_map.');
set(gca, 'YDir', 'normal')
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Frequency [MHz]', 'fontsize', Ftsz)
title('Map of interference detection', 'fontsize', Ftsz)
if flag_plot==1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_DetMap_' num2str(SNR) 'dB.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_DetMap_' num2str(SNR) 'dB.png'])
    saveas(hf, [fig_dir 'pt_DetMap_' num2str(SNR) 'dB.fig'])
end


det_map_dlted = det_map_dlted0;% | det_map_dlted1 | det_map_dlted2;

hf = figure
imagesc(t_frame*1e6,f_slice/1e6,det_map_dlted)
set(gca, 'YDir','normal')
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Frequency [MHz]', 'fontsize', Ftsz)
title('Dilated mask', 'fontsize', Ftsz)
if flag_plot==1
set(gcf, 'PaperUnits','inches','PaperSize',[4,3],'PaperPosition',[0 0 4 3])
print(hf, '-dpdf', '-r300', [fig_dir 'pt_CFAR_mask_' num2str(SNR) 'dB.pdf'])
print(hf, '-dpng', '-r300', [fig_dir 'pt_CFAR_mask_' num2str(SNR) 'dB.png'])
saveas(hf, [fig_dir 'pt_CFAR_mask_' num2str(SNR) 'dB.fig'])
end
%% Zeroing
sig_TF_zero = sig_TF;
sig_TF_zero(det_map_dlted>0.5)=0;

sig_TF_zero_norm = sig_TF_zero/max(sig_TF_zero(:));
hf = figure
imagesc(t_frame*1e6,f_slice/1e6, db(sig_TF_zero_norm))
% colormap('jet')
set(gca,'ydir','normal')
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Frequency [MHz]', 'fontsize', Ftsz)
title('t-f diagram by using CFAR-Z', 'fontsize', Ftsz)
if flag_plot==1
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
print(hf, '-dpdf', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarZ.pdf'])
print(hf, '-dpng', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarZ.png'])
saveas(hf, [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarZ.fig'])
end

%% Amplitude Correction
sig_TF_AC = ampCorrect(sig_TF, det_map_dlted);   % amplitude correction

sig_TF_AC_norm = sig_TF_AC/max(sig_TF_AC(:));
hf = figure
imagesc(t_frame*1e6,f_slice/1e6, db(sig_TF_AC_norm))
% colormap('jet')
set(gca,'ydir','normal')
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Frequency [MHz]', 'fontsize', Ftsz)
title('t-f diagram by using CFAR-AC', 'fontsize', Ftsz)
if flag_plot==1
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
print(hf, '-dpdf', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarAC.pdf'])
print(hf, '-dpng', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarAC.png'])
saveas(hf, [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarAC.fig'])
end

%% CFAR-Burg
sig_TF_Burg = BurgFusion_Mtr(sig_TF, det_map_dlted, 2);  %Burg fusion with multi-segments

sig_TF_Burg_norm = sig_TF_Burg/max(sig_TF_Burg(:));
hf = figure
imagesc(t_frame*1e6,f_slice/1e6, db(sig_TF_Burg_norm))
% colormap('jet')
set(gca,'ydir','normal')
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Frequency [MHz]', 'fontsize', Ftsz)
title('t-f diagram by using CFAR-Burg', 'fontsize', Ftsz)
if flag_plot==1
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
print(hf, '-dpdf', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarBurg.pdf'])
print(hf, '-dpng', '-r300', [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarBurg.png'])
saveas(hf, [fig_dir 'pt_TF_' num2str(SNR) 'dB_cfarBurg.fig'])
end

%% ISTFT
% zeroing
sig_TF_zero(1:mul_factor_fft*W_len/2,:)=0;  %suppress the negative spectrum
[sig_CFAR_Z,t_cfar] =  istft(sig_TF_zero,f_s,'Window',Win_len,...
                                            'OverlapLength',noverlap,...
                                            'FFTLength',mul_factor_fft*W_len);

t_cfar_d = t_cfar + tau_max;
hf = figure
plot(t_cfar_d*1e6, real(sig_CFAR_Z))
grid on
axis tight
ylim([-3.5 3.5])
xlabel('Time [\mus]',  'fontsize', Ftsz)
ylabel('Amplitude', 'fontsize', Ftsz)
title('Beat signal recovered with CFAR-Z', 'fontsize', Ftsz)
if flag_plot==1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarZ.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarZ.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarZ.fig'])
end

% amplitude correction
sig_TF_AC(1:mul_factor_fft*W_len/2,:)=0;  %suppress the negative spectrum
[sig_CFAR_AC,t_cfar] =  istft(sig_TF_AC,f_s,'Window',Win_len,...
                                            'OverlapLength',noverlap,...
                                            'FFTLength',mul_factor_fft*W_len);
hf = figure
plot(t_cfar_d*1e6, real(sig_CFAR_AC))
grid on
axis tight
ylim([-3.5 3.5])
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Amplitude', 'fontsize', Ftsz)
title('Beat signal recovered with CFAR-AC', 'fontsize', Ftsz)
if flag_plot==1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarAC.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarAC.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarAC.fig'])
end

% Burg
sig_TF_Burg(1:mul_factor_fft*W_len/2,:)=0;  %suppress the negative spectrum
[sig_CFAR_Burg,t_cfar] =  istft(sig_TF_Burg,f_s,'Window',Win_len,...
                                            'OverlapLength',noverlap,...
                                            'FFTLength',mul_factor_fft*W_len);
hf = figure
plot(t_cfar_d*1e6, real(sig_CFAR_Burg))
grid on
axis tight
ylim([-3.5 3.5])
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Amplitude', 'fontsize', Ftsz)
title('Beat signal recovered with CFAR-Burg', 'fontsize', Ftsz)
if flag_plot==1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarBurg.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarBurg.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_cfarBurg.fig'])
end
                                        
%======================================================
%% Interference mitigation by wavelet denoise
sig_wd = IM_Wdenoise(sig_full, 5);

% close all

hf = figure
plot(t*1e6, real(sig_wd));
grid on
axis tight
ylim([-3.5 3.5])
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Amplitude', 'fontsize', Ftsz)
title('Beat signal recovered with WD', 'fontsize', Ftsz)
if flag_plot == 1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_WD.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_WD.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_WD.fig'])
end

% recovered signals by all the methods
hf = figure
plot(t*1e6, real(sigN_rx),'k-',...
     t*1e6, real(sig_wd), 'm-',...
     t_cfar_d*1e6, real(sig_CFAR_Z),'b-',...
     t_cfar_d*1e6, real(sig_CFAR_AC),'r--',...
     t_cfar_d*1e6, real(sig_CFAR_Burg),'g-.');
grid on
axis tight
ylim([-7, 8])
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Amplitude', 'fontsize', Ftsz)
title('Recovered beat signals', 'fontsize', Ftsz)
legend('ref','WD','CFAR-Z','CFAR-AC','CFAR-Burg','Location','south','NumColumns',3)

ind_inset = t*1e6>=44 & t*1e6<=45;
t_inset = t(ind_inset);
ind_inset_cfar = t_cfar_d*1e6>=44 & t_cfar_d*1e6<=45;
t_cfar_d_inset = t_cfar_d(ind_inset_cfar);
axes('position', [0.3 0.72 0.4 0.2])
plot(t_inset*1e6, real(sigN_rx(ind_inset)),'k-',...
     t_inset*1e6, real(sig_wd(ind_inset)), 'm-',...
     t_cfar_d_inset*1e6, real(sig_CFAR_Z(ind_inset_cfar)),'b-',...
     t_cfar_d_inset*1e6, real(sig_CFAR_AC(ind_inset_cfar)),'r--',...
     t_cfar_d_inset*1e6, real(sig_CFAR_Burg(ind_inset_cfar)),'g-.')
grid on
ylim([-3.2 3.2])

if flag_plot == 1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_Recv.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_Recv.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_Recv.fig'])
end

%% Range profile
n_ft = 2^(nextpow2(length(sig_wd)) + 2);
r_axis_Aft = f_s/n_ft*(0:n_ft-1)/abs(sweep_slope)*c/2;

%================ANC=================
[rp_ANC, sig_fft_pos,sig_fft_neg] ...
    = AdaptiveNoiseCanceller(sig_full, 10, 80, 1.8, n_ft);
%


rp_CFAR_Z  = fft(sig_CFAR_Z, n_ft);
rp_CFAR_AC = fft(sig_CFAR_AC, n_ft);
rp_CFAR_Burg = fft(sig_CFAR_Burg, n_ft);
rp_wd      = fft(sig_wd, n_ft);

max_rp = max(abs( [ rp_CFAR_Z; rp_CFAR_AC; rp_wd ; rp_ANC] ));

rp_CFAR_Z_norm  = db( abs(rp_CFAR_Z) / max_rp );
rp_CFAR_AC_norm = db( abs(rp_CFAR_AC) / max_rp );
rp_wd_norm      = db( abs(rp_wd) / max_rp ); 
rp_ANC_norm     = db( abs(rp_ANC) / max_rp);
rp_CFAR_Burg_norm = db( abs(rp_CFAR_Burg) / max_rp);

hf = figure
plot(r_axis_Aft(1:n_ft/2), rp_ANC_norm, 'g--',...
     r_axis_Aft, rp_wd_norm,'k-',...   
     r_axis_Aft, rp_CFAR_Z_norm,'r--',...
     r_axis_Aft, rp_CFAR_AC_norm,'b-.',...
     r_axis_Aft, rp_CFAR_Burg_norm,'m:')
axis([0 250 -50 10]);
grid on
legend('ANC','WD','CFAR-Z','CFAR-AC','CFAR-Burg','Location','south','NumColumns',2)
xlabel('Range [km]', 'fontsize', Ftsz)
ylabel('Normalized amplitude [dB]', 'fontsize', Ftsz)
title('RP after IM', 'fontsize', Ftsz)


ind_inset = r_axis_Aft>=29 & r_axis_Aft<=31;
r_axis_inset = r_axis_Aft(ind_inset);
axes('position', [0.3 0.65 0.25 0.25])
plot(r_axis_inset, rp_ANC_norm(ind_inset),'g--',...
     r_axis_inset, rp_wd_norm(ind_inset),'k-',...
     r_axis_inset, rp_CFAR_Z_norm(ind_inset),'r--',...
     r_axis_inset, rp_CFAR_AC_norm(ind_inset),'b-.',...
     r_axis_inset, rp_CFAR_Burg_norm(ind_inset),'m:')
grid on
ylim([-20 2])

ind_inset = r_axis_Aft>=148.5 & r_axis_Aft<=154.5;
r_axis_inset = r_axis_Aft(ind_inset);
axes('position', [0.66 0.65 0.23 0.25])
plot(r_axis_inset, rp_ANC_norm(ind_inset),'g--',...
     r_axis_inset, rp_wd_norm(ind_inset),'k-',...
     r_axis_inset, rp_CFAR_Z_norm(ind_inset),'r--',...
     r_axis_inset, rp_CFAR_AC_norm(ind_inset),'b-.',...
     r_axis_inset, rp_CFAR_Burg_norm(ind_inset),'m:')
grid on
ylim([-25 0])

if flag_plot==1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RP_AftSup_' num2str(SNR) 'dB.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RP_AftSup_' num2str(SNR) 'dB.png'])
    saveas(hf, [fig_dir 'pt_RP_AftSup_' num2str(SNR) 'dB.fig'])
end

%% Performance evaluation
len_cfar_z = length(sig_CFAR_Z);
% SINR
SINR_cfar_z = db(norm(sig_Rx_trc(1:len_cfar_z)) ...
                 / norm(sig_CFAR_Z - sig_Rx_trc(1:len_cfar_z).' ) );
SINR_cfar_ac = db( norm(sig_Rx_trc(1:len_cfar_z) ...
                 / norm(sig_CFAR_AC - sig_Rx_trc(1:len_cfar_z).' ) ) );
SINR_cfar_burg = db( norm(sig_Rx_trc(1:len_cfar_z) ...
                 / norm(sig_CFAR_Burg - sig_Rx_trc(1:len_cfar_z).' ) ) );             
SINR_wd = db( norm(sig_Rx_trc) / norm(sig_wd - sig_Rx_trc.') );

sig_ANC = ifft(rp_ANC, n_ft);
sig_ANC = sig_ANC(1:N_sig);

SINR_anc = db( norm(sig_Rx_trc) / norm(sig_ANC - sig_Rx_trc.' ) );

% correlation coefficient
rho_cfar_z = sig_CFAR_Z' * sig_Rx_trc(1:len_cfar_z).'/norm(sig_CFAR_Z)/norm(sig_Rx_trc(1:len_cfar_z));
rho_cfar_ac = sig_CFAR_AC' * sig_Rx_trc(1:len_cfar_z).'/norm(sig_CFAR_AC) /norm(sig_Rx_trc(1:len_cfar_z));
rho_cfar_burg = sig_CFAR_Burg' * sig_Rx_trc(1:len_cfar_z).'/norm(sig_CFAR_Burg) /norm(sig_Rx_trc(1:len_cfar_z));
rho_wd      = sig_wd' * sig_Rx_trc.'/norm(sig_wd) / norm(sig_Rx_trc);
rho_anc     = sig_ANC' * sig_Rx_trc.' / norm(sig_ANC) / norm(sig_Rx_trc);

[SINR_cfar_z, SINR_cfar_ac,SINR_cfar_burg, SINR_wd, SINR_anc]
rho =[rho_cfar_z, rho_cfar_ac, rho_cfar_burg, rho_wd, rho_anc];
rho_amp = abs(rho)
rho_phase = angle(rho)
% rho = rho_amp .* exp(ii*rho_phase)

close all

hf = figure
plot(t*1e6, real(sig_ANC));
grid on
axis tight
% ylim([-2.5 2.5])
xlabel('Time [\mus]', 'fontsize', Ftsz)
ylabel('Amplitude', 'fontsize', Ftsz)
title('Beat signal recovered with ANC', 'fontsize', Ftsz)
if flag_plot == 1
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
    print(hf, '-dpdf', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_ANC.pdf'])
    print(hf, '-dpng', '-r300', [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_ANC.png'])
    saveas(hf, [fig_dir 'pt_RawSig_RealPart_AftSup_' num2str(SNR) 'dB_ANC.fig'])
end



if flag_plot==1
save([fig_dir 'simu_pt_data_' num2str(SNR) '.mat'],...
    'sig_Rx_trc','sig_CFAR_AC','sig_CFAR_Z', 'sig_ANC','sig_wd',...
    'SINR_cfar_ac','SINR_cfar_z','SINR_anc','SINR_wd',...
    'rho_cfar_ac', 'rho_cfar_z', 'rho_anc', 'rho_wd',...
    'rp_CFAR_AC', 'rp_CFAR_Z', 'rp_ANC', 'rp_wd',...
    'SNR', 'SINR_0', 'r_axis_Aft', 'n_ft', 't');
end
%% restore the search directory
path = oldpath;










