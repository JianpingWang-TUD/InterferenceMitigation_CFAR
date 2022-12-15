function sig_AftIM = IM_Wdenoise(sig, level)
% wavelet denoise based interference mitigation for FMCW radars
%
% Parameter:
%   sig   --- Interference-contaminated signal
%   level --- Level of wavelet decomposition
%   sig_AftIM --- Signal after interference mitigation
%
% Reference:
%    S. Lee, J. Lee, and S. Kim, "Mutual Interference Suppression Using
%       Wavelet Denoising in Automotive FMCW Radar Systems," IEEE 
%       T-ITS, pp. 1-11, 2019.
%
% Author: Jianping Wang
% Email: J.Wang-4@tudelft.nl
%

if size(sig,1) ==1
    sig = sig.' ;
end

xd_real = wdenoise(real(sig),level,'Wavelet','haar',...
                            'DenoisingMethod','SURE',...
                            'ThresholdRule','hard',...
                            'NoiseEstimate','LevelDependent');
xd_imag = wdenoise(imag(sig),level,'Wavelet','haar',...
                            'DenoisingMethod','SURE',...
                            'ThresholdRule','hard',...
                            'NoiseEstimate','LevelDependent');   
sig_AftIM = sig - ( xd_real + 1i*xd_imag );


