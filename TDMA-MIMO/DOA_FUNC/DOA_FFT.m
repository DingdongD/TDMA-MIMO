function [PoutFFT] = DOA_FFT(arrData, cfgDOA)
    %% 本文件为基于FFT的DOA/AOA估计
    %% By Xuliang, 20230412
    
    doa_fft = fftshift(fft(arrData, cfgDOA.FFTNum)) * 2 / cfgDOA.FFTNum ;
    PoutFFT = (doa_fft);
    
end