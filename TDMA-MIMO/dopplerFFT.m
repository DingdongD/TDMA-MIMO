function fftOut = dopplerFFT(rangeFFTOut)
    %% 本文件用于实现多普勒维FFT
    %% By Xuliang,20230412
    
    ADCNum = size(rangeFFTOut, 1);
    ChirpNum = size(rangeFFTOut, 2);
    arrNum = size(rangeFFTOut,3);
    fftOut = {};
    % 多普勒维FFT 
    dopplerWin = hanning(ChirpNum)'; % 汉宁窗
    dopplerWin3D = repmat(dopplerWin, ADCNum, 1, arrNum); % 扩充与dopplerData数据一致
    dopplerData = rangeFFTOut .* dopplerWin3D; % 多普勒加窗
    dopplerFFTOut = fftshift(fft(dopplerData, [], 2),2) * 2 * 2 / ChirpNum; % 对多普勒维做FFT【FFT补偿+汉宁窗补偿】 
    fftOut.dopplerFFT = dopplerFFTOut;
end