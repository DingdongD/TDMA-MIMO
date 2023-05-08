function fftOut = rdFFT(adcData, IQFlag)
    %% 本文件用于实现距离维、多普勒维FFT
    %% By Xuliang,20230412
    
    ADCNum = size(adcData, 1);
    ChirpNum = size(adcData, 2);
    arrNum = size(adcData,3);
    fftOut = {};
    
    if IQFlag
        % 距离维FFT
        rangeWin = hanning(ADCNum); % 距离加窗
        rangeWin3D = repmat(rangeWin, 1, ChirpNum, arrNum); % 扩充与rangeData数据一致
        rangeData = adcData .* rangeWin3D ; % 距离维加窗
        rangeFFTOut = fft(rangeData, [], 1) * 2 * 2 / ADCNum; % 对距离维做FFT【FFT补偿+汉宁窗补偿】
        
        % 多普勒维FFT 
        dopplerWin = hanning(ChirpNum)'; % 汉宁窗
        dopplerWin3D = repmat(dopplerWin, ADCNum, 1, arrNum); % 扩充与dopplerData数据一致
        dopplerData = rangeFFTOut .* dopplerWin3D; % 多普勒加窗
        dopplerFFTOut = fftshift(fft(dopplerData, [], 2),2) * 2 * 2 / ChirpNum; % 对多普勒维做FFT【FFT补偿+汉宁窗补偿】 
    else
        % 采样I路时需要注意目标距离不能超出最大约束目标距离 否则会出现距离模糊
        % 单路信号需要抛掉一半信号

        % 距离维FFT
        rangeWin = hanning(ADCNum); % 汉宁窗
        rangeWin3D = repmat(rangeWin, 1, ChirpNum, arrNum); % 扩充与rangeData数据一致
        rangeData = adcData .* rangeWin3D ; % 距离维加窗
        rangeFFTOut = fft(rangeData, [], 1) * 2 * 2 / ADCNum; % 对距离维做FFT【FFT补偿+汉宁窗补偿】
        rangeFFTOut = rangeFFTOut(1:end/2, :, :);
        
        % 多普勒维FFT 
        dopplerWin = hanning(ChirpNum)'; % 汉宁窗
        dopplerWin3D = repmat(dopplerWin, ADCNum / 2, 1, arrNum); % 扩充与dopplerData数据一致
        dopplerData = rangeFFTOut .* dopplerWin3D; % 多普勒加窗
        dopplerFFTOut = fftshift(fft(dopplerData, [], 2),2) * 2 * 2 / ChirpNum; % 对多普勒维做FFT【FFT补偿+汉宁窗补偿】

    end
    fftOut.rangeFFT = rangeFFTOut;
    fftOut.dopplerFFT = dopplerFFTOut;
end