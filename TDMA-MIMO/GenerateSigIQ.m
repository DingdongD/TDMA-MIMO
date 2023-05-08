function [signal] = GenerateSigIQ(tarOut, cfgOut)
    %% 本文件用于双路ADC-IQ信号
    %% By Xuliang, 20230411
    tarNum = tarOut.nums; % 目标数目
    numTx = cfgOut.numTx; % 发射天线数目
    numRx = cfgOut.numRx; % 接收天线数目
    
    ADCNum = cfgOut.ADCNum; % ADC采样数目
    Frame = cfgOut.Frame; % 帧数
    ChirpNum = cfgOut.ChirpNum; % 每帧发射Chirp数目
    TotalChirpNum = Frame * ChirpNum; % 总Chirp数目
    
    % 信号参数
    fc = cfgOut.fc; % 载频 Hz
    fs = cfgOut.fs; % ADC采样频率 Hz 
    Ramptime = cfgOut.Ramptime; % 工作时间
    Idletime = cfgOut.Idletime; % 空闲时间
    Slope = cfgOut.Slope; % Chirp斜率
    validB = cfgOut.validB; % 实际有效带宽
    totalB = cfgOut.totalB; % 完整带宽
    Tc = cfgOut.Tc; % 单Chirp周期
    ValidTc = cfgOut.ValidTc; % 单Chirp内ADC有效采样时间
    
    % 硬件参数
    Pt = cfgOut.Pt; % dbm 发射功率
    Fn = cfgOut.Fn; % 噪声系数
    Ls = cfgOut.Ls;  % 系统损耗
    
    % 天线阵列
    antennaPhase = cfgOut.antennaPhase; % 天线相位
    virtual_array = cfgOut.virtual_array; % 虚拟阵列 这里的虚拟阵列是struct virtual_arr中包含了发射方位相位和俯仰相位 接收天线 发射天线顺序 
    arr = virtual_array.virtual_arr;
    
    % 初阶版本顺序发射可行
%     arr = cfgOut.array; % 虚拟阵元排列
    arrNum = size(arr, 1); % 虚拟阵列数目
    arrdelay = zeros(1, arrNum) + Tc * reshape(repmat([0:numTx-1], numRx,1), 1, arrNum); % [0 0 0 0 1 1 1 1 2 2 2 2]
%     arrdelay = zeros(1, arrNum); % 这一项相当于常数项 本来考虑的是天线间时延
    
    
    delayTx = Tc * numTx; % 发射天线维度积累的时延
    arrdx = cfgOut.arrdx; % 归一化阵元间距
    arrdy = cfgOut.arrdy; 
    
    c = physconst('LightSpeed'); % 光速
    lambda = c / fc; % 波长
    Ts = 1 / fs; % ADC采样时间
    
    % 信号模型
    signal = zeros(ADCNum, TotalChirpNum, arrNum); % 信号
    
    noiseI = normrnd(0, 1, ADCNum, TotalChirpNum, arrNum); % 正态噪声
    noiseQ = normrnd(0, 1, ADCNum, TotalChirpNum, arrNum); % 正态噪声
%     noiseI = exprnd(1, ADCNum, TotalChirpNum, arrNum); % 指数分布噪声
%     noiseQ = exprnd(1, ADCNum, TotalChirpNum, arrNum); % 指数分布噪声
    for tarId = 1 : tarNum
        disp(strcat(['正在获取目标',num2str(tarId),'的数据']));
        % 目标参数
        targetR = tarOut.range(tarId);    % 距离 m
        targetV = tarOut.velocity(tarId); % 速度 m/s
        targetAzi = tarOut.Azi(tarId);    % 目标方位角 度
        targetEle = tarOut.Ele(tarId);    % 目标俯仰角 度
        targetRCS = tarOut.RCS(tarId);    % 目标RCS值
        targetGt  = tarOut.Gt(tarId);     % 发射增益 
        targetGr  = tarOut.Gr(tarId);     % 接收增益 
        [tarSNR] = CalculateSNR(targetR, targetRCS, targetGt, targetGr, lambda, Pt, Fn, Ls, validB); % 信噪比
        A = sqrt(2 * db2pow(tarSNR)); % 信号幅度
        
        targPhi0 = rand * 2 * pi; % 产生随机相位[0 2*pi]

        tempSigI = zeros(ADCNum, TotalChirpNum, arrNum); % I路缓存信号
        tempSigQ = zeros(ADCNum, TotalChirpNum, arrNum); % Q路缓存信号
        tempSig = zeros(ADCNum, TotalChirpNum, arrNum); % I路缓存信号
        
        for channelId = 1 : arrNum
            for chirpId = 1 : TotalChirpNum
                for adcId = 1 : ADCNum
                    % 动目标时延
                    tarDelay = 2 * targetR / c + 2 * targetV * ((chirpId - 1) * delayTx + arrdelay(channelId) + adcId * Ts) / c; 
                    
                    % 控制来波方向
                    % 初阶版本的目标相位表达
%                     tarPhi = targPhi0 + 2 * pi * (arr(1, channelId) * sind(targetAzi) * arrdx + ...
%                             arr(2, channelId) * sind(targetEle) * arrdy) + deg2rad(antennaPhase(channelId)); % 考虑阵元初始、随机相位和目标波达角

                    tarPhi = targPhi0 + 2 * pi * (arr(channelId) * sind(targetAzi) * arrdx + ...
                            arr(arrNum+channelId) * sind(targetEle) * arrdy) + deg2rad(antennaPhase(channelId)); % 这里考虑了不同发射顺序天线相位

                    % 中频信号：exp[1j * 2 *pi * fc * tau + 2 * pi * S * t * tau - pi * tau * tau]  t的范围为0-Tc tau的范围为t+nTc
                    tempSigI(adcId, chirpId, channelId) = A * cos(2 * pi * (fc * tarDelay + Slope * tarDelay * adcId * Ts - Slope * tarDelay * tarDelay / 2) + tarPhi);
                    tempSigQ(adcId, chirpId, channelId) = A * sin(2 * pi * (fc * tarDelay + Slope * tarDelay * adcId * Ts - Slope * tarDelay * tarDelay / 2) + tarPhi);
                    tempSig(adcId, chirpId, channelId) = tempSigI(adcId, chirpId, channelId) + 1j * tempSigQ(adcId, chirpId, channelId);
                end
            end
        end
        signal = signal + tempSig; % 多目标信号相加
    end
    signal = signal - (noiseI + noiseQ); % 考虑噪声
%     signal = reshape(signal, size(signal,1), size(signal,2), numRx, numTx);
end