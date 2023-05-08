function [signal] = GenerateSigI(tarOut, cfgOut)
    %% 本文件用于单路ADC信号
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
    arr = cfgOut.array; % 虚拟阵元排列
    antennaPhase = cfgOut.antennaPhase; % 天线相位
    arrNum = size(arr, 2); % 虚拟阵列数目
    arrdelay = zeros(1, arrNum) + Tc * reshape(repmat([0:numTx-1], numRx,1), 1, arrNum) * Tc; 
    % [0 0 0 0 1 1 1 1 2 2 2 2]
    
    delayTx = Tc * numTx; % 发射天线维度积累的时延
    arrdx = cfgOut.arrdx; % 归一化阵元间距
    arrdy = cfgOut.arrdy; 
    
    c = physconst('LightSpeed'); % 光速
    lambda = c / fc; % 波长
    Ts = 1 / fs; % ADC采样时间
    
    % 信号模型
    signal = zeros(ADCNum, TotalChirpNum, arrNum); % 信号
    noise = normrnd(0, 1, ADCNum, TotalChirpNum, arrNum); % 正态噪声
    
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
        [tarSNR] = CalculateSNR(targetR, targetRCS, targetGt, targetGr, lambda, Pt, Fn, Ls, fs); % 信噪比
        A = sqrt(2 * db2pow(tarSNR)); % 信号幅度
        
        targPhi0 = rand * 2 * pi; % 产生随机相位[0 2*pi]
        tempSig = zeros(ADCNum, TotalChirpNum, arrNum); % 缓存信号
        
        for channelId = 1 : arrNum
            for chirpId = 1 : TotalChirpNum
                for adcId = 1 : ADCNum
                    tarDelay = 2 * targetR / c + 2 * targetV * ((chirpId - 1) * delayTx + arrdelay(channelId) + adcId * Ts) / c; % 动目标时延
                    tarPhi = targPhi0 + 2 * pi * (arr(1, channelId) * sind(targetAzi) * arrdx + ...
                             arr(2, channelId) * sind(targetEle) * arrdy) + deg2rad(antennaPhase(channelId)); % 考虑阵元初始、随机相位和目标波达角
                         
                    % 中频信号：exp[1j * 2 *pi * fc * tau + 2 * pi * S * t * tau - pi * tau * tau]  t的范围为0-Tc tau的范围为t+nTc
                    tempSig(adcId, chirpId, channelId) = A * cos(2 * pi * (fc * tarDelay + Slope * tarDelay * adcId * Ts - Slope * tarDelay * tarDelay / 2) + tarPhi);

                end
            end
        end
        signal = signal + tempSig; % 多目标信号相加
    end
    signal = signal + noise; % 考虑噪声
    
end