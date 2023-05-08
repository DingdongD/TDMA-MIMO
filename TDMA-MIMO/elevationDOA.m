function [doaOut] = elevationDOA(arrData, cfgDOA)
    %% 本文件用于俯仰维度的DOA/AOA估计
    %% By Xuliang, 20230417
    %% arrData: 阵列数据
    %% cfgDOA: 参数配置

    addpath('./DOA_FUNC'); % 加入DOA函数路径
    
    doaMethod = cfgDOA.EleMethod; % 读取doa估计方法 每种方法返回空间谱
    
    if strcmp(doaMethod,'FFT')
        P = cfgDOA.ElesigNum; % P 为信源数目
        [Pout] = DOA_FFT(arrData, cfgDOA);
        [peakVal, peakIdx] = findpeaks(abs(Pout)); % 搜索空间谱峰
        [sortVal, sortIdx] = sort(peakVal); % 对峰值进行排序
        
        ampVal = sortVal(end-P+1:end); % 选择最大值
        Idx = peakIdx(sortIdx(end)); % 选择最大值
        
        freqGrids = linspace(-pi, pi, cfgDOA.FFTNum); % 频域网格
        freq = freqGrids(Idx);
        angleVal = asind(freq / 2 / pi * 2); % 角度值
        
    elseif strcmp(doaMethod,'MUSIC')
        P = cfgDOA.ElesigNum; % P 为信源数目
        thetaGrids = cfgDOA.thetaGrids; % 网格划分
        [Pout] = DOA_MUSIC(arrData, P, thetaGrids); % MUSIC的arrData数据维度可以是M*snap M阵元 snap快拍
        
        [peakVal, peakIdx] = findpeaks(abs(Pout)); % 搜索空间谱峰
        [sortVal, sortIdx] = sort(peakVal); % 对峰值进行排序
        
        ampVal = sortVal(end-P+1:end); % 选择最大值
        Idx = peakIdx(sortIdx(end-P+1:end)); % 选择最大值
        angleVal = thetaGrids(Idx); % 角度值
        
    end
    
    
    doaOut = {};
    doaOut.peakVal = ampVal; % 峰值
    doaOut.angleVal = angleVal; % 角度值
    doaOut.angleIdx = Idx; % 谱峰索引
    doaOut.spectrum = Pout; % 复空间谱

end