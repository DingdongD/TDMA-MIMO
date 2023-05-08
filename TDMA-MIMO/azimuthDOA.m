function [doaOut] = azimuthDOA(arrData, cfgDOA)
    %% 本文件用于方位维度的DOA/AOA估计
    %% By Xuliang, 20230417
    %% arrData: 阵列数据 传入方位维阵列*1
    %% cfgDOA: 参数配置

    addpath('./DOA_FUNC'); % 加入DOA函数路径
    doaMethod = cfgDOA.AziMethod; % 读取doa估计方法 每种方法返回空间谱
    doaOut = {};
    if strcmp(doaMethod,'FFT')
        P = cfgDOA.AzisigNum; % P 为信源数目
        [Pout] = DOA_FFT(arrData, cfgDOA);
        if sum(abs(Pout)) == 0
            doaOut.spectrum = Pout; % 空间谱
        else
            [peakVal, peakIdx] = findpeaks(abs(Pout)); % 搜索空间谱峰
            [sortVal, sortIdx] = sort(peakVal); % 对峰值进行排序 级联雷达中存在3列元素全为0 这种情况不存在峰值
            ampVal = sortVal(end-P+1:end); % 选择最大值
            Idx = peakIdx(sortIdx(end)); % 选择最大值
            
            freqGrids = linspace(-pi, pi, cfgDOA.FFTNum); % 频域网格
            freq = freqGrids(Idx);
            angleVal = asind(freq / 2 / pi * 2); % 角度值
            
            doaOut.peakVal = ampVal; % 峰值
            doaOut.angleVal = angleVal; % 角度值
            doaOut.angleIdx = Idx; % 谱峰索引
            doaOut.spectrum = Pout; % 空间谱
        end
        
    elseif strcmp(doaMethod,'MUSIC')
        P = cfgDOA.AzisigNum; % P 为信源数目
        thetaGrids = cfgDOA.thetaGrids; % 网格划分
        
        [Pout] = DOA_MUSIC(arrData, P, thetaGrids); % MUSIC的arrData数据维度可以是M*snap M阵元 snap快拍
        
        if sum(abs(Pout)) == 0
            doaOut.spectrum = Pout; % 空间谱
        else
            [peakVal, peakIdx] = findpeaks(abs(Pout)); % 搜索空间谱峰
            [sortVal, sortIdx] = sort(peakVal); % 对峰值进行排序

            ampVal = sortVal(end-P+1:end); % 选择最大值
            Idx = peakIdx(sortIdx(end-P+1:end)); % 选择最大值
            angleVal = thetaGrids(Idx); % 角度值
            
            doaOut.peakVal = ampVal; % 峰值
            doaOut.angleVal = angleVal; % 角度值
            doaOut.angleIdx = Idx; % 谱峰索引
            doaOut.spectrum = Pout; % 空间谱
        end
    end
    
end