function RAM = GenerateRAM(rangeFFTOut, cfgOut, cfgDOA)
    %% 本文件用于实现距离-方位图生成
    %% 这里丢弃掉了多普勒维信息来存储快拍
    %% By Xuliang,20230421
    
    rangeNum = size(rangeFFTOut, 1); % 获取距离单元数目
    thetaGrids = cfgDOA.thetaGrids; % 网格
    virtual_array = cfgOut.virtual_array; % 虚拟阵列
    
    RAM = zeros(rangeNum, length(thetaGrids)); % 初始化RAM
    for rangeIdx = 1 : rangeNum
       arrData = squeeze(rangeFFTOut(rangeIdx, :, :)); % 阵列数据 ChirpNum * arrNum
       sig = reshape(arrData, size(arrData,1), cfgOut.numRx, cfgOut.numTx); % ChirpNum * RXNum * TXNum 
       tmpSig = zeros(size(arrData,1), max(virtual_array.azi_arr)+1, max(virtual_array.ele_arr)+1); % 初始化信号子空间
       for trx_id = 1 : size(cfgOut.sigIdx,2)
           tmpSig(:, cfgOut.sigSpaceIdx(1, trx_id), cfgOut.sigSpaceIdx(2,trx_id)) = sig(:, cfgOut.sigIdx(1,trx_id), cfgOut.sigIdx(2,trx_id)); % 重排后的信号空间
       end
       aziData = squeeze(tmpSig(:, :, 1)); % 仅提取方位通道
       [doaOut] = azimuthDOA(aziData.', cfgDOA);
       RAM(rangeIdx, :) = doaOut.spectrum; % 存储每个距离单元对应的空间谱

    end
end
