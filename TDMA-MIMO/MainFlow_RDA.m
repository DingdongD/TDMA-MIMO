%% 本文件用于处理TDMA-MIMO雷达信号
%% By Xuliang,20230411
clc;clear;close all;

dataPath = './dataset/adc_raw_dat.mat'; % 保存文件路径
ShowIQ = 1; % 是否显示IQ信号
ShowRange = 1; % 是否显示RangePorfile
ShowRD = 1; % 是否显示RD
ShowCFAR = 1; % 是否显示CFAR结果
ShowPeak = 1; % 是否显示聚合峰值结果
IQFlag = 1; % 是否选择IQ路信号 0-单路 1-双路
saveFlag = 0; % 是否保存文件

%% 目标和系统参数更新
disp(strcat(['=====','刷新目标和系统参数','====='])); % 单板雷达模式
tarOut = ConfigureTarget;     % 生成目标信息
cfgOut = ConfigureParameter;  % 生成毫米波雷达系统参数

%% 雷达采集数据
disp(strcat(['=====','雷达进入工作状态','=====']));     
[RawData] = GenerateAdcData(tarOut, cfgOut, IQFlag, saveFlag, dataPath); % 开始模拟采数据

%% 物理参数更新
c = physconst('LightSpeed'); % 光速
fc = cfgOut.fc; % 载频 Hz
lambda = c / fc; % 波长

ADCNum = cfgOut.ADCNum; % ADC采样数目
ChirpNum = cfgOut.ChirpNum; % 每帧发射Chirp数目

numTx = cfgOut.numTx; % 发射天线数目
numRx = cfgOut.numRx; % 接收天线数目
% arr = cfgOut.array; % 阵元排列[初阶版本的阵列]
virtual_array = cfgOut.virtual_array; % 虚拟阵列struct
arr = virtual_array.virtual_arr;

arrNum = numTx * numRx; % 阵元数目
arrDx = cfgOut.arrdx; % 方位向阵元间距
arrDy = cfgOut.arrdx; % 俯仰向阵元间距

validB = cfgOut.validB; % 有效带宽
range_res = c / (2 * validB); % 距离分辨率

TF = cfgOut.Tc * (numTx * ChirpNum);
doppler_res = lambda / (2 * TF); % 多普勒分辨率
Frame = cfgOut.Frame; % 帧数

% 距离索引和速度索引
if IQFlag
    velocityIndex = [-ChirpNum / 2 : 1 : ChirpNum / 2 - 1] * doppler_res;
    rangeIndex = (0 : ADCNum - 1) * range_res;
else
    velocityIndex = [-ChirpNum / 2 : 1 : ChirpNum / 2 - 1] * doppler_res;
    rangeIndex = (-ADCNum / 2 : 1 : ADCNum / 2 - 1) * range_res;
    rangeIndex = rangeIndex(end/2+1:end);  % 取后面一半
end

%% IQ平衡分析
if ShowIQ
    if IQFlag % IQ两路信号
        figure(1);
        set(gcf,'unit','centimeters','position',[10,12,10,10])
        plot(1:ADCNum, real(RawData(:,1,1))); hold on;
        plot(1:ADCNum, imag(RawData(:,1,1))); hold off;
        xlabel('ADC Points');ylabel('Amplitude');title('IQ-signal Analysis');
        legend('I-Chain','Q-Chain');grid minor;
    else
        figure(1);
        set(gcf,'unit','centimeters','position',[10,12,10,10])
        plot(1:ADCNum, (RawData(:,1,1)));
        xlabel('ADC Points');ylabel('Amplitude');title('IQ-signal Analysis');grid minor;
    end
end
RawData = reshape(RawData, ADCNum, ChirpNum, [], arrNum);
disp(strcat(['雷达信号的维度为：',num2str(size(RawData))]));

for frame_id = 1 : Frame % Frame
    adcData = squeeze(RawData(:, :, frame_id, :));
    
    %% 距离维FFT和多普勒FFT
    disp(strcat(['=====','RD-MAP生成','====='])); 
    tic
    fftOut = rdFFT(adcData, IQFlag);
    toc
    rangeFFTOut = fftOut.rangeFFT;
    dopplerFFTOut = fftOut.dopplerFFT;
    if ShowRange
        figure(2);
        set(gcf,'unit','centimeters','position',[10,0,10,10])
        plot(rangeIndex, db(rangeFFTOut(:,1,1))); 
        xlabel('Range(m)');ylabel('Amplitude(dB)'); 
        title(strcat(['第',num2str(frame_id),'帧-目标距离分布']));grid minor;  
        pause(0.1);
    end
    if ShowRD
        figure(3);
        set(gcf,'unit','centimeters','position',[20,12,10,10])
        imagesc(rangeIndex, velocityIndex, db(dopplerFFTOut(:,:,1).'));
        xlabel('Range(m)');ylabel('Velocity(m/s)'); colormap('jet');
        title(strcat(['第',num2str(frame_id),'帧目标-距离多普勒分布']));
        grid minor; axis xy;
        pause(0.1);
    end

    %% 非相干积累
    disp(strcat(['=====','非相干积累','====='])); 
    RDM = dopplerFFTOut;
    tic
    [accumulateRD] = incoherent_accumulation(RDM);
    toc
    
    %% CFAR检测器
    disp(strcat(['=====','恒虚警检测','====='])); 
    Pfa = 1e-3; % 虚警概率
    TestCells = [8, 8]; % 参考窗
    GuardCells = [2, 2]; % 保护窗
    
    tic
    [cfarOut] = CFAR_2D(accumulateRD, Pfa, TestCells, GuardCells);
    toc
    cfarMap = cfarOut.cfarMap; % 检测点输出
    noiseOut = cfarOut.noiseOut; % 噪声输出
    snrOut = cfarOut.snrOut; % 信噪比输出
    
    if ShowCFAR
        figure(4);
        set(gcf,'unit','centimeters','position',[20,0,10,10])
        imagesc(rangeIndex, velocityIndex, cfarMap.');
        xlabel('Range(m)');ylabel('Velocity(m/s)'); colormap('jet');
        title(strcat(['第',num2str(frame_id),'帧目标-CFAR检测结果']));
        grid minor;axis xy;
        pause(0.1);
    end
    
    %% 峰值聚合-获取点目标
    disp(strcat(['=====','峰值聚合','====='])); 
    [range_idx, doppler_idx] = find(cfarMap);
    cfar_out_idx = [range_idx doppler_idx]; % 获取CFAR输出结果的行列索引
    tic
    [rd_peak_list, rd_peak] = peakFocus(db(accumulateRD), cfar_out_idx);
    toc
    peakMap = zeros(size(cfarMap)); % 峰值聚合结果矩阵
    for peak_idx = 1 :size(rd_peak_list, 2)
        peakMap(rd_peak_list(1,peak_idx), rd_peak_list(2,peak_idx)) = 1;
    end
    
    if ShowPeak
        figure(5);
        set(gcf,'unit','centimeters','position',[30,12,10,10])
        imagesc(rangeIndex, velocityIndex, peakMap.');
        xlabel('Range(m)');ylabel('Velocity(m/s)'); colormap('jet');
        title(strcat(['第',num2str(frame_id),'帧目标-峰值聚合结果']));
        grid minor;axis xy;
        pause(0.1);
    end

    %% DOA/AOA估计
    disp(strcat(['=====','DOA/AOA估计','====='])); 
    cfgDOA.FFTNum = 180; % FFT点数
    % 这里虽然封装了不同的DOA算法 但是需要注意的是 可用的算法有限 在本套代码里 建议使用的是FFT-FFT和FFT-MUSIC等
    % 因为MUSIC-MUSIC的使用会导致在方位维度空间谱估计时破坏相位信号 有损俯仰维度的相位估计
    
    cfgDOA.AziMethod = 'FFT'; % 方位维度DOA估计方法
    cfgDOA.EleMethod = 'MUSIC'; % 俯仰维度DOA估计方法
    
    cfgDOA.thetaGrids = linspace(-90, 90, cfgDOA.FFTNum); % 空间网格
    cfgDOA.AzisigNum = 1; % 约束每个RD-CELL上的信源数目
    cfgDOA.ElesigNum = 1; % 约束每个方位谱峰上的信源数目
    
    aziNum = length(virtual_array.noredundant_aziarr); % 方位天线数目
    
%     Aset = exp(-1j * 2 * pi * arrDx * [0:aziNum]' * sind(cfgDOA.thetaGrids)); % 稀疏字典设计
    targetPerFrame = {}; 
    targetPerFrame.rangeSet = [];
    targetPerFrame.velocitySet = [];
    targetPerFrame.snrSet = [];
    targetPerFrame.azimuthSet = [];
    targetPerFrame.elevationSet = [];
    
    if ~isempty(rd_peak_list) % 非空表示检测到目标
        rangeVal = (rd_peak_list(1, :) - 1) * range_res; % 目标距离
        speedVal = (rd_peak_list(2, :) - ChirpNum / 2 - 1) * doppler_res; % 目标速度
        
        doaInput = zeros(size(rd_peak_list, 2), arrNum);
        for tar_idx = 1 :size(rd_peak_list, 2)
            doaInput(tar_idx, :) = squeeze(dopplerFFTOut(rd_peak_list(1, tar_idx), rd_peak_list(2, tar_idx), :)); % tarNum * arrNum
        end
        doaInput = reshape(doaInput, [], numRx, numTx);
        
        % 方位角估计前需要考虑多普勒补偿
        [com_dopplerFFTOut] = compensate_doppler(doaInput, cfgOut, rd_peak_list(2, :), speedVal, rangeVal); 
        
        tic
        for peak_idx = 1 : size(rd_peak_list, 2) % 遍历检测到的每个目标
            snrVal = mag2db(snrOut(rd_peak_list(1, peak_idx), rd_peak_list(2, peak_idx))); % 信噪比的提升是由于chirpNum*ADCNum的积累
            tarData = squeeze(com_dopplerFFTOut(peak_idx, :,:));

    %         aziarr = unique(arr(1,arr(2,:)==0)); % 初阶版本获取方位维度天线排列
    %         aziArrData = arrData(aziarr+1); % 获取方位维度信号

            % 方位角解析
           sig = tarData;
           sig_space = zeros(max(virtual_array.azi_arr)+1,max(virtual_array.ele_arr)+1); % 初始化信号子空间
           for trx_id = 1 : size(cfgOut.sigIdx,2)
               sig_space(cfgOut.sigSpaceIdx(1, trx_id), cfgOut.sigSpaceIdx(2,trx_id)) = sig(cfgOut.sigIdx(1,trx_id), cfgOut.sigIdx(2,trx_id)); % 重排后的信号空间
           end
           % 至此我们生成的信号子空间维度为 方位虚拟天线数目 * 俯仰虚拟天线数目 

           eleArrData = zeros(cfgDOA.FFTNum, size(sig_space,2)); % 俯仰维度数据
            for ele_idx = 1 : size(sig_space, 2) % 这里采取遍历是为了适配不同的空间谱估计方法
                tmpAziData = sig_space(:, ele_idx);
                [azidoaOut] = azimuthDOA(tmpAziData, cfgDOA); % 提取第一列方位维度天线信息进行方位角估计 
                eleArrData(:, ele_idx) = azidoaOut.spectrum(:); % 复空间谱
            end

            for azi_peak_idx = 1 : length(azidoaOut.angleVal) % 对方位维度检测的谱峰进行检索
                tmpEleData = eleArrData(azidoaOut.angleIdx(azi_peak_idx), :).'; % 获取与方位维目标关联的信号
                [eledoaOut] = elevationDOA(tmpEleData, cfgDOA); % 进行俯仰角估计

                % 关联目标的距离、多普勒、信噪比、方位和俯仰信息
                aziVal = azidoaOut.angleVal; 
                eleVal = eledoaOut.angleVal;
                targetPerFrame.rangeSet = [targetPerFrame.rangeSet, rangeVal(peak_idx)];
                targetPerFrame.velocitySet = [targetPerFrame.velocitySet, speedVal(peak_idx)];
                targetPerFrame.snrSet = [targetPerFrame.snrSet, snrVal];
                targetPerFrame.azimuthSet = [targetPerFrame.azimuthSet,aziVal];
                targetPerFrame.elevationSet = [targetPerFrame.elevationSet,eleVal];
            end
        end
        toc
    
        %% 点云生成 
        disp(strcat(['=====','点云生成','====='])); 
        tic
        pcd_x = targetPerFrame.rangeSet .* cosd(targetPerFrame.elevationSet) .* sind(targetPerFrame.azimuthSet);
        pcd_y = targetPerFrame.rangeSet .* cosd(targetPerFrame.elevationSet) .* cosd(targetPerFrame.azimuthSet);
        pcd_z = targetPerFrame.rangeSet .* sind(targetPerFrame.elevationSet);
        PointSet = [pcd_x.', pcd_y.', pcd_z.', targetPerFrame.velocitySet.', targetPerFrame.snrSet.'];
        toc

        %% 点云聚类
        eps = 1.1; % 邻域半径
        minPointsInCluster = 3; % 簇内最小点数阈值
        xFactor = 1;   % 变大控制距离变大，变小分类距离变小 椭圆
        yFactor = 1;   % 变大控制角度变大，变小分类距离变小 椭圆 
        figure(6);
        set(gcf,'unit','centimeters','position',[30,0,10,10])
        disp(strcat(['=====','点云聚类','====='])); 
        tic
        [sumClu] = dbscanClustering(eps,PointSet,xFactor,yFactor,minPointsInCluster,frame_id); % DBSCAN聚类
        toc
    end

end
