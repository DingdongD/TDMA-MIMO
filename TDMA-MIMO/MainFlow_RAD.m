%% 本文件用于测试Range-Azimuth-Doppler流程下的信号处理栈
%% By Xuliang,20230414
clc;clear;close all;

dataPath = './dataset/adc_raw_dat.mat'; % 保存文件路径
ShowIQ = 1; % 是否显示IQ信号
ShowRange = 1; % 是否显示RangePorfile
ShowRA = 1; % 是否显示RAM
ShowRD = 1; % 是否显示RDM
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


for frame_id = 1 : Frame
    adcData = squeeze(RawData(:, :, frame_id, :));
    
    %% 距离维FFT
    disp(strcat(['=====','Range-Profile生成','====='])); 
    tic
    fftOut1 = rangeFFT(adcData, IQFlag);
    toc
    rangeFFTOut = fftOut1.rangeFFT;
    
    if ShowRange
        figure(2);
        set(gcf,'unit','centimeters','position',[20,12,10,10])
        plot(rangeIndex, db(rangeFFTOut(:,1,1))); 
        xlabel('Range(m)');ylabel('Amplitude(dB)'); 
        title(strcat(['第',num2str(frame_id),'帧-目标距离分布']));grid minor;  
        pause(0.1);
    end
    
    %% 距离-方位图生成
    cfgDOA.FFTNum = 180; % FFT点数
    cfgDOA.AziMethod = 'MUSIC'; % 由于本套架构目的是利用多快拍信息 仅支持超分算法
    cfgDOA.AzisigNum = 1; % 约束每个CELL上的信源数目
    cfgDOA.thetaGrids = linspace(-90, 90, cfgDOA.FFTNum); % 空间网格
    
    disp(strcat(['=====','Range-Azimuth-Map生成','====='])); 
    tic
    RAM = GenerateRAM(rangeFFTOut, cfgOut, cfgDOA);
    toc
    if ShowRA
        figure(3);
        set(gcf,'unit','centimeters','position',[30,12,10,10])
        [theta,rho] = meshgrid(cfgDOA.thetaGrids, rangeIndex); % 网格化
        xaxis = rho .* cosd(theta); % 横坐标
        yaxis = rho .* sind(theta); % 纵坐标
        surf(yaxis,xaxis,db(abs(RAM)),'EdgeColor','none'); 
        view(2);colormap('jet');
        xlabel('Range(m)','fontsize',15,'fontname','Times New Roman');ylabel('Range(m)','fontsize',15,'fontname','Times New Roman');grid on; axis xy
        title(strcat(['第',num2str(frame_id),'帧-距离方位谱图']));grid minor;  
        set(gca,'GridLineStyle','- -');
        set(gca,'GridAlpha',0.2);
        set(gca,'LineWidth',1.5);
        set(gca,'xminortick','on'); 
        set(gca,'ygrid','on','GridColor',[0 0 0]);
        colorbar;
    end
    
    %% CFAR检测
    disp(strcat(['=====','恒虚警检测','====='])); 
    Pfa = 1e-2; % 虚警概率
    TestCells = [8, 8]; % 参考从
    GuardCells = [4, 4]; % 保护窗
    
    tic
    [cfarOut] = CFAR_2D(abs(RAM).^2, Pfa, TestCells, GuardCells);
    toc
    cfarMap = cfarOut.cfarMap; % 检测点输出
    noiseOut = cfarOut.noiseOut; % 噪声输出
    snrOut = cfarOut.snrOut; % 信噪比输出 RAM的信噪比不能表示真实信噪比(部分空谱方法中)
    
    if ShowCFAR
        figure(4);
        set(gcf,'unit','centimeters','position',[30,0,10,10]);
        imagesc(cfgDOA.thetaGrids, rangeIndex, cfarMap);
        ylabel('Range(m)');xlabel('Angle(deg)'); colormap('jet');
        title(strcat(['第',num2str(frame_id),'帧目标-CFAR检测结果']));
        grid minor;axis xy;
        pause(0.1);
    end
    
    %% 峰值聚合-获取点目标
    disp(strcat(['=====','峰值聚合','====='])); 
    [range_idx, doppler_idx] = find(cfarMap);
    cfar_out_idx = [range_idx doppler_idx]; % 获取CFAR输出结果的行列索引
    tic
    [rd_peak_list, rd_peak] = peakFocus(db(RAM), cfar_out_idx);
    toc
    peakMap = zeros(size(cfarMap)); % 峰值聚合结果矩阵
    for peak_idx = 1 :size(rd_peak_list, 2)
        peakMap(rd_peak_list(1,peak_idx), rd_peak_list(2,peak_idx)) = 1;
    end
    
    if ShowPeak
        figure(5);
        set(gcf,'unit','centimeters','position',[40,12,10,10])
        imagesc(cfgDOA.thetaGrids, rangeIndex, peakMap);
        ylabel('Range(m)');xlabel('Angle(deg)'); colormap('jet');
        title(strcat(['第',num2str(frame_id),'帧目标-峰值聚合结果']));
        grid minor;axis xy;xlim([-60 60]);
        pause(0.1);
    end
    
    % 这里需要考虑CFAR检测的目标点[存储检测到的目标距离索引、方位索引]
    targetPerFrame = {}; 
    targetPerFrame.rangeSet = []; % 存储距离
    targetPerFrame.velocitySet = []; % 存储速度[为下面初始化]
    targetPerFrame.azimuthSet = []; % 存储方位
    
    [cfar_rid, cfar_cid] = find(cfarMap); % rid表示距离单元索引 cid表示角度单元索引
    targetPerFrame.rangeSet = [targetPerFrame.rangeSet, rangeIndex(cfar_rid)];
    targetPerFrame.azimuthSet = [targetPerFrame.azimuthSet, cfgDOA.thetaGrids(cfar_cid)];
    
    %% 多普勒FFT估计[和传统流程不一样，这个地方需要使用缓存数据rangeFFT]
    fftOut2 = dopplerFFT(rangeFFTOut); % 缓存数据大小为256*128*12
    dopplerFFTOut = fftOut2.dopplerFFT; % 获取多普勒维FFT数据 
    
    %% 非相干积累
    disp(strcat(['=====','非相干积累','====='])); 
    RDM = dopplerFFTOut;
    tic
    [accumulateRD] = incoherent_accumulation(RDM);
    toc
    
    if ShowRD
        figure(6);
        set(gcf,'unit','centimeters','position',[20,0,10,10])
        imagesc(rangeIndex, velocityIndex, db(dopplerFFTOut(:,:,1).'));
        xlabel('Range(m)');ylabel('Velocity(m/s)'); colormap('jet');
        title(strcat(['第',num2str(frame_id),'帧目标-距离多普勒分布']));
        grid minor; axis xy;
        pause(0.1);
    end
    
    for tarRangeIdx = 1 : length(targetPerFrame.rangeSet) % 遍历检测到的距离单元 搜索谱峰估计速度
        [velVal, velIdx] = findpeaks(accumulateRD(cfar_rid(tarRangeIdx),:)); % 搜索速度谱峰
        [sortVal, sortIdx] = sort(velVal); % 对峰值进行排序
        Idx = velIdx(sortIdx(end)); % 选择最大值
        targetPerFrame.velocitySet = [targetPerFrame.velocitySet, velocityIndex(Idx)];
    end
    
end