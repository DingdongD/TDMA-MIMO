function [cfarOut] = CFAR_2D(RDM, Pfa, TestCells, GuardCells)
    % rd_map：输入数据
    % Pfa：虚警概率
    % TestCells：测试单元 [4 4]
    % GuardCells：保护单元 [4 4]
    
    % CFAR检测器
    detector = phased.CFARDetector2D('TrainingBandSize',TestCells, ...
    'ThresholdFactor','Auto','GuardBandSize',GuardCells, ...
    'ProbabilityFalseAlarm',Pfa,'Method','SOCA','ThresholdOutputPort',true, 'NoisePowerOutputPort',true);

    N_x = size(RDM,1); % 距离单元
    N_y = size(RDM,2); % 多普勒单元
    % 获取二维滑动窗口的size
    Ngc = detector.GuardBandSize(2); % 保护窗列向
    Ngr = detector.GuardBandSize(1); % 保护窗行向
    Ntc = detector.TrainingBandSize(2); % 参考窗行向
    Ntr = detector.TrainingBandSize(1); % 参考窗行向
    cutidx = [];
    colstart = Ntc + Ngc + 1; % 列窗首
    colend = N_y + ( Ntc + Ngc); % 列窗尾
    rowstart = Ntr + Ngr + 1; % 行窗首
    rowend = N_x + ( Ntr + Ngr); % 行窗尾
    for m = colstart:colend
        for n = rowstart:rowend
            cutidx = [cutidx,[n;m]]; % 完整二维窗
        end
    end
    
    ncutcells = size(cutidx,2); % 获取窗口遍历区间
    
    rd_map_padding = repmat(RDM, 3, 3); % 对RDM补零
    chosen_rd_map = rd_map_padding(N_x+1-Ntr-Ngr:2*N_x+Ntr+Ngr,N_y+1-Ntc-Ngc:2*N_y+Ntc+Ngc);
    
    [dets, ~, noise] = detector(chosen_rd_map,cutidx); % 输出CFAR结果
    
    cfar_out = zeros(size(chosen_rd_map)); % 检测点输出
    noise_out = zeros(size(chosen_rd_map)); % 噪声
    snr_out = zeros(size(chosen_rd_map)); % 信噪比
    for k = 1:size(dets,1)
        if dets(k) == 1
            cfar_out(cutidx(1,k),cutidx(2,k)) = dets(k); 
            noise_out(cutidx(1,k),cutidx(2,k)) = noise(k);
            snr_out(cutidx(1,k),cutidx(2,k)) = chosen_rd_map(cutidx(1,k),cutidx(2,k));
        end
    end
    
    cfarOut = {};
    cfarOut.cfarMap = cfar_out(Ntr+Ngr+1:Ntr+Ngr+N_x,Ntc+Ngc+1:Ntc+Ngc+N_y);
    cfarOut.snrOut = snr_out(Ntr+Ngr+1:Ntr+Ngr+N_x,Ntc+Ngc+1:Ntc+Ngc+N_y);
    cfarOut.noiseOut = noise_out(Ntr+Ngr+1:Ntr+Ngr+N_x,Ntc+Ngc+1:Ntc+Ngc+N_y);
    cfarOut.snrOut = cfarOut.snrOut ./ (eps + cfarOut.noiseOut);
end