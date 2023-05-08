function [com_dopplerFFTOut] = compensate_doppler(doaInput, cfgOut, dopplerIdx, speedVal, rangeVal)
    %% 本文件用于生成由TDMA-MIMO导致的多普勒补偿
    %% By Xuliang, 20230418
    tarNum = length(dopplerIdx); % 目标数目
    numTx = cfgOut.numTx;
    numRx = cfgOut.numRx;
    ChirpNum = cfgOut.ChirpNum;
    angleFFT_size = 128;
%     dopplerIdx = rd_peak_list(2, :);

    if cfgOut.applyVmaxExtend 
        sig_bin = [];
        info_overlaped_diff1tx = cfgOut.virtual_array.info_overlaped_diff1tx; 
        if mod(numTx, 2) == 1 % 天线数目为奇数
            tmpDopplerIdx = reshape(dopplerIdx, [], 1); % 临时多普勒索引
            dopplerInd_unwrap = tmpDopplerIdx + ((1:numTx) - ceil(numTx / 2)) * ChirpNum; % tarNum * TXNum
        else % 天线数目为偶数
            tmpDopplerIdx = reshape(dopplerIdx, [], 1); % 临时多普勒索引
            if speedVal > 0
                dopplerInd_unwrap = tmpDopplerIdx + ((1:numTx)- (numTx / 2 + 1)) * ChirpNum; % tarNum * TXNum
            else
                dopplerInd_unwrap = tmpDopplerIdx + ((1:numTx)- numTx / 2) * ChirpNum; % tarNum * TXNum
            end
        end
        
        sig_bin_org = doaInput; % TARNUM * RXNUM * TXNUM
        deltaPhi = 2 * pi * (dopplerInd_unwrap - ChirpNum / 2) / (numTx * ChirpNum); % 多普勒相位修正 tarNum * TXNUM 
        deltaPhis(:, 1, 1, :) = deltaPhi; % tarNum * 1 * 1 * TXNUM
        tmpTX = (0 : numTx - 1);
        tmpTXs(1,1,:,1) = tmpTX; % tarNum * 1 * TXNUM
        correct_martrix = exp(-1j * tmpTXs .* deltaPhis); % tarNum * 1 * TXNUM * TXNUM
        sig_bin = sig_bin_org .* correct_martrix; % tarNum * RXNUM * TXNUM * TXNUM

        % 使用重叠天线来做最大速度的解包 
        nore_tx = [];
        nore_rx = [];
        re_tx = [];
        re_rx = [];
        for iid = 1 : length(info_overlaped_diff1tx)
            nore_rx_set = info_overlaped_diff1tx(:,3);
            nore_tx_set = info_overlaped_diff1tx(:,4);
            re_rx_set = info_overlaped_diff1tx(:,7);
            re_tx_set = info_overlaped_diff1tx(:,8);
            nore_rx = [nore_rx, find(cfgOut.PosRX_BOARD_ID == nore_rx_set(iid))];
            nore_tx = [nore_tx, find(cfgOut.PosTX_Trans_ID == nore_tx_set(iid))];
            re_rx = [re_rx, find(cfgOut.PosRX_BOARD_ID == re_rx_set(iid))];
            re_tx = [re_tx, find(cfgOut.PosTX_Trans_ID == re_tx_set(iid))];
        end
        index_overlaped_diff1tx = cat(2, nore_rx.', nore_tx.', re_rx.', re_tx.'); % 32 * 4
        
        % 获取未冗余阵元和冗余阵元对应的信号
        sig_overlap1 = zeros(tarNum, length(index_overlaped_diff1tx));
        sig_overlap2 = zeros(tarNum, length(index_overlaped_diff1tx));
        for iid = 1 : length(index_overlaped_diff1tx)
            tmppos1 = index_overlaped_diff1tx(:, 1);
            tmppos2 = index_overlaped_diff1tx(:, 2);
            tmppos3 = index_overlaped_diff1tx(:, 3);
            tmppos4 = index_overlaped_diff1tx(:, 4);
            sig_overlap1(:,iid) = sig_bin_org(:, tmppos1(iid), tmppos2(iid));
            sig_overlap2(:, iid) = sig_bin_org(:, tmppos1(iid), tmppos2(iid));
        end
        sig_overlap = cat(3, sig_overlap1, sig_overlap2); % tarNum * 32 * 2
        
        angle_sum_test = zeros(size(sig_overlap,1), size(sig_overlap,2), size(deltaPhi,2)); % 检查每个假设的重叠天线对相位差
        
        for sig_id = 1 : size(angle_sum_test,2)
            deltaPhiss(:, 1, :) = deltaPhi; % tarNum * 1 * TXnum
            tmpSigs = sig_overlap(:, 1:sig_id, 2); % tarNum * sig_id * 1
            signal2 = tmpSigs .* exp(-1j * deltaPhiss); % tarNum * 32 * TXnum
            angle_sum_test(:, sig_id, :) = angle(sum(sig_overlap(:,1:sig_id,1) .* conj(signal2),2)); % tarNum * 32 * TXnum
        end
        [~, doppler_unwrap_integ_overlap_index] = min(abs(angle_sum_test),[],3); % tarNum * 1 * 32
        doppler_unwrap_integ_overlap_index = squeeze(doppler_unwrap_integ_overlap_index); % 选择相位差最小的假设来估计解系数
        
        doppler_unwrap_integ_index = zeros(size(doppler_unwrap_integ_overlap_index, 1), 1);
        for i = 1:size(doppler_unwrap_integ_overlap_index, 1)
            doppler_unwrap_integ_index(i) = mode(doppler_unwrap_integ_overlap_index(i, :));
        end
        
        correct_sig = []; % tarNum * numRX * numTX
        for i = 1 : length(doppler_unwrap_integ_index)
            correct_sig(i, :, :) = squeeze(sig_bin(i, :, :, doppler_unwrap_integ_index(i))); 
        end
        
        index_noapply = find(rangeVal <= cfgOut.min_dis_apply_vmax_extend); % 未超出最小距离约束
        delta_phi_noapply = reshape(2 * pi * (dopplerIdx - ChirpNum / 2) / (numTx * ChirpNum),[],1);
        sig_bin_noapply = doaInput .* (reshape(exp(-1j * delta_phi_noapply * tmpTX), size(delta_phi_noapply, 1), 1, size(tmpTX,2))); % 单一相位约束
        correct_sig(index_noapply, :, :) = sig_bin_noapply(index_noapply, :, :);
        
        com_dopplerFFTOut = correct_sig;
        
%         % 下面这部分暂时没有用到 先不做考虑 后面有时间再维护
%         % 角度维FFT信噪比求解
%         noredundant_arr = cfgOut.virtual_array.noredundant_aziarr; % 获取方位维度天线
%         noredundant_rows0 = [];
%         noredundant_rows1 = [];
%         rx_pos_set = noredundant_arr(:,3); % 接收天线位置集合
%         tx_pos_set = noredundant_arr(:,4); % 发射天线位置集合
%         for rx_id = 1 :length(virtual_array.noredundant_arr(rx_pos_set))
%             rx_arr = rx_pos_set(rx_id);
%             noredundant_rows0 = [noredundant_rows0, find(cfgOut.PosRX_BOARD_ID == rx_arr)]; % 确定接收天线的序号
%         end
%         for tx_id = 1 :length(virtual_array.noredundant_arr(tx_pos_set))
%             tx_arr = tx_pos_set(tx_id);
%             noredundant_rows1 = [noredundant_rows1, find(cfgOut.PosTX_Trans_ID == tx_arr)]; % 确定发射天线的顺序序号 不是板子上排列
%         end
%         noredundant_rows = [noredundant_rows0;noredundant_rows1]; 
%         
%         sig_bin_row1 = zeros(tarNum, max(virtual_array.azi_arr)+1, max(virtual_array.ele_arr)+1, numTx); 
%                 
%         for trx_id = 1 : size(noredundant_rows1,1) % 86个阵列
%             sig_bin_row1(:, noredundant_rows(1,trx_id), noredundant_rows(2,trx_id), :) = sig_bin(:, noredundant_rows(1,trx_id), noredundant_rows(2,trx_id), :); % 重排后的信号空间 1*12
%         end
%         sig_bin_row1 = squeeze(sig_bin_row1(:,:,1,:)); % tarNum * 86 * txNum
%         sig_bin_row1_fft = fftshift(fft(sig_bin_row1, angleFFT_size, 2), 2); % 方位维FFT
%         
%         angle_bin_skip_left = 4;
%         angle_bin_skip_right = 4;
%         sig_bin_row1_fft_cut = abs(sig_bin_row1_fft(:, angle_bin_skip_left + 1 : angleFFT_size - angle_bin_skip_right, :));
%         [max_val, max_idx] = max(max(abs(sig_bin_row1_fft_cut), [], 2), [], 3);
        
    else
        sig_bin_org = doaInput; % TARNUM * RXNUM * TXNUM
        deltaPhi = 2 * pi * (dopplerIdx - ChirpNum / 2) / (numTx * ChirpNum); % 多普勒相位修正 tarNum * 1
        deltaPhi = deltaPhi.';
        tmpTX = (0 : numTx - 1); % 1 * TXNUM
        correct_martrix = exp(-1j * deltaPhi * tmpTX ); % TARNUM * TXNUM 
        correct_martrixs(:, 1, :) = correct_martrix;
        com_dopplerFFTOut = sig_bin_org .* correct_martrixs; % TARNUM * RXNUM * TXNUM
    end
end