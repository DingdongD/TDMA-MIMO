function [cfgOut] = ConfigureParameter()
    %% 本文件用于仿真配置毫米波雷达系统参数
    %% By Xuliang, 20230411
    cfgOut.Mode = 1; % 仿真模式
    cfgOut.ADCNum = 256; % ADC采样数目
    cfgOut.ChirpNum = 128; % 每帧发射Chirp数目
    cfgOut.Frame = 1; % 总帧数
    cfgOut.applyVmaxExtend = 1; % 允许速度扩展
    if cfgOut.applyVmaxExtend
        cfgOut.min_dis_apply_vmax_extend = 10; % 允许速度扩展的最小距离
    end
    
    if cfgOut.Mode == 1 
        disp(strcat(['=====','已选择AWR1243-单板雷达模式','====='])); % 单板雷达模式3T4R
        cfgOut.numTx = 3;
        cfgOut.numRx = 4;  
        
        cfgOut.PosTX_X = [0 2 4];
        cfgOut.PosTX_Y = [0 1 0]; 
        cfgOut.PosRX_X = [0 1 2 3];
        cfgOut.PosRX_Y = [0 0 0 0];
        cfgOut.PosTX_BOARD_ID = [1 2 3]; % 发射天线板上顺序 后面没用到
        cfgOut.PosTX_Trans_ID = [1 3 2]; % 发射天线发射顺序 
        cfgOut.PosRX_X = [0 1 2 3];
        cfgOut.PosRX_Y = zeros(1,cfgOut.numRx); 
        cfgOut.PosRX_BOARD_ID = [1 2 3 4]; % 接收天线板上顺序
        
    elseif cfgOut.Mode == 2
        disp(strcat(['=====','已选择AWR2243-级联雷达模式','====='])); % 级联雷达模式12T16R
        cfgOut.numTx = 12;
        cfgOut.numRx = 16;
        cfgOut.PosTX_X = [11 10 9 32 28 24 20 16 12 8 4 0];
        cfgOut.PosTX_Y = [6 4 1 0 0 0 0 0 0 0 0 0]; 
        cfgOut.PosTX_BOARD_ID = [12 11 10 3 2 1 9 8 7 6 5 4]; % 发射天线板上顺序 后面没用到
        cfgOut.PosTX_Trans_ID = [12:-1:1]; % 发射天线发射顺序
        
        cfgOut.PosRX_X = [11 12 13 14 50 51 52 53 46 47 48 49 0 1 2 3];
        cfgOut.PosRX_Y = zeros(1,cfgOut.numRx); 
        cfgOut.PosRX_BOARD_ID = [13 14 15 16 1 2 3 4 9 10 11 12 5 6 7 8]; % 接收天线板上顺序
    elseif cfgOut.Mode == 3
        disp(strcat(['=====','已选择特斯拉6T8R-级联雷达模式','=====']));
        cfgOut.numTx = 6;
        cfgOut.numRx = 8;
        cfgOut.PosTX_X = [0 17 34 41 48 55];
        cfgOut.PosTX_Y = [0 0 0 0 4 8 12]; 
        cfgOut.PosTX_BOARD_ID = [1 2 3 4 5 6]; % 发射天线板上顺序 后面没用到
        cfgOut.PosTX_Trans_ID = [1 2 3 4 5 6]; % 发射天线发射顺序
        
        cfgOut.PosRX_X = [0 2 3 5 8 11 14 17];
        cfgOut.PosRX_Y = zeros(1,cfgOut.numRx); 
        cfgOut.PosRX_BOARD_ID = [1 2 3 4 5 6 7 8]; % 接收天线板上顺序
    end
    
    virtual_azi_arr = repmat(cfgOut.PosTX_X(cfgOut.PosTX_Trans_ID), cfgOut.numRx, 1) + repmat(cfgOut.PosRX_X(cfgOut.PosRX_BOARD_ID), cfgOut.numTx, 1).';
    virtual_ele_arr = repmat(cfgOut.PosTX_Y(cfgOut.PosTX_Trans_ID), cfgOut.numRx, 1) + repmat(cfgOut.PosRX_Y(cfgOut.PosRX_BOARD_ID), cfgOut.numTx, 1).';
    
    virtual_array.azi_arr = reshape(virtual_azi_arr,[],1);
    virtual_array.ele_arr = reshape(virtual_ele_arr,[],1);
    virtual_array.tx_id = reshape(repmat(cfgOut.PosTX_Trans_ID, cfgOut.numRx, 1),[],1);
    virtual_array.rx_id = reshape(repmat(cfgOut.PosRX_BOARD_ID, cfgOut.numTx, 1).',[],1);
    virtual_array.virtual_arr = cat(2,virtual_array.azi_arr,virtual_array.ele_arr,virtual_array.rx_id,virtual_array.tx_id);
    
   % 获取非冗余阵元对
    [~, noredundant_idx] = unique(virtual_array.virtual_arr(:,1:2),'rows'); % 对元素去重
    virtual_array.noredundant_arr = virtual_array.virtual_arr(noredundant_idx,:); % 获取非冗余阵列的方位、俯仰、收发ID
    virtual_array.noredundant_aziarr = virtual_array.noredundant_arr(virtual_array.noredundant_arr(:,2)==0,:); % 获取方位维度天线
    
    % 获取冗余阵元对
    redundant_idx = setxor([1:cfgOut.numRx*cfgOut.numTx],noredundant_idx); % 58 * 1 返回不在 A 和 B 的交集中的数据（对称差集），不包括重复项。
    virtual_array.redundant_arr = virtual_array.virtual_arr(redundant_idx,:); % 获取冗余阵列的方位、俯仰、收发ID
    
    % 查找并关联重叠的TX-RX对
    if ~isempty(redundant_idx) % 单板不存在冗余阵元
        info_overlaped_associate_arr = [];
        for re_arr_id = 1 : size(virtual_array.redundant_arr,1)
            % 这里是对冗余阵元和非冗余阵元进行匹配 输出结果为[[T,F,F,T],...]等形式 维度为134*4 生成掩码
            mask = virtual_array.noredundant_arr == virtual_array.redundant_arr(re_arr_id, :); 
            mask = mask(:,1) & mask(:,2); % 这里是对冗余阵元和非冗余阵元进行匹配 输出结果为[[T,F,F,T],...]等形式 维度为134*4
            info_associate = virtual_array.noredundant_arr(mask,:); % 找到冗余和非冗余中接收和发射天线重叠的元素
            info_overlaped = virtual_array.redundant_arr(re_arr_id,:); % 当前冗余阵元值
            % 将有关联的非冗余阵元和冗余阵元加入列表 重叠的意思指的是由不同TX-RX天线对形成的相位差是相同的
            info_overlaped_associate_arr = [info_overlaped_associate_arr; [info_associate, info_overlaped]];
        end
        diff_tx = abs(info_overlaped_associate_arr(:,8) - info_overlaped_associate_arr(:,4)); % 计算发射天线的位置差
        info_overlaped_diff1tx_arr = info_overlaped_associate_arr(diff_tx==1,:); % 32*8 发射天线差为1的重叠阵元
        [sotedVal, sortedIdx] = sort(info_overlaped_diff1tx_arr(:, 1)); % 按方位相位差排序
        virtual_array.info_overlaped_diff1tx = info_overlaped_diff1tx_arr(sortedIdx, :);
    else
        cfgOut.applyVmaxExtend = 0; % 允许速度扩展
    end
    
    sig_space_idx0 = virtual_array.noredundant_arr(:,1)+1;
    sig_space_idx1 = virtual_array.noredundant_arr(:,2)+1;
    sig_idx0 = [];
    sig_idx1 = [];
    rx_pos_set = virtual_array.noredundant_arr(:,3); % 接收天线位置集合
    tx_pos_set = virtual_array.noredundant_arr(:,4); % 发射天线位置集合
    for rx_id = 1 :length(virtual_array.noredundant_arr(rx_pos_set))
        rx_arr = rx_pos_set(rx_id);
        sig_idx0 = [sig_idx0, find(cfgOut.PosRX_BOARD_ID == rx_arr)]; % 确定接收天线的序号
    end
    for tx_id = 1 :length(virtual_array.noredundant_arr(tx_pos_set))
        tx_arr = tx_pos_set(tx_id);
        sig_idx1 = [sig_idx1, find(cfgOut.PosTX_Trans_ID == tx_arr)]; % 确定发射天线的顺序序号 不是板子上排列
    end
    cfgOut.sigSpaceIdx = [sig_space_idx0.'; sig_space_idx1.'];
    cfgOut.sigIdx = [sig_idx0; sig_idx1];
    
    cfgOut.virtual_array = virtual_array; % 将virtual_array参数存入更新
    
    cfgOut.arrdx = 0.5; % 归一化阵元间距
    cfgOut.arrdy = 0.5; 
    
    % 下面的是未考虑发射顺序前的版本
%     arrX = [];
%     arrY = [];
%     for id = 1 : cfgOut.numTx
%        arrX = [arrX PosRX_X + PosTX_X(id)]; % MIMO阵元排列
%        arrY = [arrY PosRX_Y + PosTX_Y(id)];
%     end   
%     cfgOut.array(1,:) = arrX;
%     cfgOut.array(2,:) = arrY;
    
    cfgOut.antennaPhase = zeros(1,cfgOut.numRx*cfgOut.numTx); % 初始化相位[可以结合幅相误差校准]
    
    
    % 信号参数
    cfgOut.startFreq = 76.5e9; % 起始频率
    cfgOut.fs = 6874e3; % ADC采样频率 Hz 6874e3
    cfgOut.Ramptime = 80e-6; % 工作时间 80e-6
    cfgOut.Idletime = 30e-6; % 空闲时间 30e-6
    cfgOut.Slope = 46.397e12; % Chirp斜率 46.397
    cfgOut.validB = 1 / cfgOut.fs * cfgOut.ADCNum * cfgOut.Slope; % 实际有效带宽
    cfgOut.totalB = cfgOut.Ramptime * cfgOut.Slope; % 完整带宽
    cfgOut.Tc = cfgOut.Idletime + cfgOut.Ramptime; % 单Chirp周期
    cfgOut.ValidTc = 1 / cfgOut.fs * cfgOut.ADCNum; % 单Chirp内ADC有效采样时间
    cfgOut.fc = cfgOut.totalB / 2 + cfgOut.startFreq + 1 / cfgOut.fs * cfgOut.ADCNum; % 载频 Hz 后面两项为有效的起始频率 
    
    % 硬件参数
    cfgOut.Pt = 12; % dbm 发射功率
    cfgOut.Fn = 12; % 噪声系数
    cfgOut.Ls = 3;  % 系统损耗
    
%     cfgOut.Pa = 48; % dB 接收链路增益
    
end