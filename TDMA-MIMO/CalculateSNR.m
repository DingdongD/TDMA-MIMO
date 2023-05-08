function [tarSNR] = CalculateSNR(targetR, targetRCS, Gt, Gr, lambda, Pt, Fn, Ls, FsAdc)
%% 本文件用于仿真配置目标参数
%% By Xuliang, 20230411

    Pt_dB = Pt - 30; % 发射功率 dB
    R = db(targetR) / 2; % 距离
    lambda = db(lambda ^ 2) / 2; % 波长
    const = db((4 * pi) ^ 3) / 2; % 常数项
    
    K = 1.380649e-23; % 玻尔兹曼常数 J/k 
    B = FsAdc; % 采样带宽 HZ
    T0 = 290; % 开尔文系数
    KBT = db(K * B * T0) / 2; 
    
    tarSNR = (Pt_dB + Gt + Gr + targetRCS + lambda) - (KBT + Fn + const + 4 * R + Ls);
    
end