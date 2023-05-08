function [tarOut] = ConfigureTarget()
    %% 本文件用于仿真配置目标参数
    %% By Xuliang, 20230411
    
    tarOut.nums = 4; % 目标数目
    
    % 目标参数: 距离 速度 方位 俯仰 RCS 目标增益
    tarOut.range = [5 10 15 18]; % 距离 m
    tarOut.velocity = [0.4 -0.3 0 0.6]; % 速度 m/s
    tarOut.Azi = [15 -2 30 -25]; % 目标方位角度
    tarOut.Ele = [0 15 -5 8];     % 目标俯仰角度
    tarOut.RCS = [20 20 10 20];  % 目标RCS值
    tarOut.Gt  = [12 12 12 12];  % 发射增益
    tarOut.Gr  = [12 12 12 13];  % 接收增益
    tarOut.trueX = tarOut.range .* cosd(tarOut.Ele) .* sind(tarOut.Azi);
    tarOut.trueY = tarOut.range .* cosd(tarOut.Ele) .* cosd(tarOut.Azi);
    tarOut.trueZ = tarOut.range .* sind(tarOut.Ele);
    
end