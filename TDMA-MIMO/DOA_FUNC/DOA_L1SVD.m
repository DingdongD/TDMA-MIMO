function PoutSVD = DOA_L1SVD(X, A, P)
    % 本程序为L1-SVD的函数实现文件
    % X ：基带信号
    % A ：过完备基
    % P ： 信源数目
    
    [M, snap] = size(X); % M 阵元 snap 快拍
    thetaNum = size(A, 2); % 扫描网格点数
    DK1 = eye(P);
    DK2 = zeros(P, snap-P);
    DK = [DK1, DK2].'; % SNAP * P的选择矩阵
    [U, Sigm, V] = svd(X); % 奇异值分解
    Xsv = X * V * DK; % 获取新的接收矩阵
    
    % 低信噪比测试方案
    cvx_begin quiet
        variables p q
        variables r(thetaNum)
        variable SSV1(thetaNum, P) complex
        expressions Rn(M, P) complex

        minimize(p + 2.7 * q); % 优化目标
        subject to
            Rn = Xsv - A * SSV1; % 求残差
            Rvec = vec(Rn); % 矩阵转换为向量
            norm(Rvec) <= p; % 第一个不等式约束
            sum(r) <= q; % 第二个不等式约束
            for i = 1 : thetaNum % 第三个不等式约束
                norm(SSV1(i, :)) <= r(i);
            end
    cvx_end
    
    % 高信噪比测试方案
    % confidence_interval = 0.9; % 置信值
    % noise = X - X0; % 噪声估计
    % noise_var = var(noise(:)); % 噪声方差估计
    % regulari_param = compute_regulariParam(confidence_interval, noise_var, M, P); % 根据卡方分布反演门限值
    % cvx_begin quiet
    %     variable SSV1(length(theta_grids), P) complex
    %     minimize(sum(norms(SSV1, 2, 2)))
    %     subject to
    %         norm(Xsv - A * SSV1, 'fro') <= regulari_param 
    % cvx_end
    PoutSVD = abs(SSV1(:, :)  / max(SSV1(:, :))); % 求解功率谱
    
end