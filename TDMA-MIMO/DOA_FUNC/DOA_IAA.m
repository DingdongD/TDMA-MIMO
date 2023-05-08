function PoutIAA = DOA_IAA(X, A, params)
    % 本程序用于实现IAA算法的3种变体形式
    % 参考文献:AReduced ComplexityApproach to IAA Beamforming for
    % Efficient DOA Estimation of Coherent Sources
    % By Xuliang, 20230301
    % X : 输入基带信号
    % A : 字典矩阵
    % params : (mode, iter_num1, iter_num2, beta_thres)
    % mode: 选择算法运行的模式 IAA-APES/IAA-ML/IAA-RC
    % iter_num1 : 第一轮迭代次数  iter_num2 第二轮迭代次数
    % beta_res : 门限系数【针对RC/RCML】
    
    [M, snap] = size(X); % 阵列数目 * 快拍
    thetaNum = size(A, 2); % 原子数目
    threshold = params.threshold;
    mode = params.mode;
    if strcmp(mode, "APES")
        iter_num = params.iter_num;
        P_old = (sum(abs(A' * X / M), 2) / snap).^2;
        for iter_id = 1 : iter_num
            R = A * spdiags(P_old, 0, thetaNum, thetaNum) * A'; 
            invR = inv(R); 

            P = zeros(thetaNum, 1); 
            for snapIdx = 1 :snap
                x = X(:, snapIdx); 
                invRx = invR * x; % M * 1
                invRa = invR * A; % M * thetaNum
                ainvRa = sum(conj(A) .* invRa, 1).'; 
                coeff = A' * invRx ./ real(ainvRa); 
                P = P + abs(coeff).^2;
            end
            P = P / snap;

            if norm(P_old - P) < threshold 
                break;
            end
            P_old = P;
        end

    elseif strcmp(mode, "ML")
        Ru = X * X' / snap;
        for k = 1 : thetaNum
            P(k) = pinv(A(:, k)' * pinv(Ru) * A(:, k));
        end
        
        Pold = diag(P);
        Rold = pinv(A * Pold * A');
        iter_num = params.iter_num;
        for iter_id = 1 : iter_num
            [~, idx] = sort(P);
            for k = 1 : thetaNum
               P_prev(idx(k)) = P(idx(k));
               P(idx(k)) = max(0, P(idx(k)) + (A(:, idx(k))' * Rold * (Ru - pinv(Rold)) * Rold * A(:, idx(k))) / (A(:, idx(k))' * Rold * A(:, idx(k)))^2);
               Rold = Rold - ((P(idx(k)) - P_prev(idx(k))) * Rold * A(:, idx(k)) * A(:, idx(k))' * Rold) * pinv(1 + (P(idx(k)) - P_prev(idx(k))) * A(:, idx(k))' * Rold * A(:, idx(k)));
            end
            
            if norm(P - P_prev) < threshold 
                break;
            end
        end
   
    elseif strcmp(mode, "RC")
        iter_num1 = params.iter_num1;
        iter_num2 = params.iter_num2;
        beta_thres = params.beta_thres;

        P = (sum(abs(A' * X / M), 2) / snap).^2;
        Pold = diag(P);
        Ru = X * X' / snap;
        for iter_id = 1 : iter_num1
            Rold = (A * Pold * A'); 
            for k = 1 : thetaNum
               w(k, :) = A(:, k)' * pinv(Rold) / (A(:, k)' * pinv(Rold) * A(:, k));
               P(k) = real(w(k, :) * Ru * w(k, :)');
            end
            Pold = diag(P);
        end

        [Pval, idx] = sort(P); % 升序
        [Rrow] = find(Pval > beta_thres * P(idx(1))); % 易误点: 不等式右边为排序后最小元素值乘以门限 左边为排序元素值 目的是找到第一个满足阈值的元素列标
        Ridx = Rrow(1);

        Q_left = (A(:, idx(1:Ridx-1)) * diag(P(idx(1 : Ridx-1))) * A(:, idx(1:Ridx-1))');
        for iter_id = 1 : iter_num2
            Q_right = (A(:, idx(Ridx:end)) * diag(P(idx(Ridx:end))) * A(:, idx(Ridx:end))'); 
            RR = Q_left + Q_right;
            for k = 1 : length(Rrow)
                w(idx(Rrow(k)), :) = A(:, idx(Rrow(k)))' * pinv(RR) / (A(:, idx(Rrow(k)))' * pinv(RR) * A(:, idx(Rrow(k))));
                P(idx(Rrow(k))) = w(idx(Rrow(k)), :) * Ru * w(idx(Rrow(k)), :)';
            end
        end
        
    else
        disp(["Check the input mode again!"]);
    end
           
   PoutIAA = P;
end
    
    