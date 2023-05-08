% function [PoutANM,u_vec,T] = DOA_ANM(Y, P)
%     % 本程序为vanilla-原子范数的函数实现文件
%     % Y ：基带信号
%     % A ：过完备基
%     % P ： 信源数目
%     f0 = 77e9; % 频率
%     c = 3e8; % 光速
%     lambda = c / f0; % 波长
%     d = lambda / 2; % 阵元间距
%     [M, snap] = size(Y); % 阵元 快拍
%     
%     if snap == 1 % 单快拍模型
%         sigma = 1; % 噪声参数
%         regular_param = sqrt(M * log(M * sigma));
%         cvx_begin sdp quiet
%         cvx_solver sdpt3
%             variable T(M, M) hermitian toeplitz
%             variable x
%             variable z(M,1) complex
%             minimize (regular_param * 0.5 *(x + T(1,1)) + 0.5 * norm(Y-z))
%             [x Y'; Y T] >= 0;
%         cvx_end
%         [Phi, Val] = rootmusic(T, P, 'corr');
%         Phis = Phi / 2 / pi ;
%         estimated_theta = asind(-Phis * lambda / d);
%     else % 多快拍模型
%         regular_param = sqrt(M * (snap + log(M) + sqrt(2 * snap * log(M))));
%         cvx_begin sdp quiet
%         cvx_solver sdpt3
%             variable T(M,M) hermitian toeplitz
%             variable X(snap, snap) hermitian
%             variable Z(M, snap) complex
%             minimize (regular_param * (trace(X) + trace(T)) + 1 / 2 * sum_square_abs(vec(Y - Z)));
%             [X Y';Y T] >= 0;
%         cvx_end
%         
%         [Phi, Val] = rootmusic(T, P, 'corr'); % 从u中恢复出相位值
%         Phis = Phi / 2 / pi ;
%         estimated_theta = asind(-Phis * lambda / d);
%     end
%     PoutANM = estimated_theta.';
%     u_vec = zeros(M, 1); % 构建下对角
%     for iid = 1 : M % 使用i去控制每个对角线
%         for pid = iid : M
%             u_vec(iid) = u_vec(iid) + T(pid , pid - iid + 1); % pid 控制列 iid=3 pid=3:8 pid-iid+1=1:6 
%         end
%         u_vec(iid) = u_vec(iid) / (M - iid + 1);
%     end
% end
%    

function PoutANM = DOA_ANM(Y, P)
    % 本程序为vanilla-原子范数的函数实现文件
    % Y ：基带信号
    % A ：过完备基
    % P ： 信源数目
    f0 = 77e9; % 频率
    c = 3e8; % 光速
    lambda = c / f0; % 波长
    d = lambda / 2; % 阵元间距
    [M, snap] = size(Y); % 阵元 快拍
    
    if snap == 1 % 单快拍模型
        sigma = 1; % 噪声参数
        regular_param = sqrt(M * log(M * sigma));
        cvx_begin sdp quiet
        cvx_solver sdpt3
            variable T(M, M) hermitian toeplitz
            variable x
            variable z(M,1) complex
            minimize (regular_param * 0.5 *(x + T(1,1)) + 0.5 * norm(Y-z))
            [x Y'; Y T] >= 0;
        cvx_end
%         [Phi, Val] = rootmusic(T, P, 'corr');
%         Phis = Phi / 2 / pi ;
%         estimated_theta = asind(-Phis * lambda / d);
    else % 多快拍模型
        regular_param = sqrt(M * (snap + log(M) + sqrt(2 * snap * log(M))));
        cvx_begin sdp quiet
        cvx_solver sdpt3
            variable T(M,M) hermitian toeplitz
            variable X(snap, snap) hermitian
            variable Z(M, snap) complex
            minimize (regular_param * (trace(X) + trace(T)) + 1 / 2 * sum_square_abs(vec(Y - Z)));
            [X Y';Y T] >= 0;
        cvx_end
%         [Phi, Val] = rootmusic(T, P, 'corr'); % 从u中恢复出相位值
%         Phis = Phi / 2 / pi ;
%         estimated_theta = asind(-Phis * lambda / d);
    end
%     PoutANM = estimated_theta;
    PoutANM = T;
end
    