function  [sumClu]=dbscanClustering(eps,Obj,xFactor,yFactor,minPointsInCluster,FrameIndx)
   
    % 聚类判定准则：        (x-r0)^2/xFactor^2  + (y-y0)^2/yFactor^2 + (v-v0)^2 /vFactor^2 < eps^2   
    % maxClusters：        最大聚类数
    % minPointsInCluster： 一个cluster中最小点数
    % maxPoints：          一个cluster中允许最多点数
    % Obj标准数据结构       = [x,y,R,v，peakVal，SNR，aoaVar];
  
    %参数
    epsilon2_= eps*eps; % 邻域半径      
    numPoints =	size(Obj,1); % 点数
    visited = zeros(numPoints,1); % 访问矩阵，用于记录是否访问
    clusterId = 0; % 目标序号
    sumClu=[]; 
    
    colors = 'bgrcm';
    for i=1:numPoints
        
        if visited(i) == 0 %未标记点 将所有的点标记为0
            visited(i) = 1; %取第一个核心点q，标记为1
            tempIdx = i; % 用于标记当前点   
            x = Obj(i,1); % 用于获取当前点的x坐标
            y = Obj(i,2); % 用于获取当前点的y坐标

            numInEps=1; 
            for k=1:numPoints % 遍历所有点
                if visited(k) == 0 % 状态0  
                    summ = (Obj(k,1)-x)^2/xFactor^2+...
                         (Obj(k,2)-y)^2 /yFactor^2;
                    if summ < epsilon2_  % 如果满足距离准则，则该点属于此cluster
                        numInEps = numInEps+ 1; % 邻居节点数+1                  
                        tempIdx = [tempIdx k]; % 存储目标点
                    end
                end
            end        
            
            if numInEps > minPointsInCluster % 如果邻居节点数>要求的最小邻居节点数
                visited(i) = 1;   % 将第i点标记为核心点
                for in = 1:numInEps   % 遍历满足距离规则的点，并标记为已访问
                    visited(tempIdx(in)) = 1;   
                end
                
                next = 2;
                while next <= length(tempIdx) % 依次访问标记的节点
                    point_ref = tempIdx(next);
                    x = Obj(point_ref,1);
                    y = Obj(point_ref,2);
                    tempInd = [];
                    for ind=1:numPoints % 将标记节点作为核心节点，对所有点遍历
                        if visited(ind) == 0 % 状态0  
                            summ = (Obj(ind,1)-x)^2/xFactor^2+...
                                 (Obj(ind,2)-y)^2/yFactor^2;                            
                            if summ < epsilon2_ % 查看是否所有点都满足距离规则
                               tempInd = [tempInd ind]; % 满足则将该点加入
                            end
                        end
                    end

                    if length(tempInd) > minPointsInCluster % 如果核心点的簇内节点邻居也满足最小聚类点数
                        visited(point_ref) = 1;  % 标记簇内节点为核心节点，认为同类
                        numInEps = numInEps+ length(tempInd); % 计算簇内点数目
                        tempIdx = [tempIdx tempInd]; 
                        for kk = 1:length(tempInd)
                            visited(tempInd(kk)) = 1;  % 标记为已访问
                        end
                    else
                        visited(point_ref) = -1; %不满足最小聚类点数则认为是边界点
                    end
                    next = next+1;
                end
                
                tempClu = Obj(tempIdx,:); % obj = [ X,Y,h ,objSpeed,snr]; 获取核心点的属性信息
                cluLength = size(tempIdx,2);
              
                for pp = 1:cluLength % 对每个簇循环
                    ind = tempIdx(pp); 
                    if  visited(ind) == 1 % 绘制核心点
                         plot(tempClu(pp,1),tempClu(pp,2),'.','color', colors(mod(clusterId,length(colors))+1));
                         hold on;
                    elseif  visited(ind) == -1  % 绘制边界点
                           plot(tempClu(pp,1),tempClu(pp,2),'*','color', colors(mod(clusterId,length(colors))+1));                    
                         hold on;
                    end
                end
                title(['聚类结果：（*:边界点 , o：核心点 , .:噪声点)帧数：',num2str(FrameIndx)]);
                xlabel('点迹水平位置 ： m');
                ylabel('点迹垂直位置 ： m');
                hold on;
            else
                visited(i) = -2; %噪声点  
                cluLength = 1;
                plot(Obj(i,1),Obj(i,2),'r.');
                hold on;
            end
             
            if cluLength > 1 
                clusterId = clusterId+1; % New cluster ID
                output_IndexArray(tempIdx) = clusterId;
                sumClu(clusterId).numPoints = cluLength; % 簇内长度
                sumClu(clusterId).x_mean = mean(tempClu(:,1));  % 簇内点平均x坐标
                sumClu(clusterId).y_mean = mean(tempClu(:,2));  % 簇内点平均y坐标
%                 sumClu(clusterId).z_mean = mean(tempClu(:,3));  % 簇内点平均z坐标

               % 通过SNR加权
               sumClu(clusterId).x_SNR = (1./tempClu(:,end)')*tempClu(:,1)/sum(1./tempClu(:,end));   
               sumClu(clusterId).y_SNR = (1./tempClu(:,end)')*tempClu(:,2)/sum(1./tempClu(:,end)); 
%                sumClu(clusterId).z_SNR = (1./tempClu(:,end)')*tempClu(:,3)/sum(1./tempClu(:,end));

                % 仅取峰值
                [ ~,I] = max(1./tempClu(:,end)); % snr峰值
                tempx = tempClu(:,1);
                tempy = tempClu(:,2);
%                 tempz = tempClu(:,3);
          
                sumClu(clusterId).x_peak = tempx(I);    
                sumClu(clusterId).y_peak = tempy(I);
%                 sumClu(clusterId).z_peak = tempz(I);

                % y取最近位置 x取均值
                sumClu(clusterId).x_edge = mean(tempClu(:,1));   
                sumClu(clusterId).y_edge = min(tempClu(: ,2));
%                 sumClu(clusterId).z_edge = mean(tempClu(:,3));
     
                sumClu(clusterId).x = sumClu(clusterId).x_mean; % 提取平均x坐标  
                sumClu(clusterId).y = sumClu(clusterId).y_mean; % 提取平均y坐标  
%                 sumClu(clusterId).z = sumClu(clusterId).z_mean; % 提取平均z坐标  
        
                sumClu(clusterId).xsize = max(abs(tempClu(:,1) - sumClu(clusterId).x));
                sumClu(clusterId).ysize = max(abs(tempClu(:,2) - sumClu(clusterId).y));
%                 sumClu(clusterId).zsize = max(abs(tempClu(:,3) - sumClu(clusterId).z));
                sumClu(clusterId).v = mean(tempClu(:,end-1)); % 速度
                sumClu(clusterId).head =  sumClu(clusterId).y - sumClu(clusterId).ysize/2;

                if sumClu(clusterId).xsize < 1e-3
                    sumClu(clusterId).xsize = 1;
                end
                if sumClu(clusterId).ysize < 1e-3
                    sumClu(clusterId).ysize = 1;
                end

                if cluLength>1
                   sumClu(clusterId).centerRangeVar = var(sqrt(tempClu(:,1).^2+tempClu(:,2).^2)); % 距离方差
                   sumClu(clusterId).centerAngleVar = mean(deg2rad(atand(tempClu(:,2)/(tempClu(:,1)+eps))).^2); 
                   sumClu(clusterId).centerDopplerVar = var(tempClu(:,end-1)); % 多普勒方差

                else
                   sumClu(clusterId).centerRangeVar = 1;
                   sumClu(clusterId).centerAngleVar = 1;
                   sumClu(clusterId).centerDopplerVar = 1;
                end
            else
               continue;
            end
 
         end
    end   

    grid on
%     xlim([-100,100])
%     ylim([0,100])
    hold off
end