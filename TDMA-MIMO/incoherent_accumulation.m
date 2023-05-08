function [accumulateRD] = incoherent_accumulation(RDM)
    %% 本文件用于实现非相干积累
    %% By Xuliang,20230412
    % RDM：ADCNum * ChirpNum * ARRNum
    
    accumulateRD = squeeze(sum(abs(RDM), 3))/sqrt(size(RDM,3));
end