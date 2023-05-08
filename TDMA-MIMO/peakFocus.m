function [rd_peak_list, rd_peak] = peakFocus(RDM, cfar_out_list)
    rd_peak = zeros(size(RDM)); % 用于存放峰值能量
    rd_peak_list = []; % 用于存放峰值索引
    data_length = size(cfar_out_list, 1); % 数据长度
    
    for target_idx = 1 :data_length 
       range_idx = cfar_out_list(target_idx,1);
       doppler_idx = cfar_out_list(target_idx,2);
       
       if range_idx > 1 && range_idx < size(RDM,1) && ...
               doppler_idx > 1 && doppler_idx < size(RDM,2)
           if RDM(range_idx,doppler_idx) > RDM(range_idx-1,doppler_idx) && ...
                   RDM(range_idx, doppler_idx) > RDM(range_idx+1,doppler_idx) && ...
                   RDM(range_idx, doppler_idx) > RDM(range_idx,doppler_idx-1) && ...
                   RDM(range_idx, doppler_idx) > RDM(range_idx,doppler_idx+1) 
               
               rd_peak(range_idx, doppler_idx) = RDM(range_idx,doppler_idx);
               rdlist = [range_idx; doppler_idx];
               rd_peak_list = [rd_peak_list rdlist];         
           end
       end
    end
end