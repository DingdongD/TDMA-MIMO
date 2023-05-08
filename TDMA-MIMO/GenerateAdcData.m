function [adcData] = GenerateAdcData(tarOut, cfgOut, IQFlag, saveFlag, dataPath)
    %% IQFlag: 生成信号模式标记
    %% saveFlag: 将生成数据进行保存
    %% 本文件用于双路ADC-IQ信号
    %% By Xuliang, 20230411
    
    % ADC数据生成
    if IQFlag
        tic
        adcData = GenerateSigIQ(tarOut, cfgOut);
        toc
    else
        tic
        [adcData] = GenerateSigI(tarOut, cfgOut);
        toc
    end
    disp(strcat(['=====','ADC信号生成完毕','====='])); 
    
    if saveFlag
        save(dataPath, 'adcData');
    end
end