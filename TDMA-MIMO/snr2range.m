%% 雷达方程：画出在某一目标RCS值下不同距离下的SNR
clear all; close all; clc;


%参数设置以及单位转换
P_t = 12;       %dBm   10*log10(P/1mw)    [这里说的是从PA出来的功率，其实它会经过馈电网络的损耗再给到天线，实际发出的可能会小于该值]
G_t = 12;       %dBi   10*log10(P/1)
G_r = 12;       %dBi    
RCS = 10;       %dBsm  10*log10(RCS)      

C      = 3e8;
fre    = 77e9;     %Hz
lambda = C/fre;

nums_rangefft   = 256;
nums_dopplerfft = 128;
G_s             = nums_dopplerfft*nums_rangefft;   
% G_s = 1;
 
K   = 1.380649e-23;   %J/k      玻尔兹曼常数
T_0 = 290;            %k       
B   = 1.72e9;           %Hz      20e6 6874e3 
F   = 12;             %dB
L   = 3;              %dB     


%做单位转换：需要统一将dB转换成幅度
P_t = (10^(P_t/10))/1e3;    %Pt应该变成w
G_t = (10^(G_t/10));
G_r = (10^(G_r/10));
RCS = (10^(RCS/10));
F   = (10^(F/10));     
L   = (10^(L/10));



%计算以及看看不同距离下的SNR
R   = (0:0.1:180);
SNR = P_t * G_t * G_r * RCS * lambda^2 * G_s ./ (K * T_0 * B * F * ((4*pi)^3) * (R.^4) * L);
SNR = 10*log10(SNR);
figure(1)
plot(R,SNR);
title(['SNR与距离关系曲线 ', 'RCS = ',num2str(10*log10(RCS)),'dBsm']);xlabel('range(m)');ylabel('SNR(dB)');
grid on;





