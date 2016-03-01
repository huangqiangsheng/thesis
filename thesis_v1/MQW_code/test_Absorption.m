clear;clc
load test_Absorption
MQW3 = Absorption(MQW2,lambda);
C = 0;%（dB） %insertion Loss 10dB
Length = 80e-6;%(cm)
To_TE = 0.25;%光限制因子
figure(1);
plot(lambda*1e6,10.^(-Absorption2dB(MQW3.TEAlpha,Length,To_TE,C)/10))
xlim([1.54,1.57])
xlabel('lambda (\mum)');
ylabel('Absorption (dB)');