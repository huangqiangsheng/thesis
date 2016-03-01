clear
load Absorption_E30_v3;
MQW_E30 = MQW;
load Absorption_E50_v3;
MQW_E50 = MQW;
load Absorption_E70_v3;
MQW_E70 = MQW;
load Absorption_E90_v3;
MQW_E90 = MQW;
load Absorption_E110_v3;
MQW_E110 = MQW;
evnm =1239.842;
figure(202)
subplot(2,1,1);
plot(MQW_E30.lambda, abs(MQW_E30.TEAlpha),...
     MQW_E50.lambda,abs(MQW_E50.TEAlpha),...
     MQW_E70.lambda,abs(MQW_E70.TEAlpha),...
     MQW_E90.lambda,abs(MQW_E90.TEAlpha),...
     MQW_E110.lambda,abs(MQW_E110.TEAlpha),...
     'k','LineWidth',2);
xlabel('\lambda (nm)');
ylabel('Absorption Coefficient (1/cm)')
title('TE Continuum');
grid on;
legend('30KV/cm','50KV/cm','70KV/cm','90KV/cm','110KV/cm','Location','NorthEastOutside');
subplot(2,1,2)
plot(MQW_E30.lambda, abs(MQW_E30.TMAlpha),...
     MQW_E50.lambda,abs(MQW_E50.TMAlpha),...
     MQW_E70.lambda,abs(MQW_E70.TMAlpha),...
     MQW_E90.lambda, abs(MQW_E90.TMAlpha),...
     MQW_E110.lambda,abs(MQW_E110.TMAlpha),...
      'k','LineWidth',2);
xlabel('\lambda (nm)');
ylabel('Absorption Coefficient (1/cm)')
title('TM Continuum');
grid on;
legend('30KV/cm','50KV/cm','70KV/cm','90KV/cm','110KV/cm','Location','NorthEastOutside');