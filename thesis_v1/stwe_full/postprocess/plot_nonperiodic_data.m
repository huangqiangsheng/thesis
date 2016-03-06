% 20100114 by TYB
% collect result for the optimized structure and do postprocess
clear;
% load the configure.
run E:/博士-课题/EM/唐永波师兄论文/source_code/stwe/stwe_setup  % setup file including ZC, ZS, data.
freq = active.freq;

% CS=0;CL=0;
% Lens = [200]*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_TW=GSVm-GSVm0;
% S11_TW=Ta;
% 
% CS=0;CL=0;
% Lens = [60,200,620]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_S1=GSVm-GSVm0;
% S11_S1=Ta;
% 
% CS=0;CL=0;
% Lens = [120,70,	555,	600,	130, 760]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_S2=GSVm-GSVm0;
% S11_S2=Ta;
% 
% CS=0;CL=0;
% Lens = [105	64	375	405	60	450	595	76	265]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_S3=GSVm-GSVm0;
% S11_S3=Ta;
% 
% subplot(2,1,1);
% plot(freq, S21_TW, freq, S21_S1, freq, S21_S2, freq, S21_S3,'LineWidth',2);
% ylim([-4,0]);xlim([0,100]*1e9);grid on
% xlabel('Frequency');
% ylabel('Modulation response [dB]');
% legend('TW','1S','2S','3S');
% 
% subplot(2,1,2);
% plot(freq, S11_TW, freq, S11_S1, freq, S11_S2, freq, S11_S3,'LineWidth',2);
% ylim([-20,0]);xlim([0,100]*1e9);grid on
% xlabel('Frequency');
% ylabel('Reflection S11 [dB]');
% legend('TW','1S','2S','3S');

% CS=0;CL=0;
% Lens = [116.7	66.7	437.5	537.5	66.7	437.5	537.5  66.7	333.3]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_S3p=GSVm-GSVm0;
% S11_S3p=Ta;
% 
% CS=0;CL=0;
% Lens = [153.1	50	384.4	543.8	50	384.4	543.8	50  384.4	543.8	50 	453.1]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_S4p=GSVm-GSVm0;
% S11_S4p=Ta;
% 
% subplot(2,1,1);
% plot(freq, S21_S3p, freq, S21_S4p,'LineWidth',2);
% ylim([-4,0]);xlim([0,100]*1e9);grid on
% xlabel('Frequency');
% ylabel('Modulation response [dB]');
% legend('TW','S3p');
% 
% subplot(2,1,2);
% plot(freq, S11_S3p, freq, S11_S4p,'LineWidth',2);
% ylim([-20,0]);xlim([0,100]*1e9);grid on
% xlabel('Frequency');
% ylabel('Reflection S11 [dB]');
% legend('TW','S3p');

CS=0;CL=0;ZL=50;
Lens = [120	70	555	600	130		760]'*1e-6;
[GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
S21_Z50=GSVm-GSVm0;
S11_Z50=Ta;

CS=0;CL=0;ZL=45;
Lens = [105 77 490 510 123 455]'*1e-6;
[GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
S21_Z45=GSVm-GSVm0;
S11_Z45=Ta;

CS=0;CL=0;ZL=40;
Lens = [95 76 365 385 124 330]'*1e-6;
[GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
S21_Z40=GSVm-GSVm0;
S11_Z40=Ta;

CS=0;CL=0;ZL=35;
Lens = [135 64 300 320 136 200]'*1e-6;
[GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
S21_Z35=GSVm-GSVm0;
S11_Z35=Ta;

subplot(2,1,1);
plot(freq, S21_Z50, freq, S21_Z45,freq, S21_Z40, freq, S21_Z35,'LineWidth',2);
ylim([-4,0.5]);xlim([0,100]*1e9);grid on
xlabel('Frequency');
ylabel('Modulation response [dB]');
legend('Z_L=50','Z_L=45','Z_L=40','Z_L=35');

subplot(2,1,2);
plot(freq, S11_Z50, freq, S11_Z45,freq, S11_Z40, freq, S11_Z35,'LineWidth',2);
ylim([-20,0]);xlim([0,100]*1e9);grid on
xlabel('Frequency');
ylabel('Reflection S11 [dB]');
legend('Z_L=50','Z_L=45','Z_L=40','Z_L=35');

% CS=0;CL=0;ZL=50;
% Lens = [116.7	66.7	437.5	537.5	66.7 437.5	537.5	66.7 333.3]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_Z50p=GSVm-GSVm0;
% S11_Z50p=Ta;
% 
% CS=0;CL=0;ZL=45;
% Lens = [70.8	66.7	350	398	66.7	350	398	66.7 350	398	66.7	295.8]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_Z45p=GSVm-GSVm0;
% S11_Z45p=Ta;
% 
% CS=0;CL=0;ZL=40;
% Lens = [33.3	66.7	266.7	289.6	66.7	266.7	289.6	66.7 266.7	289.6	66.7	270.8]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_Z40p=GSVm-GSVm0;
% S11_Z40p=Ta;
% 
% CS=0;CL=0;ZL=35;
% Lens = [66.7	66.7	247.9	268.7	66.7	247.9	268.7	66.7  247.9	268.7	66.7	266.7]'*1e-6;
% [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
% S21_Z35p=GSVm-GSVm0;
% S11_Z35p=Ta;
% 
% subplot(2,1,1);
% plot(freq, S21_Z50p, freq, S21_Z45p,freq, S21_Z40p, freq, S21_Z35p,'LineWidth',2);
% ylim([-4,0.5]);xlim([0,100]*1e9);grid on
% xlabel('Frequency');
% ylabel('Modulation response [dB]');
% legend('Z_L=50','Z_L=45','Z_L=40','Z_L=35');
% 
% subplot(2,1,2);
% plot(freq, S11_Z50p, freq, S11_Z45p,freq, S11_Z40p, freq, S11_Z35p,'LineWidth',2);
% ylim([-20,0]);xlim([0,100]*1e9);grid on
% xlabel('Frequency');
% ylabel('Reflection S11 [dB]');
% legend('Z_L=50','Z_L=45','Z_L=40','Z_L=35');
