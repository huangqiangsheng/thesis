% postprocess
% 20100131 by TYB, adapt to the new functions
% 20100114 by TYB
% collect result for the optimized structure and do postprocess
clear;
% load the configure.
stwe_setup();  % setup file including ZC, ZS, data.
global passive active optical;
freq = active.freq;

 CS=0;CL=0;ZS=50;ZL=50;
%  Lens = [0,100,0]*1e-6; %Rsm = 2
%  Lens = [185,100,570]*1e-6; %Rsm = 3 & 4
%  Lens = [240;66;400;400;84;515 ]*1e-6; %Rsm = 2 Lens = 150
%  Lens = [390;65;585;585;85;730 ]*1e-6; %Rsm = 3 Lens = 150
%  Lens = [270;45;400;400;55;530 ]*1e-6; %Rsm = 2 Lens = 100
%   Lens = [500;40;600;600;60;700 ]*1e-6; %Rsm = 3 Lens = 100
 Lens = 100e-6;
 [GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
 S21_TW=GSVm-GSVm0;
 S11_TW=Ta;

 BD3GSVm = get_f3dB_stwe(Lens,CS,CL,ZS,ZL,0,0);
 
figure(1)
plot(freq/1e9, S21_TW,'LineWidth',2);
line([0,100],[-3,-3],'Color','k')
ylim([-4,0.5]);xlim([0,100]);
xlabel('f(GHz)');
ylabel('E\O Response R_E_O [dB]');

figure(2)
plot(freq/1e9, S11_TW,'LineWidth',2);
ylim([-50,0]);
xlim([0,100]);
xlabel('f(GHz)');
ylabel('Microwave Reflection |S11| [dB]');


