%clear;
% setup
function stwe_setup()
clear passive active optical;
global passive active optical stwe_setup_flag

% % set passive Zc and gamma
% run F:\project\Hybrid_EAM\Hybrid_design\matlab\version_2\passive_part\TML_DATA
% Zcp = transpose(passive_data(:,3));
% gammap = transpose(passive_data(:,2))*1j;
% RcpDC = 1/4.54e7/1.5e-6/4e-6;
% clear passive_data
% 
% % set active Zc and gamma 
% run F:\project\Hybrid_EAM\Hybrid_design\matlab\version_2\structure_Park\PARK_DATA
% %run F:\project\Hybrid_EAM\Hybrid_design\matlab\version_2\structure_Park\Park_narrow_mesa\mesa2um\PARK_DATA_1p5um
% Zcm = transpose(active_data(:,3));
% gammam = transpose(active_data(:,2))*1j;
% clear active_data

%20100131 HFSS data
% set active Zc and gamma
flag_Rs_Ci = 0; % 0 using Rs = 3; Ci = 1pF; 1 using Calculated data
%run .\HFSS_0224\active_2013_6_3_hmetal_1
% run ucsb_ms_chapt2
run ucsb_chapt2
Zcm = transpose(active_data(:,2));
gammam = transpose(active_data(:,3));
freq = transpose(active_data(:,1))*1e9;
RcpDC = 1/4.1e7/(3e-6*8e-6+1e-6*4e-6+2e-6*2e-6);
clear active_data

run ucsb_chapt2
passive_data = active_data;
Zcp = transpose(passive_data(:,2));
gammap = transpose(passive_data(:,3));
RcpDC = 1/4.1e7/(3e-6*8e-6+1e-6*4e-6+2e-6*2e-6);
freq = transpose(passive_data(:,1))*1e9;
clear passive_data active_data

% passive parameters
passive.Zc=Zcp;% characteristic impedance
passive.Gamma = gammap;% propagation constant
passive.Rdc=RcpDC;%% DC resistance
passive.tao=0;%% ignore
passive.freq=freq;

% active parameters
active.Zc=Zcm;% characteristic impedance
active.Gamma=gammam;% propagation constant
active.Rdc=RcpDC;%% DC resistance
%Rs=1.8; % olms.mm
if flag_Rs_Ci == 1
	[Rsm,Cim] = Get_Rsm_Cim(freq,active.Zc,active.Gamma);
	Rsm = Rsm*1e3; %olms.mm
	Cim = Cim/1e3; %F/mm
elseif flag_Rs_Ci == 0
	Rsm=1.3; % olms.mm
	Cim=1e-12; % F/mm
end
active.Cim = Cim;
active.Rsm = Rsm;
active.tao=Rsm*Cim; %% Rs*Ci
%active.tao=0; %% 0
active.freq=freq;

% figure(200);
% M = 10*log10(abs((1./(1+1i.*Rsm.*Cim.*freq*2*pi))).^2);
% plot(freq/1e9,M,'LineWidth',2);
% title('Idea speed(no reflection)');xlabel('f(GHz)');ylabel('M(\omega)[dB]');
% sprintf('3dB:%.3f GHz\n',1/(2*pi*max(Rsm)*max(Cim))/1e9)
% optical parameters
optical.ng=3.6;

stwe_setup_flag=1;