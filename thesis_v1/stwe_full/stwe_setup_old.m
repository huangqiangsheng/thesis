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

run E:\Doctor\EM\电极设计\HFSS_0224\EAM_InGaAsP_v1
Zcm = transpose(active_data(:,2));
gammam = transpose(active_data(:,3));
freq = transpose(active_data(:,1))*1e9;
RcpDC = 1/4.1e7/(3e-6*4e-6+1e-6*4e-6+2e-6*2e-6);
clear active_data

% passive parameters
passive.Zc=50;% characteristic impedance
passive.Gamma = 3.5;% propagation constant
passive.Rdc=RcpDC;%% DC resistance
passive.tao=0;%% ignore
passive.freq=freq;

% active parameters
active.Zc=Zcm;% characteristic impedance
active.Gamma=gammam;% propagation constant
active.Rdc=RcpDC;%% DC resistance
%Rs=1.8; % olms.mm
Rs=3; % olms.mm
Ci=1e-12; % F/mm
active.tao=Rs*Ci; %% Rs*Ci
%active.tao=0; %% 0
active.freq=freq;

% optical parameters
optical.ng=3.5;

stwe_setup_flag=1;