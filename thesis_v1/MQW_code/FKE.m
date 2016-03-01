function MQW = FKE(MQW,Ebias)
%% Example to calculate the AlGaAs MQW parameters.
%% Combine airy funciton method and TMM method
%% Reference "physics_of_optoelectronic_devices",shun_lien_chuang,
%   P157 - P160
%% "Exact calculations of quasibound states of an isolated quantum
%   well with uniform electric field: Quantum-well Stark resonance", 
%   D.Ahn and,S.L.Chuang,1986
%% created by HQS 2011_12_4
%% modified by HQS 2012_1_16 correct TMM_EXCAL Function
% load ./creat/InGaAsP; %modified
tol=1e-32; 
sample_num = 200;
x = linspace(-MQW.t/2, MQW.t/2, sample_num); %sample point
dx = x(2) - x(1);
% Ebias = 105e5; % applied electric field
Vb = (1:sample_num)*Ebias*dx - fix(sample_num/2)*Ebias*dx; %applied bias potential V.
MQW.grid = x; %grid point
%%
v0 = MQW.barrier.Ec -MQW.well.Ec;
v = ones(1,sample_num)*v0;
m = ones(1,sample_num)*MQW.barrier.me;
for iter = 1:MQW.num_w;
    v(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = 0;
    m(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = MQW.well.me;
end
v = v + Vb;
[waveFunc, Ec] = TMM_EXCAL(x, v, m, 4,tol,v0,MQW.tw,MQW.tb, MQW.well.me,MQW.barrier.me,Ebias,0);
Ec = real(Ec); % ignore the escaple time
MQW.EcEn = Ec; % potential of electron
figure(1);

subplot(3,1,1);
[AX,H1,H2] = plotyy(x*1e9,real(waveFunc),x*1e9,v+ MQW.well.Ec);
set(get(AX(1),'Ylabel'),'String','normalize \phi_e(z)');
line([min(x*1e9),max(x*1e9)],[Ec+ MQW.well.Ec,Ec+ MQW.well.Ec],'Parent',AX(2),'Color','k')
set(H1,'Color','b');
set(get(AX(2),'Ylabel'),'String','eV');
title('Electron')
ylim(AX(1),[0,1.5e4]);
ylim(AX(2),[min(v+ MQW.well.Ec)-0.01 max(v+ MQW.well.Ec)+0.01]);
xlabel('z(nm)')
%%
v0 = MQW.barrier.EvH -MQW.well.EvH;
v = ones(1,sample_num)*v0;
m = ones(1,sample_num)*MQW.barrier.mHHz;
for iter = 1:MQW.num_w;
    v(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = 0;
    m(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) =  MQW.well.mHHz;
end
v = v + Vb;
[waveFunH, EvH] = TMM_EXCAL(x, -v, m, 5,tol,-v0,MQW.tw,MQW.tb, MQW.well.mHHz,MQW.barrier.mHHz,Ebias,1);
EvH = real(EvH); % ignore the escaple time
MQW.EvHm = EvH; % potential of heavy hole
figure(1);
subplot(3,1,2);
[AX,H1,H2] = plotyy(x*1e9,abs(waveFunH),x*1e9,v+MQW.well.EvH);
set(get(AX(1),'Ylabel'),'String','normalize \phi_v_H(z)'); 
line([min(x*1e9),max(x*1e9)],[-EvH+MQW.well.EvH,-EvH+MQW.well.EvH],'Parent',AX(2),'Color','k')
set(H1,'Color','b');
set(get(AX(2),'Ylabel'),'String','eV');
title('Heavry Hole')
ylim(AX(1),[0,1.5e4]);
ylim(AX(2),[min(v+MQW.well.EvH)-0.01 max(v+MQW.well.EvH)+0.01]);
xlabel('z(nm)')
%%
v0 = MQW.barrier.EvL -MQW.well.EvL;
v = ones(1,sample_num)*v0;
m = ones(1,sample_num)*MQW.barrier.mLHz;
for iter = 1:MQW.num_w;
    v(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = 0;
    m(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) =  MQW.well.mLHz;
end
v = v + Vb;
[waveFunL, EvL] = TMM_EXCAL(x, -v, m, 6,tol,-v0,MQW.tw,MQW.tb, MQW.well.mLHz,MQW.barrier.mLHz,Ebias,1);
EvL = real(EvL); % ignore the escaple time
MQW.EvLm = EvL;  % potential of light hole
figure(1);
subplot(3,1,3);
[AX,H1,H2] = plotyy(x*1e9,abs(waveFunL),x*1e9,v+MQW.well.EvL);
set(get(AX(1),'Ylabel'),'String','normalize \phi_v_L(z)');
line([min(x*1e9),max(x*1e9)],[-EvL+MQW.well.EvL,-EvL+MQW.well.EvL],'Parent',AX(2),'Color','k')
set(H1,'Color','b');
set(get(AX(2),'Ylabel'),'String','eV');
title('Light Hole')
ylim(AX(1),[0,1.5e4]);
ylim(AX(2),[min(v+MQW.well.EvL)-0.01 max(v+MQW.well.EvL)+0.01]);
xlabel('z(nm)')
%%
MQW.waveFunc = waveFunc;% wavefunction for electron
MQW.waveFunH = waveFunH;% wavefunction for heavy hole
MQW.waveFunL = waveFunL;% wavefunction for light hole

MQW.Ec = Ec + MQW.well.Ec;     % conduct band
MQW.EvH = -EvH + MQW.well.EvH;  % for heavy hole
MQW.EvL = -EvL + MQW.well.EvL;  % for light hole
evnm =1239.842;% 1 eV = evnm/wl; evnm = hc/q*1e9
MQW.PL_cHH = evnm./(MQW.Ec-MQW.EvH);
MQW.PL_cLH = evnm./(MQW.Ec-MQW.EvL);
MQW.Ebias = Ebias;
fprintf('Ebias:%.3fKV/cm\n',MQW.Ebias*1e-5);
fprintf('PL_cHH: %fnm\n',MQW.PL_cHH);
fprintf('PL_CHL: %fnm\n',MQW.PL_cLH);
InmH2 = abs(sum(conj(MQW.waveFunc).* MQW.waveFunH)*...
        abs(MQW.grid(2) - MQW.grid(1)))^2;
fprintf('InmH2: %.3f\n',InmH2);
InmL2 = abs(sum(conj(MQW.waveFunc).* MQW.waveFunL)*...
        abs(MQW.grid(2) - MQW.grid(1)))^2;
fprintf('InmL2: %.3f\n',InmL2);
% save FKE MQW;
end