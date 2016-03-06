clear;
Lens=[278	150	883.8]'*1e-6;
CS=15e-15;CL=15e-15;ZS=50;ZL=50;
global passive active optical stwe_setup_flag;

stwe_setup;  % setup file including ZC, ZS, data.

[GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);
figure; 
[AX,H1,H2]=plotyy(active.freq, GSVm, active.freq, Ta);
grid on;
set(AX(2),'YLim',[-20,-5]);
set(AX(1),'YLim',[-3.5,0.3]);
set(AX(1),'YTickMode', 'auto');
set(AX(2),'YTickMode', 'auto');

% S-paramters at the electrical ports
Aa=ABCD(1,:);
Ba=ABCD(2,:);
Ca=ABCD(3,:);
Da=ABCD(4,:);
[s11,s12,s21,s22] = ABCD2S(Aa,Ba,Ca,Da,ZS,ZL);
S11 = 20*log10(abs(s11)); S21 = 20*log10(abs(s21));%convert Np/m to dB/m

%calculate the equavelent chacteristic impedance for the whole electrode.
delta = (Aa-Da)./Ca/2;
ZZ = sqrt((Aa+Da).*(Aa+Da)-4)./2./Ca;
ZZ = ((real(ZZ)>0)*2-1).*ZZ;
Zc1 = delta+ZZ;
Zc2 = -(delta-ZZ);
% electrode loss
gam = acosh(Aa/2+Da/2);
gam_i = find(imag(gam)<0);
gam(gam_i:end) = gam(gam_i:end)+j*2*pi;
if length(Lens)==1
    tLen=[0,Lens,0];
else
    tLen=Lens;
end    
tLen(3:3:end-1)=[];
gam = gam./sum(tLen); 
alpha = 20*log10(1)*real(gam); % attenuation per unit
% if need microwave index, should kick off the two leading electrodes.