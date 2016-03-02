function [MQW] = Absorption(MQW,lambda)
%% Using Variational Method for Exciton Problem
%% Reference "physics_of_optoelectronic_devices",shun_lien_chuang,
%   P563 - P566
%% "Modeling of self-electro-optic-effect devices" 
%   P. J. Mares and S. L. Chuang,1993
%% created by HQS 2011_12_4
%% corrected bug int the Fhw() function by HQS
%     load ExEffect;
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    epsr0 = 8.85418782e-12;
    mu0 = 1.256637e-6;
    evnm =1239.842;%
    %32 153TE   %82 767 %106 882

    S0H = 2;%1 (106)%2;(82)%2; %(32)%Smmerfeld enhancement factor. modified according to your expriment data
    S0L = 2;%1 (106)%6;(82)%6; %(32)
%     gammaH = (6-(106-MQW.Ebias/1e5)/(106-32)*(6-2))*1e-3*q;
%     gammaH = 1e-3*q; Ebias = 0
    gammaH = (1e-3+(MQW.Ebias/4e6)^2*0.4e-3)*q;
    gammaL = (1e-3+(MQW.Ebias/4e6)^2*0.4e-3)*q;
    me_start = MQW.well.me*0.14;
    Lp = MQW.tw;
%     gammaH =2e-3*q; %16e-3*q;  (106)%7e-3*q;(82)%2e-3*q;%(32)%(2.5e-3+(MQW.Ebias-32e5)/(106e5-32e5)*5.5e-3)*q;%2.5e-3*q(32),8e-3*q;(106) % Linewidth,modified according to your expriment data
%     gammaL =0.6e-3*q;%2.2e-3*q;(106)%2.2e-3*q;(82)%0.6e-3*q;%(32)%(2e-3+(MQW.Ebias-32e5)/(106e5-32e5)*1.5e-3)*q;%2.5e-3*q(32),3.5e-3*q(106)
    
    c0 = 1/sqrt(epsr0*mu0);
    nr = sqrt(MQW.well.epsr);
   lambda = lambda * 1e9;
%     lambda = linspace(1480,1580,100);
%     lambda = linspace(MQW.PL_cHH-20,MQW.PL_cHH+20,100);             
    omega = 2*pi./(lambda*1e-9./(c0/nr));
    muH =  1/(1/MQW.well.me + 1/MQW.well.mHHt) * m0;
    muL =  1/(1/MQW.well.me + 1/MQW.well.mLHt) * m0;
    C0CH = omega * muH * q^2 /(nr*c0*epsr0*pi*m0^2*Lp);
    C0CL = omega * muL * q^2 /(nr*c0*epsr0*pi*m0^2*Lp);
    C0D = 2*omega*q^2*hbar^2/(nr*c0*epsr0*m0^2*Lp);
    energy = evnm./lambda * q;
    EcvHE = (MQW.Ec-MQW.EvH+MQW.EbH)*q ;
    EcvLE = (MQW.Ec-MQW.EvL+MQW.EbL)*q;
    EcvH = (MQW.Ec-MQW.EvH+4e-3)*q;
    EcvL = (MQW.Ec-MQW.EvL)*q;
    InmH2 = abs(sum(conj(MQW.waveFunc).* MQW.waveFunH)*...
            abs(MQW.grid(2) - MQW.grid(1)))^2;
    InmL2 = abs(sum(conj(MQW.waveFunc).* MQW.waveFunL)*...
            abs(MQW.grid(2) - MQW.grid(1)))^2;
    phiH2 = 2/pi*(1/(MQW.ExwavelenH*1e-9))^2;
    phiL2 = 2/pi*(1/(MQW.ExwavelenL*1e-9))^2;
    phiH2 = 2/pi*(exp(-MQW.ExwavelenH*1e-9)./(-MQW.ExwavelenH*1e-9)).^2;
    phiL2 = 2/pi*(exp(-MQW.ExwavelenL*1e-9)./(-MQW.ExwavelenL*1e-9)).^2;
    MAVG = m0*(m0-MQW.well.me*m0)*MQW.well.Eg*(MQW.well.Eg+MQW.well.delta)/...
             (6*me_start*m0*(MQW.well.Eg+2*3/MQW.well.delta)) * q;
    
    TE_M0H = 3/2;
    TE_M0L = 1/2;
    TM_M0H = 0;
    TM_M0L = 2;
    TE_C1H = 3/4;
    TE_C2H = 3/4;
    TE_C1L = 5/4;
    TE_C2L = -3/4;
    TM_C1H = 3/2;
    TM_C2H = -3/2;
    TM_C1L = 1/2;
    TM_C2L = 3/2;
    EEn = MQW.EcEn*q;
    EHm = MQW.EvHm*q;
    ELm = MQW.EvLm*q;


 %%
    
    TE_alphaHD = Discrete(C0D,TE_M0H,MAVG,InmH2,phiH2,gammaH,EcvHE,energy);
    TE_alphaLD = Discrete(C0D,TE_M0L,MAVG,InmL2,phiL2,gammaL,EcvLE,energy);
    TM_alphaHD = Discrete(C0D,TM_M0H,MAVG,InmH2,phiH2,gammaH,EcvHE,energy);
    TM_alphaLD = Discrete(C0D,TM_M0L,MAVG,InmL2,phiL2,gammaL,EcvLE,energy);
    TE_alphaD = TE_alphaHD + TE_alphaLD;
    TM_alphaD = TM_alphaHD + TM_alphaLD;
%     figure(200)
%     subplot(2,1,1);
%     plot(lambda, abs(TE_alphaD));
%     xlabel('\lambda (nm)');
%     ylabel('Absorption Coefficient (1/cm)')
%     title('TE Exciton');
%     subplot(2,1,2);
%     plot(lambda, abs(TM_alphaD));
%     xlabel('\lambda (nm)');
%     ylabel('Absorption Coefficient (1/cm)')
%     title('TM Exciton');
 %%
    TE_alphaHC = Continuum(C0CH,InmH2,MAVG,TE_C1H,TE_C2H,EEn,EHm,S0H,gammaH,EcvH,energy);
    TE_alphaLC = Continuum(C0CL,InmL2,MAVG,TE_C1L,TE_C2L,EEn,ELm,S0L,gammaL,EcvL,energy);
    TM_alphaHC = Continuum(C0CH,InmH2,MAVG,TM_C1H,TM_C2H,EEn,EHm,S0H,gammaH,EcvH,energy);
    TM_alphaLC = Continuum(C0CL,InmL2,MAVG,TM_C1L,TM_C2L,EEn,ELm,S0L,gammaL,EcvL,energy);
    TE_alphaC = TE_alphaHC + TE_alphaLC;
    TM_alphaC = TM_alphaHC + TM_alphaLC;
%     figure(201)
%     subplot(2,1,1);
%     plot(lambda, abs(TE_alphaC));
%     xlabel('\lambda (nm)');
%     ylabel('Absorption Coefficient (1/cm)')
%     title('TE Continuum');
%     subplot(2,1,2);
%     plot(lambda, abs(TM_alphaC));
%     xlabel('\lambda (nm)');
%     ylabel('Absorption Coefficient (1/cm)')
%     title('TM Continuum');
 %%
    TE_alpha = TE_alphaD + TE_alphaC;
    TM_alpha = TM_alphaD + TM_alphaC;
    MQW.TEAlpha = TE_alpha;
    MQW.TMAlpha = TM_alpha;
    MQW.lambda = lambda;
    figure(203)
    subplot(2,1,1);
    plot(lambda, abs(TE_alpha),'LineWidth',2);
    xlabel('\lambda (nm)');
    ylabel('Absorption Coefficient (1/m)')
    title('TE Total');grid on
    subplot(2,1,2);
    plot(lambda, abs(TM_alpha),'LineWidth',2);
    xlabel('\lambda (nm)');
    ylabel('Absorption Coefficient (1/m)')
    title('TM Total');grid on
%     save Absorption_E105_v3 MQW;
end

function alpha = Discrete(C0,M0,MAVG,Inm2,phi2,gamma,EcvE,energy)
    alpha = C0 .* M0.* MAVG .* Inm2 .* phi2.*gamma./...
           ((EcvE)^2.*((EcvE-energy).^2+gamma.^2));
end

function alpha = Continuum(C0,Inm2,MAVG,C1,C2,Een,Ehm,S0,gamma,Ecv,energy)
   alpha = zeros(size(energy));
   q = 1.6021765e-19;  %C
   Een = Een/q;
   Ehm = Ehm/q;
   gamma = gamma/q;
   Ecv = Ecv/q;
   energy = energy/q;
   for iter = 1 : length(energy)
       f = quadl(@(x)Fhw(x,C1,C2,Een,Ehm,S0,gamma,Ecv,energy(iter)),0,1);
       alpha(iter) = f;
   end
   alpha = alpha .*C0 .* Inm2 .* MAVG / q^2;
end

function f = Fhw(x,C1,C2,Een,Ehm,S0,gamma,Ecv,energy)
    Et = tan(pi/2*x);
%   Et = x;
    f = pi/2*(1+Et.^2).*(C1+C2.*(Een+abs(Ehm))./(Een+abs(Ehm)+Et)).*S0 .*gamma./...
        ((Ecv+Et).^2.*((Ecv+Et-energy).^2 + gamma.^2));
        
%     f = (C1+C2.*(Een+Ehm)./(Een+Ehm+Et)).*S0 .*gamma./...
%         ((Ecv+Et).^2.*((Ecv+Et-energy).^2 + gamma.^2));
end