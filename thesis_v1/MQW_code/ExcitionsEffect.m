function [MQW] = ExcitionsEffect(MQW)
%% Using Variational Method for Exciton Problem
%% Reference "physics_of_optoelectronic_devices",shun_lien_chuang,
%   P561 - P563
%% "Modeling of self-electro-optic-effect devices" 
%   P. J. Mares and S. L. Chuang,1993
%% created by HQS 2011_12_4

%     load FKE;
    muH =  1/(1/MQW.well.me + 1/MQW.well.mHHt);
    muL =  1/(1/MQW.well.me + 1/MQW.well.mLHt);
%     figure(1);
%     plot(MQW.grid,abs(MQW.waveFunc));hold on;
%     plot(MQW.grid,abs(MQW.waveFunH),'r');hold off;
    [EbH,wavelenH] = VarMethod(MQW.grid,MQW.grid,abs(MQW.waveFunc),abs(MQW.waveFunH),muH,MQW.well.epsr);
    [EbL,wavelenL] = VarMethod(MQW.grid,MQW.grid,abs(MQW.waveFunc),abs(MQW.waveFunL),muL,MQW.well.epsr);
    MQW.EbH = EbH;
    MQW.ExwavelenH = wavelenH*1e9;
    MQW.EbL = EbL;
    MQW.ExwavelenL = wavelenL*1e9;
    fprintf('H: Eex=%.3fmeV; wavelength = %.5f nm\n', EbH*1e3, wavelenH*1e9);
    fprintf('L: Eex=%.3fmeV; wavelength = %.5f nm\n', EbL*1e3, wavelenL*1e9);
    evnm =1239.842;% 1 eV = evnm/wl; evnm = hc/q*1e9
    MQW.PL_cHH = evnm./(MQW.Ec-MQW.EvH+EbH);
    MQW.PL_cLH = evnm./(MQW.Ec-MQW.EvL+EbL);
%     fprintf('H: Ec-HH : %f\n',MQW.Ec-MQW.EvH+EbH);
%     fprintf('L: Ec-HL:  %f\n',MQW.Ec-MQW.EvL+EbL);
    fprintf('Ebias: %.3f KV/cm\n',MQW.Ebias*1e-5);
    fprintf('PL_cHH: %fnm\n',MQW.PL_cHH);
    fprintf('PL_CHL: %fnm\n',MQW.PL_cLH);
%     save ExEffect MQW

end