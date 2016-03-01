clear
load ./creat/InAlGaAs; %modified
C = 0;%（dB） %insertion Loss 10dB
Length = 80e-4;%(cm)
To_TE = 0.38;%光限制因子
To_TM = 0.38;
 
SweepNum = 1;
lambda = linspace(1500,1600,201)*1e-9;
% Ebias = linspace(-50,120,SweepNum)*1e5;
Ebias = 0;
TEalpha = zeros(length(lambda),SweepNum);
TMalpha = TEalpha;
PL_cHH = zeros(1,SweepNum);
PL_cHL = zeros(1,SweepNum);
index = find(lambda == 1540*1e-9);
for iter = 1:SweepNum
    MQW1 = FKE(MQW,Ebias(iter));
    MQW2 = ExcitionsEffect(MQW1);
    MQW3 = Absorption(MQW2,lambda);
    TEalpha(:,iter) = MQW3.TEAlpha;
    TMalpha(:,iter) = MQW3.TMAlpha;
    PL_cHH(iter) = MQW3.PL_cHH;
    PL_cHL(iter) = MQW3.PL_cLH;
    figure(20);
    subplot(2,1,1);
    hold on;
    plot((Ebias(iter)/1e5-27)/20,Absorption2dB(TMalpha(index,iter),Length,To_TM,C),'*');
    subplot(2,1,2);
    hold on;
    plot((Ebias(iter)/1e5-27)/20,Absorption2dB(TEalpha(index,iter),Length,To_TE,C),'*');
end


lastnum = 0;
%  figure(1);
%  plot(Ebias(1:end-lastnum)/1e5,TMalpha(index,1:end-lastnum)./(max(TMalpha(index,1:end-lastnum))),'ko-',...
%       Ebias(1:end-lastnum)/1e5,TEalpha(index,1:end-lastnum)./(max(TEalpha(index,1:end-lastnum))),...
%       'r*-','LineWidth',2);
%  title('Absorption @1550nm');
%  xlabel('E field(KV/cm)');ylabel('Absorption(a.u.)');
%  legend('TM','TE','Location','NorthWest');
%  
%  Vbias = (Ebias/1e5-28.57)/21;
%  figure(2);
%   plot(Vbias(1:end-lastnum),TMalpha(index,1:end-lastnum)./(max(TMalpha(index,1:end-lastnum))),'ko-',...
%        Vbias(1:end-lastnum),TEalpha(index,1:end-lastnum)./(max(TEalpha(index,1:end-lastnum))),...
%       'r*-','LineWidth',2);
%  title('Absorption @1550nm');
%  xlabel('Voltage(V)');ylabel('Absorption(a.u.)');
%  legend('TM','TE','Location','NorthWest');
%  save DataAbsorptionSweep PL_cHH PL_cHL lambda Ebias TEalpha TMalpha;
% %  C = 10*log10(CouplingCoef*(1-R)^2);

 TE_Tran = -4.343*To_TE*TEalpha*Length+C;
 TM_Tran = -4.343*To_TM*TMalpha*Length+C;
  figure(3);
  plot(Vbias(1:end-lastnum),TM_Tran(:,1:end-lastnum),'ko-',...
       Vbias(1:end-lastnum),TE_Tran(:,1:end-lastnum),...
      'r*-','LineWidth',2);
 title('Absorption @1550nm');
 xlabel('Voltage(V)');ylabel('Transmission(dB)');
 legend('TM','TE','Location','NorthWest');
 