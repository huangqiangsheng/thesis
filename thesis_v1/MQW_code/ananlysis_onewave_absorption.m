clear
load DataAbsorptionSweep
index =  find(lambda == 1520*1e-9);
lastnum = 0;
 figure(1);
 plot(Ebias(1:end-lastnum)/1e5,TMalpha(index,1:end-lastnum)./(max(TMalpha(index,1:end-lastnum))),'ko-',...
      Ebias(1:end-lastnum)/1e5,TEalpha(index,1:end-lastnum)./(max(TEalpha(index,1:end-lastnum))),...
      'r*-','LineWidth',2);
 title('Absorption @1550');
 xlabel('E field(KV/cm)');ylabel('Absorption(a.u.)');
 legend('TM','TE','Location','NorthWest');
 
 Vbias = (Ebias/1e5-27)/20;
 figure(2);
  plot(Vbias(1:end-lastnum),TMalpha(index,1:end-lastnum)./(max(TMalpha(index,1:end-lastnum))),'ko-',...
       Vbias(1:end-lastnum),TEalpha(index,1:end-lastnum)./(max(TEalpha(index,1:end-lastnum))),...
      'r*-','LineWidth',2);
 title('Absorption @1550');
 xlabel('Voltage(V)');ylabel('Absorption(a.u.)');
 legend('TM','TE','Location','NorthWest');
%  C = 10*log10(CouplingCoef*(1-R)^2);
 C = 0;%（dB） %insertion Loss 10dB
 Length = 100e-4;%(cm)
 To_TE = 0.38;%光限制因子
 To_TM = 0.38;
 TE_Tran = -4.343*To_TE*TEalpha*Length+C;
 TM_Tran = -4.343*To_TM*TMalpha*Length+C;
  figure(3);
  plot(Vbias(1:end-lastnum),TM_Tran(index,1:end-lastnum),'ko-',...
       Vbias(1:end-lastnum),TE_Tran(index,1:end-lastnum),...
      'r*-','LineWidth',2);
 title('Absorption @1550');
 xlabel('Voltage(V)');ylabel('Transmission(dB)');
 legend('TM','TE','Location','NorthEast');grid on