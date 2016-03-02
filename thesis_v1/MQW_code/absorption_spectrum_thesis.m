clear
load ./creat/InAlGaAs_MQW_new; %modified
C = 0;%（dB） %insertion Loss 10dB
Length = 80e-6;%(cm)
To_TE = 0.25;%光限制因子
To_TM = 0.233;

 SweepNum = 3;
lambda = linspace(1400,1600,2001)*1e-9;
Ebias = [0,16,32]*1e5;
% Ebias = linspace(0,32e5,SweepNum);
% Vbias1 = -(Ebias/1e5-6.78)/11.30;
% Vbias2 = linspace(max(Vbias1),1.0,11);
% Vbias = 0;
% Ebias = 40e5;
% Ebias = 0;
TEalpha = zeros(length(lambda),length(Ebias));
TMalpha = TEalpha;
PL_cHH = zeros(1,SweepNum);
PL_cHH1 = PL_cHH;
PL_cLH = zeros(1,SweepNum);
PL_cLH1 = PL_cLH;
% index = find(lambda == 1540*1e-9);
for iter = 1:SweepNum
    MQW1 = FKE(MQW,Ebias(iter));
    PL_cHH1(iter) = MQW1.PL_cHH;
    PL_cLH1(iter) = MQW1.PL_cLH;
    MQW2 = ExcitionsEffect(MQW1);
    MQW3 = Absorption(MQW2,lambda);
    TEalpha(:,iter) = MQW3.TEAlpha;
    TMalpha(:,iter) = MQW3.TMAlpha;
    PL_cHH(iter) = MQW3.PL_cHH;
    PL_cLH(iter) = MQW3.PL_cLH;
%     figure(20);
%     subplot(2,1,1);
%     hold on;
%     plot((Ebias(iter)/1e5-27)/20,Absorption2dB(TMalpha(index,iter),Length,To_TM,C),'*');
%     subplot(2,1,2);
%     hold on;
%     plot((Ebias(iter)/1e5-27)/20,Absorption2dB(TEalpha(index,iter),Length,To_TE,C),'*');
end
figure(10);
plot(lambda*1e6, abs(TEalpha)/1e3,'LineWidth',2);
xlabel('\lambda (\mum)');
ylabel('TE Absorption Coefficient (1/mm)')
legend('0 KV/cm', '16 KV/cm', '32 KV/cm');
figure(11);
plot(lambda*1e6, abs(TMalpha)/1e6,'LineWidth',2);
xlabel('\lambda (\mum)');
ylabel('TM Absorption Coefficient (1/mm)')
legend('0 KV/cm', '16 KV/cm', '32 KV/cm');

% figure;plot(Ebias/1e5,[PL_cHH1;PL_cHH],'Linewidth',2);
% xlabel('Electric Field (KV/cm)');
% ylabel('Exicton Absorption Peak (nm)')
% xlim([0,32])
% figure;plot(Ebias/1e5,[PL_cLH1;PL_cLH],'Linewidth',2)
% xlabel('Electric Field (KV/cm)');
% ylabel('Exicton Absorption Peak (nm)')
% xlim([0,32])
% slop = dlambdadV();
% dlambda = lambda(2)-lambda(1);
% for iter = 1:(length(Vbias2)-1)
%    start_iter = round(-slop*(Vbias2(iter+1)-Vbias1(end))/dlambda);
%    TEalpha(1:end-start_iter+1,iter+SweepNum) = MQW3.TEAlpha(start_iter:end);
%    TMalpha(1:end-start_iter+1,iter+SweepNum) = MQW3.TMAlpha(start_iter:end);
% end
% TE_Tran = -4.343*To_TE*TEalpha*Length+C;
% TM_Tran = -4.343*To_TM*TMalpha*Length+C;
% figure;
% tran = 10.^(TE_Tran/10);
% h = pcolor(Vbias,lambda*1e9,TE_Tran);
% set(h, 'LineStyle','none');
% xlabel('\lambda (nm)');
% ylabel('Bias (V)')
% % xlim([-3,0.7])
% ylim([1540,1570])
% save absorption_spectrum3 tran Vbias lambda TE_Tran
% figure(3);
% plot(lambda*1e6,TE_Tran(:,1:end-lastnum),...
%       'r*-','LineWidth',2);
%  xlabel('wavelength(\mum)');ylabel('Transmission(dB)');
%  legend('TM','TE','Location','NorthWest');
%  