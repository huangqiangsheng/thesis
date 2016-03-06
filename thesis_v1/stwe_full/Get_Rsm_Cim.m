function [Rsm,Cim] = Get_Rsm_Cim(freq,Zc,Gamma)
% waveguide parameter
% still has some problem
% Wa = 0.5e-6;
% Wg = 8.5e-6;%10e-6;
% delta_InGaAs = 69183;%20
% delta_p = 1530;
% delta_n = 37000;%109920;%37000;%110400
% delta_sub_n = 37000;%109920;%37000;%110400
% hp = 1.5e-6;
% hn = 0;
% hInGaAs = 0.1e-6;
% d = 0.15e-6;

Wa = 0.6e-6;
Wg = 8.5e-6;%10e-6;
delta_InGaAs = 69183;%20
delta_p = 1530;
delta_n = 37000;%109920;%37000;%110400
delta_sub_n = 37000;%109920;%37000;%110400
hp = 0.95e-6;
hn = 0;
hInGaAs = 0.1e-6;
d = 0.15e-6;

[Rc,Lm,Rsm,Cim,Ce] = EqualCirleElement...
         (freq,Zc,Gamma,Wa,Wg,delta_InGaAs,delta_p,delta_n,delta_sub_n,...
         hInGaAs,hp,hn,d);
figure(100);
fprintf('Rsm= %.5f olms.mm\n',Rsm*1000)
subplot(2,2,1);plot(freq/1e9,Rc./1000./sqrt(freq/1e9),'LineWidth',2);
title('Rc');xlabel('f(GHz)');ylabel('Rc(\Omega/(mm*GHz^1^/^2^)');
subplot(2,2,2);plot(freq/1e9,Lm,'LineWidth',2);
title('Lm');xlabel('f(GHz)');ylabel('Lm(H/m)');
subplot(2,2,3);plot(freq/1e9,Cim/1e3/1e-12,'LineWidth',2);
title('Cim');xlabel('f(GHz)');ylabel('Cim(pF/mm)');
subplot(2,2,4);plot(freq/1e9,Ce/1e3/1e-12,'LineWidth',2);
title('Ce');xlabel('f(GHz)');ylabel('Ce(pF/mm)');
figure(200);
M = 10*log10(abs((1./(1+1i.*Rsm.*Cim.*freq*2*pi))).^2);
% plot(freq/1e9,M,'LineWidth',2);
% title('Idea speed(no reflection)');xlabel('f(GHz)');ylabel('M(\omega)[dB]');
end