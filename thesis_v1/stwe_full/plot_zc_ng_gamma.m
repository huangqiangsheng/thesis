run ucsb_chapt2
freq = active_data(:,1);
Zc = active_data(:,2);
loss = active_data(:,4);
neff = active_data(:,5);
figure(1)
plot(freq,real(Zc)); hold on;
plot(freq,imag(Zc));
hold off;
xlabel('f (GHz)');
ylabel('Characteristic Impedance Z_C(\Omega)');
axis tight
figure(2)
[AX,H1,H2] = plotyy(freq,loss,freq,real(neff)); hold on;
HH1=get(AX(1),'Ylabel');
set(HH1,'String','Loss (dB/mm)');
set(HH1,'Color','k');
HH2=get(AX(2),'Ylabel');
set(HH2,'String','n_m');
set(HH2,'Color','k');
% set(H1,'Color','k');
% set(H2,'Color','k');
% set(AX(1),'XLim',[0 100]);
% set(AX(2),'XLim',[0]);
% set(AX(1),'xtick',0:10:100);
% set(AX(2),'xtick',0);
hold off;
xlabel('f (GHz)');
% axis tight