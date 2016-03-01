% clear
load absorption_spectrum
figure(7);
% subplot(1,2,2)
h = pcolor(Vbias,lambda*1e9,10*log10(tran));
set(h, 'LineStyle','none');
% colormap gray
ylabel('\lambda (nm)');
xlabel('Bias (V)')
xlim([-3,1])
ylim([1540,1570])