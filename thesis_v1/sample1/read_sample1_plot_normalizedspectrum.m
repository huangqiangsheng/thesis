clc;clear;
filename = ['sample1_No.10_sweep_v_sweep_lambda_10_dB.csv'];
M = csvread(filename);
wl = M(1,2:end);
voltage = M(18:3:end,1);
power = M(19:3:end,2:end);
current = M(20:3:end,2:end);
% for iter = 1:length(voltage)
%     ind = find(power(iter,:)<-110);
%     power(iter,ind) = (power(iter,ind-1)+power(iter,ind+1))/2;
% end
zpower = max(power,[],1);
% normal_power = power;
normal_power = power-repmat(zpower,length(voltage),1);
fnormal_power= normal_power;
for iter = 1:length(voltage)
    fnormal_power(iter,:) = smooth(wl,normal_power(iter,:),0.1,'rloess');
end
% for iter = 1:6
%     figure(iter);
%     plot(wl,fnormal_power((iter-1)*7+1:iter*7,:)); 
%     xlabel('\lambda (nm)');
%     ylabel('Transmission variation(dB)')
%     legend(num2str(voltage((iter-1)*7+1:iter*7,:)),'Location','SouthEastOutside');
%     axis([1540,1570,-25,0])
% end
figure(7);
% subplot(1,2,1);
% h = pcolor(voltage,wl,fnormal_power'); 
h = pcolor(voltage,wl,(fnormal_power'));
colormap gray
set(h, 'LineStyle','none');
ylabel('\lambda (nm)');
xlabel('Bias (V)')
xlim([-3,1.0])
ylim([1540,1570])
% title('normalized transmission without filter')
% tv = [0.8, 0.4, -0.01, -1.0, -2.0];
% for iter = 1:length(tv);
%     ind = find(voltage >= tv(iter));
%     figure(8); hold on;
%     plot(wl,10.^(fnormal_power(ind(end),:)/10));
% end