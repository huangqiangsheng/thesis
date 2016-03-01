% plot(Vbias(1:end),PL_cHH,'ko-','LineWidth',2);
% title('Absorption Peak V.S. Applied Voltage');
% xlabel('Voltage(V)');ylabel('Absorption Peak \lambda(nm)');
figure(1)
plot(lambda*1e9,TE_Tran(:,1:1:9),'LineWidth',2); 
xlabel('\lambda(nm)');ylabel('Transmission(dB)');
figure(2)
plot(lambda*1e9,TE_Tran(:,1:2:9),'LineWidth',2); 
xlabel('\lambda(nm)');ylabel('Transmission(dB)');
for iter = 1:length(Vbias)
    sV(iter,:) = [num2str(Vbias(iter),'%.1f') 'V'];
end
legend(sV(1,:),sV(3,:),sV(5,:),sV(7,:),sV(9,:));