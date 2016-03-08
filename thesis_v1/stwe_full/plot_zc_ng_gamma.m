run ucsb_ms_chapt2
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
plot(freq,loss); hold on;
plot(freq,real(neff));
hold off;
xlabel('f (GHz)');
axis tight