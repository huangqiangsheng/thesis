clear
kB = 1.38065e-23;
Tem = 300; % temperature
h = 6.62606896e-34; % J s       Planck constant
hbar = h/2/pi;      % J s       reduced Planck constant
q = 1.6021765e-19;  %C
m0 = 9.10938215e-31; %kg
Vin = 1; % ‰»ÎµÁ—π
P_Cladding = AlGaAs_params(0.3,1550,300);
N_Cladding = AlGaAs_params(0,1550,300);
% P_Cladding = InP_params();
% N_Cladding = InP_params();
P_Cladding.Doping = 2e17;
% P_Cladding.Doping = 1e18;
P_Cladding.Nv = 2*(P_Cladding.mHH*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;
P_Cladding.Nc = 2*(P_Cladding.me*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;
N_Cladding.Doping = -4e16;
% N_Cladding.Doping = -3e18;
N_Cladding.Nv = 2*(N_Cladding.mHH*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;
N_Cladding.Nc = 2*(N_Cladding.me*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;

P_Cladding.Ec0 = N_Cladding.Ec0 - (N_Cladding.Eg0 - P_Cladding.Eg0)*0.67;
P_Cladding.Ev0 = N_Cladding.Ev0 + (N_Cladding.Eg0 - P_Cladding.Eg0)*0.33;

V0 = P_Cladding.Eg0 + (N_Cladding.Ec0 - P_Cladding.Ec0) + ...
     kB*Tem/q*(log(P_Cladding.Doping/P_Cladding.Nv * ...
     -N_Cladding.Doping/N_Cladding.Nc));
 V0 = V0 + Vin;
x = [0];
epsr = [P_Cladding.epsr,N_Cladding.epsr];
N = [P_Cladding.Doping,N_Cladding.Doping];
Ecband = [P_Cladding.Ec0,N_Cladding.Ec0];
Evband = [P_Cladding.Ev0,N_Cladding.Ev0];

tol = 1e-16;
xp = linspace(0.01,0.1,50);
for iter = 1 :length(xp)
    DeltaV(iter) = DepletionWidth(xp(iter),x,epsr,N,V0);
end
plot(xp,DeltaV);
xlabel('x(um)');hold on;
[gx,gy] = ginput(1);
plot([gx(1) gx(1)],[0 gy(1)],'r');
[xp0,fval,exitflag] = fminsearch(@(xp)DepletionWidth(xp,x,epsr,N,V0),gx(1), optimset('TolX',tol));
plot(xp0,fval,'r*');hold off
xn0 = xp0 * P_Cladding.Doping / N_Cladding.Doping;
fprintf('xp = %.3eum;\txn = %.3eum\n',xp0,-xn0);

interp_num = 20;
x0 =  [-xp0*1e-6 x -xn0*1e-6+x(end)];
[xc,Ec] = Draw_DepletionWidth(x0,epsr,N,Ecband,interp_num);
[xv,Ev] = Draw_DepletionWidth(x0,epsr,N,Evband,interp_num);
plot(xc*1e6,Ec);hold on;plot(xv*1e6,Ev,'r');
xlabel('um');ylabel('eV');