%2012/2/29 黄强盛
function [Vin,Vslop] = creat_PIN()
    SweepNum = 5;
%     Vin = linspace(0,3.6,SweepNum);
    Vin = linspace(0,-4,SweepNum);
    Vslop = zeros(size(Vin));
    for iter = 1:SweepNum
        Vslop(iter) = PIN(Vin(iter));
    end
    figure;
    plot(Vin,Vslop/1e5,'LineWidth',2);
    Vin = Vin';
    Vslop = Vslop';
    save Vin_2_Ein Vin Vslop
end
function Vslop = PIN(Vin)
%Vin输入电压强度
%Vslop内部电场强度
kB = 1.38065e-23;
Tem = 300; % temperature
h = 6.62606896e-34; % J s       Planck constant
hbar = h/2/pi;      % J s       reduced Planck constant
q = 1.6021765e-19;  %C
m0 = 9.10938215e-31; %kg
epsr0 = 8.85418782e-12;
P_Cladding = InP_params();
P_Cladding.len = -1.5e-6;
P_Cladding.Doping = 1e18;
P_SCH = InGaAsP_params(0.87,0.28,1550);
P_SCH.len = 0.15e-6;

load InGaAsP_MQW;

N_SCH = InGaAsP_params(0.87,0.28,1550);
N_SCH.len = 0.1e-6;
N_Cladding = InP_params();
N_Cladding.len = 1e-6;
N_Cladding.Doping = -3e18;

P_Cladding.Nv = 2*(P_Cladding.mHH*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6; %单位/cm^3
P_Cladding.Nc = 2*(P_Cladding.me*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;
N_Cladding.Nv = 2*(N_Cladding.mHH*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;
N_Cladding.Nc = 2*(N_Cladding.me*m0*kB*Tem/(2*pi*hbar^2))^(3/2)*1e-6;

V0 = P_Cladding.Eg0 + (N_Cladding.Ec0 - P_Cladding.Ec0) + ...
     kB*Tem/q*(log(P_Cladding.Doping/P_Cladding.Nv * ...
     -N_Cladding.Doping/N_Cladding.Nc));
% Vin = 0.5;%(V)输入电压强度
Ein = Vin/(MQW.t+N_SCH.len+P_SCH.len); % 输入电场强度;
V0 = V0 + Vin;
x = zeros(1,(1+2+MQW.num_w+MQW.num_b));  %xp是待求的解
x(2) = P_SCH.len;
epsr = zeros(1,(2+2+MQW.num_w+MQW.num_b));
epsr(1) = P_Cladding.epsr;
epsr(2) = P_SCH.epsr;
N = zeros(1,(2+2+MQW.num_w+MQW.num_b));
N(1) = P_Cladding.Doping;
N(2) = 0;
Ecband = zeros(1,(2+2+MQW.num_w+MQW.num_b));
Ecband(1) = P_Cladding.Ec0;
Ecband(2) = P_SCH.Ec0;
Evband = zeros(1,(2+2+MQW.num_w+MQW.num_b));
Evband(1) = P_Cladding.Ev0;
Evband(2) = P_SCH.Ev0;
for iter = 1: (MQW.num_w+MQW.num_b)
    if mod(iter,2) == 0
        x(iter+2) = x(iter+1) + MQW.tw;
        epsr(iter+2) = MQW.well.epsr;
        N(iter+2) = 0;
        Ecband(iter+2) = MQW.well.Ec0;
        Evband(iter+2) = MQW.well.Ev0;
    else
        x(iter+2) = x(iter+1) + 2*MQW.tb; % ATT!!!
        epsr(iter+2) = MQW.barrier.epsr;
        N(iter+2) = 0;
        Ecband(iter+2) = MQW.barrier.Ec0;
        Evband(iter+2) = MQW.barrier.Ev0;
    end
end
x(end) = x(end-1) + N_SCH.len;
epsr(end-1) = N_SCH.epsr;
epsr(end) = N_Cladding.epsr;
N(end-1) = 0;
N(end) = N_Cladding.Doping;
Ecband(end-1) = N_SCH.Ec0;
Ecband(end) = N_Cladding.Ec0;
Evband(end-1) = N_SCH.Ev0;
Evband(end) = N_Cladding.Ev0;
tol = 1e-16;
xp = linspace(0.001,0.1,200);
for iter = 1 :length(xp)
    DeltaV(iter) = DepletionWidth(xp(iter),x,epsr,N,V0);
end
figure(100);plot(xp,DeltaV);
xlabel('x(um)');hold on;
[gx,gy] = ginput(1);
plot([gx(1) gx(1)],[0 gy(1)],'r');
[xp0,fval,exitflag] = fminsearch(@(xp)DepletionWidth(xp,x,epsr,N,V0),gx(1), optimset('TolX',tol));
plot(xp0,fval,'r*');hold off
xn0 = xp0 * P_Cladding.Doping / N_Cladding.Doping;
fprintf('xp = %.3eum;\txn = %.3eum\n',xp0,-xn0);

interp_num = 4;
x0 =  [-xp0*1e-6 x -xn0*1e-6+x(end)];
[xc,Ec] = Draw_DepletionWidth(x0,epsr,N,Ecband,interp_num);
[xv,Ev] = Draw_DepletionWidth(x0,epsr,N,Evband,interp_num);
 figure(101);plot(xc*1e6,Ec,'LineWidth',2);hold on;plot(xv*1e6,Ev,'r','LineWidth',2);
 xlabel('x \mum');ylabel('eV');hold off;
Vslop = q*P_Cladding.Doping*1e6/(MQW.well.epsr*epsr0)*xp0*1e-6;
fprintf('Vin = %.3fV E = %.3fKV/cm\n',Vin,Vslop/1e5);
end
