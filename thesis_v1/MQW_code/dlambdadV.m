function slop = dlambdadV()
const = init_const();
well = InAlGaAs_params(0.69,0.09,1550);
well.me = 0.05;
dx = 110e-10;
V = 100e-6*2.3e-6*0.5e-6;
dIdV = 0.015;
eta = 0.3;
tau = 1e-9; %carrier lifetime
lambda0 = 1.55e-6;
slop = -lambda0^2*dx*const.hbar*eta*tau*dIdV...
            /(const.c0*well.me*const.m0*const.q*V);
end