%# calculate the parameters for InGaAsP
%% E. Herbert Li, Physica E (2000) 215-273

function mat = InGaAsP_params(w,v, wl)
%strained In_w Ga_1-w As_v P_1-v
%wl: wavelength [nm]
%flag:1 well 0 barrier
if nargin <3
    wl=1550;%nm
end

mat.c11 = InGaAsP_c11(w,v);
mat.c12 = InGaAsP_c12(w,v);

mat.gamma1= InGaAsP_gamma1(w,v);
mat.gamma2= InGaAsP_gamma2(w,v);

mat.a = InGaAsP_a(w,v);
mat.b = InGaAsP_b(w,v);
mat.d = InGaAsP_d(w,v);
mat.ac = InGaAsP_ac(w,v);
mat.av = InGaAsP_av(w,v);

%% lattice length
mat.a0 = InGaAsP_a0(w,v);
%% band parameters ignoring strain
mat.delta = InGaAsP_delta(w,v);
mat.Evav0 = InGaAsP_Evav(w,v);
mat.Eg0 = InGaAsP_Eg(w,v);
[mat.Ec0, mat.Ev0, mat.Eso0]=InGaAsP_band(w,v);
%% band parameters with strain
[mat.Ec,mat.Ev, mat.EvH, mat.EvL,mat.Eso]=InGaAsP_band_strain(w,v);
mat.Eg = mat.Ec - mat.Ev;
%% strain
[mat.exx,mat.eyy,mat.ezz]=InGaAsP_strain(w,v);
%% effective mass [m0]
mat.me = InGaAsP_me(w,v);
mat.mHH = InGaAsP_mHH(w,v);
mat.mLH = InGaAsP_mLH(w,v);
mat.mHHz = InGaAsP_mHHz(w,v);
mat.mLHz = InGaAsP_mLHz(w,v);
mat.mHHt = InGaAsP_mHHt(w,v);
mat.mLHt = InGaAsP_mLHt(w,v);
%% low frequency relative epsilon
mat.epsr = InGaAsP_epsr(w,v);
%% lifetime [s]
% [mat.taun, mat.taup] = InGaAsP_tau(w,v);
%% n, refractive index with strain
%% n0, refractive index without considering strain
mat.n0 = InGaAsP_n(w,v,wl);

end

function c11 = InGaAsP_c11(w,v)
    c11 = 11.9*(1 - w)*v + 8.329*w*v + 14.05*(1 - w)*(1 - v) + 10.11*w*(1 - v);
end
function c12 = InGaAsP_c12(w,v)
    c12 = 5.38*(1 - w)*v + 4.526*w*v + 6.203*(1 - w)*(1 - v) + 5.61*w*(1 - v);
end

function gamma1= InGaAsP_gamma1(w,v)
gamma1_GaAs = 6.8;
gamma1_GaP = 4.05;
gamma1_InAs = 20.4;
gamma1_InP = 4.95;
% gamma1 = 6.68*(1 - w)*v + 21.8*w*v + 4.05*(1 - w)*(1 - v) + 5.060*w*(1 - v);
gamma1 = gamma1_GaAs*(1-w)*v + gamma1_GaP*(1-w)*(1-v)+gamma1_InAs*w*v+gamma1_InP*w*(1-v);
% gamma1 = w*gamma1_InAs + (1-w-v)*gamma1_GaAs + v*gamma1_AlAs;
end
function gamma2= InGaAsP_gamma2(w,v)
gamma2_GaAs = 1.9;
gamma2_GaP = 0.49;
gamma2_InAs = 8.3;
gamma2_InP = 1.65;
gamma2 = gamma2_GaAs*(1-w)*v + gamma2_GaP*(1-w)*(1-v)+gamma2_InAs*w*v+gamma2_InP*w*(1-v);
end

function a = InGaAsP_a(w,v)
a_GaAs = -8.33;
a_GaP = -8.83;
a_InAs = -6.08;
a_InP = -6.31;

a = a_GaAs*(1-w)*v + a_GaP*(1-w)*(1-v)+a_InAs*w*v+a_InP*w*(1-v);

end

function b = InGaAsP_b(w,v)
b = -1.7*(1 - w)*v + (-1.8)*w*v + (-1.8)*(1 - w)*(1 - v) + (-2.0)*w*(1 - v);
end

function d = InGaAsP_d(w,v)
d_GaAs = -4.55;
d_GaP = -4.5;
d_InAs = -3.6;
d_InP = -5.6;

d = d_GaAs*(1-w)*v + d_GaP*(1-w)*(1-v)+d_InAs*w*v+d_InP*w*(1-v);
end
function a0 = InGaAsP_a0(w,v)
a0 = 5.6533*(1 - w)*v + 6.0583*w*v + 5.4505*(1 - w)*(1 - v) + 5.8687*w*(1 - v);
end
function Eg = InGaAsP_Eg(w,v)
% E. Herbert Li, Physica E (2000) 215-273
% Eg =1.35 - 1.17*v + 0.668*(1 - w) + 0.069*v*(1 - w) + 0.15*v*v +...
%     0.03*(1 - w)*v*v+0.785*(1 - w)*(1-w) + 0.758*(1 - w)*(1 - w) - 0.322*v*(1 - w)*(1 - w);
% Zhang PhD dissertation
%Eg = 1.35 + 1.09*(1-w) - v + 0.33 * (1-w)*(1-v) - (0.73-0.28*v)*(1-w)*w -...
%     (0.101+0.109*(1-w))*v*(1-v)+0.05*sqrt((1-w)*v*w*(1-v));
% Physics of opoelectronic devices shun lien chang
  Eg = 1.35 + 0.668*(1-w) - 1.068*v + 0.758*(1-w)*(1-w) + 0.078*v*v-...
       0.069*(1-w)*v - 0.332*(1-w)*(1-w)*v + 0.03*(1-w)*v*v;
end
function ac = InGaAsP_ac(w,v)
ac_GaAs = -7.17;
ac_GaP = -7.14;
ac_InAs = -5.08;
ac_InP = -5.04;

ac =  ac_GaAs*(1-w)*v + ac_GaP*(1-w)*(1-v)+ac_InAs*w*v+ac_InP*w*(1-v);
end
function av = InGaAsP_av(w,v)
av_GaAs = 1.16;
av_GaP = 1.70;
av_InAs = 1.00;
av_InP = 1.27;

av = av_GaAs*(1-w)*v + av_GaP*(1-w)*(1-v)+av_InAs*w*v+av_InP*w*(1-v);
end

function delta = InGaAsP_delta(w,v)
delta_GaAs = 0.34;
delta_GaP = 0.10;
delta_InAs = 0.371;
delta_InP = 0.10;
delta = delta_GaAs*(1-w)*v + delta_GaP*(1-w)*(1-v)+delta_InAs*w*v+delta_InP*w*(1-v);
end

function Evav = InGaAsP_Evav(w,v)
Evav_GaAs = -6.92;
Evav_GaP = -7.40;
Evav_InAs = -6.67;
Evav_InP = -7.04;

Evav = Evav_GaAs*(1-w)*v + Evav_GaP*(1-w)*(1-v)+Evav_InAs*w*v+Evav_InP*w*(1-v);
end
function [Ec, Ev, Eso]=InGaAsP_band(w,v)
Evav = InGaAsP_Evav(w,v);
delta = InGaAsP_delta(w,v);
Eg = InGaAsP_Eg(w,v);

Ev = Evav + delta/3;
Eso = Ev - delta;
Ec = Ev + Eg;
end

function [exx,eyy,ezz]=InGaAsP_strain(w,v)
a0_InP = 5.8688;%A
a0 = InGaAsP_a0(w,v);
exx = (a0_InP-a0)/a0;
eyy = exx;
c11 = InGaAsP_c11(w,v);
c12 = InGaAsP_c12(w,v);
ezz = -2*c12/c11*exx;
end

function [Ec,Ev, EvH, EvL,Eso]=InGaAsP_band_strain(w,v)
[Ec0, Ev0, Eso0]=InGaAsP_band(w,v);
Evav0 = InGaAsP_Evav(w,v);
[exx,eyy,ezz]=InGaAsP_strain(w,v);
delta = InGaAsP_delta(w,v);

Eg = InGaAsP_Eg(w,v);
av = InGaAsP_av(w,v);
b = InGaAsP_b(w,v);
ac = InGaAsP_ac(w,v);

Pe = -av*(exx+eyy+ezz);
Qe = -b/2*(exx+eyy-2*ezz);
Pc = ac*(exx+eyy+ezz);

Evav = Evav0 - Pe;
Ev = Evav + delta/3;
Ec = Ev0+Eg+Pc;
delta_so = sqrt(delta.^2+2*delta.*Qe+9*Qe.*Qe);
EvH = Ev0 - Pe -Qe;
EvL = Ev0 - Pe +0.5*(Qe - delta + delta_so);
Eso = Ev0 - Pe +0.5*(Qe - delta - delta_so);
end

function me = InGaAsP_me(w,v)
me_GaAs = 0.067;
me_GaP = 0.254;
me_InAs = 0.027;
me_InP = 0.077;
%  interp
% me = me_GaAs*(1-w)*v + me_GaP*(1-w)*(1-v)+me_InAs*w*v+me_InP*w*(1-v); %mo
% E. Herbert Li, Physica E (2000) 215-273
me = 0.0632*(1 - w)*v + 0.0213*w*v + 0.158*(1 - w)*(1 - v) + 0.077*w*(1 - v);
% Physics of opoelectronic devices shun lien chang
% me = 0.08 - 0.116*(1-w) + 0.026*v - 0.059*(1-w)*v + (0.064 - 0.02*(1-w))*v*v+...
%      (0.06+0.032*v)*(1-w)*(1-w);
end
function mHH = InGaAsP_mHH(w,v)
mHH_GaAs = 0.38;
mHH_GaP = 0.67;
mHH_InAs = 0.34;
mHH_InP = 0.61;

% mHH = mHH_GaAs*(1-w)*v + mHH_GaP*(1-w)*(1-v)+mHH_InAs*w*v+mHH_InP*w*(1-v); %mo
mHH = 0.5*(1 - w)*v + 0.5172*w*v + 0.54*(1 - w)*(1 - v) + 0.56*w*(1 - v);
end
function mLH = InGaAsP_mLH(w,v)
mLH_GaAs = 0.09;
mLH_GaP = 0.17;
mLH_InAs = 0.027;
mLH_InP = 0.12;

% mLH = mLH_GaAs*(1-w)*v + mLH_GaP*(1-w)*(1-v)+mLH_InAs*w*v+mLH_InP*w*(1-v); 
mLH = 0.088*(1 - w)*v + 0.024*w*v + 0.16*(1 - w)*(1 - v) + 0.12*w*(1 - v);
end

function mHHz = InGaAsP_mHHz(w,v)
gamma1= InGaAsP_gamma1(w,v);
gamma2= InGaAsP_gamma2(w,v);

mHHz = 1/(gamma1 - 2* gamma2);
end

function mHHt = InGaAsP_mHHt(w,v)
gamma1= InGaAsP_gamma1(w,v);
gamma2= InGaAsP_gamma2(w,v);

mHHt = 1/(gamma1 + gamma2);
end

function mLHz = InGaAsP_mLHz(w,v)
gamma1= InGaAsP_gamma1(w,v);
gamma2= InGaAsP_gamma2(w,v);

mLHz = 1/(gamma1 + 2* gamma2);
end

function mLHt = InGaAsP_mLHt(w,v)
gamma1= InGaAsP_gamma1(w,v);
gamma2= InGaAsP_gamma2(w,v);

mLHt = 1/(gamma1 - gamma2);
end

function epsr = InGaAsP_epsr(w,v)
epsr = 13.18*(1 - w)*v + 14.55*w*v + 11.1*(1 - w)*(1 - v) + 12.35*w*(1 - v);
end
function [taun, taup] = InGaAsP_tau(w,v)
taun_InAs = 3e-8;
taun_GaAs = 1e-9;
taun_AlAs = 1e-9;
taun = w*taun_InAs + (1-w-v)*taun_GaAs + v*taun_AlAs;

taup_InAs = 3e-6;
taup_GaAs = 2e-8;
taup_AlAs = 2e-8;
taup = w*taup_InAs + (1-w-v)*taup_GaAs + v*taup_AlAs;
end

function [n] = InGaAsP_n(w,v,wl)
%C. H. Henry, L. F. Johnson, R. A. Logan, and D. P. Clarke, 
%¡°Determination of the refractive index of InGaAsP epitaxial layers by mode line luminescence spectroscopy,¡± 
%IEEE J. Quantum Electron., vol. QE-21, pp. 1887-1892, 1985.
evnm =1239.842;% 1 eV = evnm/wl; evnm = hc/q*1e9
E = evnm./wl;
[Ec,Ev, EvH, EvL,Eso]=InGaAsP_band_strain(w,v);
E_PL = Ec-Ev;
A1 = 13.3510 - 5.4554*E_PL + 1.2332*E_PL.^2;
A2 = 0.7140 - 0.3606*E_PL;
E1 = 2.5048;
E2 = 0.1638;
n = sqrt(1+A1./(1-(E./(E_PL+E1)).^2)+A2./(1-(E./(E_PL+E2)).^2));
end