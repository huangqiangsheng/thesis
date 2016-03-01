%# Modeling of self_electro_optic_effect devices 1993 P.J.Mares and S.L.Chuang
%% E. Herbert Li, Physica E (2000) 215-273

function mat = AlGaAs_params(w,wl)
%strained Al_w_Ga_1-w_As
%wl: wavelength [nm]
if nargin <3
    wl=1550;%nm
end

mat.c11 = AlGaAs_c11(w);
mat.c12 = AlGaAs_c12(w);

mat.gamma1= AlGaAs_gamma1(w);
mat.gamma2= AlGaAs_gamma2(w);

mat.a = AlGaAs_a(w);
mat.b = AlGaAs_b(w);
mat.d = AlGaAs_d(w);
mat.ac = AlGaAs_ac(w);
mat.av = AlGaAs_av(w);

%% lattice length
mat.a0 = AlGaAs_a0(w);
%% band parameters ignoring strain
mat.delta = AlGaAs_delta(w);
mat.Evav0 = AlGaAs_Evav(w);
mat.Eg0 = AlGaAs_Eg(w);
[mat.Ec0, mat.Ev0, mat.Eso0]=AlGaAs_band(w);
%% band parameters with strain
[mat.Ec,mat.Ev, mat.EvH, mat.EvL,mat.Eso]=AlGaAs_band_strain(w);
mat.Eg = mat.Ec - mat.Ev;
%% strain
[mat.exx,mat.eyy,mat.ezz]=AlGaAs_strain(w);
%% effective mass [m0]
mat.me = AlGaAs_me(w);
mat.mHH = AlGaAs_mHH(w);
mat.mLH = AlGaAs_mLH(w);
mat.mHHz = AlGaAs_mHHz(w);
mat.mLHz = AlGaAs_mLHz(w);
mat.mHHt = AlGaAs_mHHt(w);
mat.mLHt = AlGaAs_mLHt(w);
%% low frequency relative epsilon
mat.epsr = AlGaAs_epsr(w);
%% lifetime [s]
% [mat.taun, mat.taup] = AlGaAs_tau(w);
%% n, refractive index with strain
%% n0, refractive index without considering strain
% [mat.n, mat.n0] = AlGaAs_n(w,wl);

end

function c11 = AlGaAs_c11(w)
    c11 = 11.90+0.12*w;
end
function c12 = AlGaAs_c12(w)
    c12 = 5.38 - 0.08*w;
end

function gamma1= AlGaAs_gamma1(w)
GaAs_gamma1 = 6.79;
AlAs_gamma1 = 3.45;
gamma1 = AlAs_gamma1*w + GaAs_gamma1*(1-w);
% gamma1 = w*gamma1_InAs + (1-w-v)*gamma1_GaAs + v*gamma1_AlAs;
end
function gamma2= AlGaAs_gamma2(w)
GaAs_gamma2 = 1.92;
AlAs_gamma2 = 0.68;
gamma2 = AlAs_gamma2*w + GaAs_gamma2*(1-w);
end

function a = AlGaAs_a(w)
a_GaAs = -8.33;
a_AlAs = -8.11;

a = a_AlAs * w + a_GaAs*(1-w);

end

function b = AlGaAs_b(w)
b = -1.7 + 0.2*w;
end

function d = AlGaAs_d(w)
d_GaAs = -4.55;
d_AlAs = -3.4;

d = d_AlAs * w + d_GaAs*(1-w);
end
function a0 = AlGaAs_a0(w)
a0 = 5.6533 + 0.0078*w;
end
function Eg = AlGaAs_Eg(w)
Eg = 1.426+1.247*w;
end
function ac = AlGaAs_ac(w)
ac_GaAs = -7.17;
ac_AlAs = -5.64;

ac =  ac_AlAs * w + ac_GaAs*(1-w);
end
function av = AlGaAs_av(w)
av_GaAs = 1.16;
av_AlAs = 2.46;

av = av_AlAs * w + av_GaAs*(1-w);
end

function delta = AlGaAs_delta(w)
delta_GaAs = 0.34;
delta_AlAs = 0.28;
delta =delta_AlAs * w + delta_GaAs*(1-w);
end

function Evav = AlGaAs_Evav(w)
Evav_GaAs = -6.92;
Evav_AlAs = -7.49;

Evav = Evav_AlAs * w + Evav_GaAs*(1-w);
end
function [Ec, Ev, Eso]=AlGaAs_band(w)
Evav = AlGaAs_Evav(w);
delta = AlGaAs_delta(w);
Eg = AlGaAs_Eg(w);

Ev = Evav + delta/3;
Eso = Ev - delta;
Ec = Ev + Eg;
end

function [exx,eyy,ezz]=AlGaAs_strain(w)
a0_InP = 5.8688;%A
a0 = AlGaAs_a0(w);
exx = (a0_InP-a0)/a0;
eyy = exx;
c11 = AlGaAs_c11(w);
c12 = AlGaAs_c12(w);
ezz = -2*c12/c11*exx;
end

function [Ec,Ev, EvH, EvL,Eso]=AlGaAs_band_strain(w)
[Ec0, Ev0, Eso0]=AlGaAs_band(w);
Evav0 = AlGaAs_Evav(w);
[exx,eyy,ezz]=AlGaAs_strain(w);
delta = AlGaAs_delta(w);

Eg = AlGaAs_Eg(w);
av = AlGaAs_av(w);
b = AlGaAs_b(w);
ac = AlGaAs_ac(w);

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

function me = AlGaAs_me(w)
me_GaAs = 0.0665;
me_AlAs = 0.15;

me =me_AlAs * w + me_GaAs*(1-w);
end
function mHH = AlGaAs_mHH(w)
mHH_GaAs = 0.38;
mHH_AlAs = 0.79;

mHH = mHH_AlAs * w + mHH_GaAs*(1-w);
end
function mLH = AlGaAs_mLH(w)
mLH_GaAs = 0.087;
mLH_AlAs = 0.12;

mLH = mLH_AlAs * w + mLH_GaAs*(1-w);
end

function mHHz = AlGaAs_mHHz(w)
gamma1= AlGaAs_gamma1(w);
gamma2= AlGaAs_gamma2(w);

mHHz = 1/(gamma1 - 2* gamma2);
end

function mHHt = AlGaAs_mHHt(w)
gamma1= AlGaAs_gamma1(w);
gamma2= AlGaAs_gamma2(w);

mHHt = 1/(gamma1 + gamma2);
end

function mLHz = AlGaAs_mLHz(w)
gamma1= AlGaAs_gamma1(w);
gamma2= AlGaAs_gamma2(w);

mLHz = 1/(gamma1 + 2* gamma2);
end

function mLHt = AlGaAs_mLHt(w)
gamma1= AlGaAs_gamma1(w);
gamma2= AlGaAs_gamma2(w);

mLHt = 1/(gamma1 - gamma2);
end

function epsr = AlGaAs_epsr(w)
epsr = 13.18 - 3.12*w;
end

% function [taun, taup] = AlGaAs_tau(w)
% taun_InAs = 3e-8;
% taun_GaAs = 1e-9;
% taun_AlAs = 1e-9;
% taun = w*taun_InAs + (1-w-v)*taun_GaAs + v*taun_AlAs;
% 
% taup_InAs = 3e-6;
% taup_GaAs = 2e-8;
% taup_AlAs = 2e-8;
% taup = w*taup_InAs + (1-w-v)*taup_GaAs + v*taup_AlAs;
% end
% 
% function [nst, nlm] = AlGaAs_n(w,wl)
% %%%%(* S. Adachi, JAP, vol .53, no .8 (1982) pp .5863-5869. *)
% %%%%(* Adachi's Parameter*)
% %%%% wl in nm
% H =@(z) z>=0;
% %AA =9.29 *w+25.3 *v+5.14 *(1-w-v);
% %BB =7.86 *w-0.8 *v+10.15  *(1-w-v);
% % from atlas manual 2010
% AA =5.14 *w+ 25.3 *v+ 6.3 *(1-w-v);
% BB =10.15 *w-0.8 *v+9.4  *(1-w-v);
% 
% f =@(x) x.^-2*(2-(1+x).^0.5-((1-x)>=0)*(1-x).^0.5);%% chi -> x
% 
% evnm =1239.842;% 1 eV = evnm/wl; evnm = hc/q*1e9
% wl_eV = evnm/wl;
% [Ec, Ev, Eso]=AlGaAs_band(w);
% Chi0 =  wl_eV/(Ec-Ev);
% Chiso = wl_eV/(Ec-Eso);
% nlm = (AA.*(f(Chi0)+0.5*(Chiso./Chi0).^1.5*f(Chiso))+BB).^0.5;
% 
% [Ec,Ev, EvH, EvL,Eso]=AlGaAs_band_strain(w);
% 
% Chi1 = wl_eV/(Ec-EvH);
% Chi2 = wl_eV/(Ec-EvL);
% Chi3 = wl_eV/(Ec-Eso);
% 
% %%%%(* P.K. Bhattacharya et al., Solid State Electron., vol .29, no .2 (1986) pp. 261-267.*)
% nst = (0.5*AA.*(((Ec-Ev)./(Ec-EvH)).^1.5.*f(Chi1)+((Ec-Ev)./(Ec-EvL)).^1.5*f(Chi2)+((Ec-Ev)./(Ec-Eso)).^1.5*f(Chi3))+BB).^0.5;
% end