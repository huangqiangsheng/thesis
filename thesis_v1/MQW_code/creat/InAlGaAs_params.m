%# calculate the parameters for InAlGaAs
%% E. Herbert Li, Physica E (2000) 215-273

function mat = InAlGaAs_params(w,v, wl)
%strained In_w Al_v Ga_1-w-v As
%wl: wavelength [nm]
if nargin <3
    wl=1550;%nm
end

mat.c11 = InAlGaAs_c11(w,v);
mat.c12 = InAlGaAs_c12(w,v);
mat.c44 = InAlGaAs_c44(w,v);
mat.gamma1= InAlGaAs_gamma1(w,v);
mat.gamma2= InAlGaAs_gamma2(w,v);
mat.gamma3= InAlGaAs_gamma3(w,v);
mat.a = InAlGaAs_a(w,v);
mat.b = InAlGaAs_b(w,v);
mat.d = InAlGaAs_d(w,v);
mat.ac = InAlGaAs_ac(w,v);
mat.av = InAlGaAs_av(w,v);

%% lattice length
mat.a0 = InAlGaAs_a0(w,v);
%% band parameters ignoring strain
mat.delta = InAlGaAs_delta(w,v);
mat.Evav0 = InAlGaAs_Evav(w,v);
mat.Eg0 = InAlGaAs_Eg(w,v);
[mat.Ec0, mat.Ev0, mat.Eso0]=InAlGaAs_band(w,v);
%% band parameters with strain
[mat.Ec,mat.Ev, mat.EvH, mat.EvL,mat.Eso]=InAlGaAs_band_strain(w,v);
mat.Eg = mat.Ec - mat.Ev;
%% strain
[mat.exx,mat.eyy,mat.ezz]=InAlGaAs_strain(w,v);
%% effective mass [m0]
mat.me = InAlGaAs_me(w,v);
mat.mHH = InAlGaAs_mHH(w,v);
mat.mLH = InAlGaAs_mLH(w,v);
mat.mHHz = InAlGaAs_mHHz(w,v);
mat.mLHz = InAlGaAs_mLHz(w,v);
mat.mHHt = InAlGaAs_mHHt(w,v);
mat.mLHt = InAlGaAs_mLHt(w,v);
%% low frequency relative epsilon
mat.epsr = InAlGaAs_epsr(w,v);
%% lifetime [s]
[mat.taun, mat.taup] = InAlGaAs_tau(w,v);
%% n, refractive index with strain
%% n0, refractive index without considering strain
[mat.n, mat.n0] = InAlGaAs_n(w,v,wl);

end

function c11 = InAlGaAs_c11(w,v)
c11_InAs = 8.329;
c11_GaAs = 11.879; %11.9;
c11_AlAs = 12.02;

% c11 = 11.9-3.571*w+0.12*v;
c11 = c11_InAs*w + c11_AlAs*v + (1-w-v)*c11_GaAs;
end
function c12 = InAlGaAs_c12(w,v)
c12_InAs = 4.526;
c12_GaAs = 5.376;
c12_AlAs = 5.34;

% c12 = 5.38-0.854*w-0.08*v;
c12 = c12_InAs*w + c12_AlAs*v + (1-w-v)*c12_GaAs;
end

function c44 = InAlGaAs_c44(w,v)
c44_InAs = 3.96;
c44_GaAs = 5.94;
c44_AlAs = 5.42;

c44 = w*c44_InAs+v*c44_AlAs+(1-w-v)*c44_GaAs;
end
function gamma1= InAlGaAs_gamma1(w,v)
gamma1_InAs = 20.4;
gamma1_GaAs = 6.85;
gamma1_AlAs = 3.45;

gamma1 = 6.68+12.22*w-2.52*v;
% gamma1 = w*gamma1_InAs + (1-w-v)*gamma1_GaAs + v*gamma1_AlAs;
end
function gamma2= InAlGaAs_gamma2(w,v)
gamma2_InAs = 8.3;
gamma2_GaAs = 2.1;
gamma2_AlAs = 0.68;

gamma2 = 2.34 + 6.03*w-0.97*v;
end
function gamma3= InAlGaAs_gamma3(w,v)
gamma3_InAs = 9.1;
gamma3_GaAs = 2.9;
gamma3_AlAs = 1.29;

gamma3 = w*gamma3_InAs + (1-w-v)*gamma3_GaAs + v*gamma3_AlAs;
end
function a = InAlGaAs_a(w,v)
a_InAs = -6.08;
a_GaAs = -8.33;
a_AlAs = -8.11;

a = w*a_InAs + (1-w-v)*a_GaAs + v*a_AlAs;
%a =-1/3*(InAlGaAs_c11(w,v)+2*InAlGaAs_c12(w,v))
end
function b = InAlGaAs_b(w,v)
b_InAs = -1.8;
b_GaAs = -1.7;
b_AlAs = -1.5;

b = w*b_InAs + (1-w-v)*b_GaAs + v*b_AlAs;
% b = -1.7-0.1 *w+0.2 *v;
end
function d = InAlGaAs_d(w,v)
d_InAs = -3.6;
d_GaAs = -4.55;
d_AlAs = -3.4;

d = w*d_InAs + (1-w-v)*d_GaAs + v*d_AlAs;
end
function a0 = InAlGaAs_a0(w,v)
a0_InAs = -6.0583;
a0_GaAs = -5.6533;
a0_AlAs = -5.6611;

%a0 = w*a0_InAs + (1-w-v)*a0_GaAs + v*a0_AlAs;
a0 = 5.653325+0.404975*w+0.007775*v;
end
function Eg = InAlGaAs_Eg(w,v)
Eg_InAs = 0.354;
Eg_GaAs = 1.424;
Eg_AlAs = 3.017;

Eg = 0.36+2.093*v+0.629*(1-w-v)+0.577*v.*v+0.436*(1-w-v).*(1-w-v)+1.013*v.*(1-w-v)-2*w.*v.*(1-w-v);
end
function ac = InAlGaAs_ac(w,v)
ac_InAs = -5.08;
ac_GaAs = -7.17;
ac_AlAs = -5.64;

ac = w*ac_InAs + (1-w-v)*ac_GaAs + v*ac_AlAs;
end
function av = InAlGaAs_av(w,v)
av_InAs = 1;
av_GaAs = 1.16;
av_AlAs = 2.47;

av = w*av_InAs + (1-w-v)*av_GaAs + v*av_AlAs;
end

function delta = InAlGaAs_delta(w,v)
delta_InAs = 0.38;%eV
delta_GaAs = 0.34;
delta_AlAs = 0.28;

% delta = w*0.41 + (1-w-v)*0.34 + v*0.275;
delta = w*delta_InAs + v * delta_AlAs + (1-w-v)*delta_GaAs;
end

function Evav = InAlGaAs_Evav(w,v)
Evav_InAs = -6.67;%eV
Evav_GaAs = -6.92;
Evav_AlAs = -7.49;

Evav = w*Evav_InAs + (1-w-v)*Evav_GaAs + v*Evav_AlAs;
end
function [Ec, Ev, Eso]=InAlGaAs_band(w,v)
Evav = InAlGaAs_Evav(w,v);
delta = InAlGaAs_delta(w,v);
Eg = InAlGaAs_Eg(w,v);

Ev = Evav + delta/3;
Eso = Ev - delta;
Ec = Ev + Eg;
end

function [exx,eyy,ezz]=InAlGaAs_strain(w,v)
a0_InP = 5.8688;%A
a0 = InAlGaAs_a0(w,v);
exx = (a0_InP-a0)/a0;
eyy = exx;
c11 = InAlGaAs_c11(w,v);
c12 = InAlGaAs_c12(w,v);
ezz = -2*c12/c11*exx;
end

function [Ec,Ev, EvH, EvL,Eso]=InAlGaAs_band_strain(w,v)
[Ec0, Ev0, Eso0]=InAlGaAs_band(w,v);
Evav0 = InAlGaAs_Evav(w,v);
[exx,eyy,ezz]=InAlGaAs_strain(w,v);
delta = InAlGaAs_delta(w,v);

Eg = InAlGaAs_Eg(w,v);
av = InAlGaAs_av(w,v);
b = InAlGaAs_b(w,v);
ac = InAlGaAs_ac(w,v);

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

function me = InAlGaAs_me(w,v)
me_InAs = 0.0213; %m0
me_GaAs = 0.0632;
me_AlAs = 0.1719;
me = w*me_InAs + (1-w-v)*me_GaAs + v*me_AlAs;
end
function mHH = InAlGaAs_mHH(w,v)
mHH_InAs = 0.5172;%0.4; %m0
mHH_GaAs = 0.5;
mHH_AlAs = 0.7;
mHH = w*mHH_InAs + (1-w-v)*mHH_GaAs + v*mHH_AlAs;
end
function mLH = InAlGaAs_mLH(w,v)
mLH_InAs = 0.024;%0.026; %m0
mLH_GaAs = 0.089;
mLH_AlAs = 0.1415;
mLH = w*mLH_InAs + (1-w-v)*mLH_GaAs + v*mLH_AlAs;
end

function mHHz = InAlGaAs_mHHz(w,v)
gamma1= InAlGaAs_gamma1(w,v);
gamma2= InAlGaAs_gamma2(w,v);

mHHz = 1/(gamma1 - 2* gamma2);
end

function mHHt = InAlGaAs_mHHt(w,v)
gamma1= InAlGaAs_gamma1(w,v);
gamma2= InAlGaAs_gamma2(w,v);

mHHt = 1/(gamma1 + gamma2);
end

function mLHz = InAlGaAs_mLHz(w,v)
gamma1= InAlGaAs_gamma1(w,v);
gamma2= InAlGaAs_gamma2(w,v);

mLHz = 1/(gamma1 + 2* gamma2);
end

function mLHt = InAlGaAs_mLHt(w,v)
gamma1= InAlGaAs_gamma1(w,v);
gamma2= InAlGaAs_gamma2(w,v);

mLHt = 1/(gamma1 - gamma2);
end

function epsr = InAlGaAs_epsr(w,v)
epsr_InAs = 15.15;
epsr_GaAs = 13.18;
epsr_AlAs = 10.06;

epsr = 13.18 + 1.37*w - 3.12*v;
end
function [taun, taup] = InAlGaAs_tau(w,v)
taun_InAs = 3e-8;
taun_GaAs = 1e-9;
taun_AlAs = 1e-9;
taun = w*taun_InAs + (1-w-v)*taun_GaAs + v*taun_AlAs;

taup_InAs = 3e-6;
taup_GaAs = 2e-8;
taup_AlAs = 2e-8;
taup = w*taup_InAs + (1-w-v)*taup_GaAs + v*taup_AlAs;
end

function [nst, nlm] = InAlGaAs_n(w,v,wl)
%%%%(* S. Adachi, JAP, vol .53, no .8 (1982) pp .5863-5869. *)
%%%%(* Adachi's Parameter*)
%%%% wl in nm
H =@(z) z>=0;
%AA =9.29 *w+25.3 *v+5.14 *(1-w-v);
%BB =7.86 *w-0.8 *v+10.15  *(1-w-v);
% from atlas manual 2010
AA =5.14 *w+ 25.3 *v+ 6.3 *(1-w-v);
BB =10.15 *w-0.8 *v+9.4  *(1-w-v);

f =@(x) x.^-2*(2-(1+x).^0.5-((1-x)>=0)*(1-x).^0.5);%% chi -> x

evnm =1239.842;% 1 eV = evnm/wl; evnm = hc/q*1e9
wl_eV = evnm/wl;
[Ec, Ev, Eso]=InAlGaAs_band(w,v);
Chi0 =  wl_eV/(Ec-Ev);
Chiso = wl_eV/(Ec-Eso);
nlm = (AA.*(f(Chi0)+0.5*(Chiso./Chi0).^1.5*f(Chiso))+BB).^0.5;

[Ec,Ev, EvH, EvL,Eso]=InAlGaAs_band_strain(w,v);

Chi1 = wl_eV/(Ec-EvH);
Chi2 = wl_eV/(Ec-EvL);
Chi3 = wl_eV/(Ec-Eso);

%%%%(* P.K. Bhattacharya et al., Solid State Electron., vol .29, no .2 (1986) pp. 261-267.*)
nst = (0.5*AA.*(((Ec-Ev)./(Ec-EvH)).^1.5.*f(Chi1)+((Ec-Ev)./(Ec-EvL)).^1.5*f(Chi2)+((Ec-Ev)./(Ec-Eso)).^1.5*f(Chi3))+BB).^0.5;
end