function mat = InP_params(wl)
    %strained In_w Ga_1-w As_v P_1-v
%wl: wavelength [nm]
%flag:1 well 0 barrier
if nargin <1
    wl=1550;%nm
end

mat.c11 = InP_c11();
mat.c12 = InP_c12();

mat.gamma1= InP_gamma1();
mat.gamma2= InP_gamma2();

mat.a = InP_a();
mat.b = InP_b();
mat.d = InP_d();
mat.ac = InP_ac();
mat.av = InP_av();

%% lattice length
mat.a0 = InP_a0();
%% band parameters ignoring strain
mat.delta = InP_delta();
mat.Evav0 = InP_Evav();
mat.Eg0 = InP_Eg();
[mat.Ec0, mat.Ev0, mat.Eso0]=InP_band();
%% band parameters with strain
[mat.Ec,mat.Ev, mat.EvH, mat.EvL,mat.Eso]=InP_band_strain();
mat.Eg = mat.Ec - mat.Ev;
%% strain
[mat.exx,mat.eyy,mat.ezz]=InP_strain();
%% effective mass [m0]
mat.me = InP_me();
mat.mHH = InP_mHH();
mat.mLH = InP_mLH();
mat.mHHz = InP_mHHz();
mat.mLHz = InP_mLHz();
mat.mHHt = InP_mHHt();
mat.mLHt = InP_mLHt();
%% low frequency relative epsilon
mat.epsr = InP_epsr();
%% lifetime [s]
% [mat.taun, mat.taup] = InP_tau(w,v);
%% n, refractive index with strain
%% n0, refractive index without considering strain
% [mat.n, mat.n0] = InP_n(w,v,wl);

end

function c11 = InP_c11()
    c11 = 10.11;
end
function c12 = InP_c12()
    c12 = 5.61;
end

function gamma1= InP_gamma1()
gamma1 = 4.95;
end
function gamma2= InP_gamma2()
gamma2 = 1.65;
end

function a = InP_a()
a = -6.31;
end

function b = InP_b()
b = -1.7;
end

function d = InP_d()
d = -5.6;
end
function a0 = InP_a0()
a0 = 5.8688;
end
function Eg = InP_Eg()
Eg = 1.344;
end
function ac = InP_ac()
ac= -5.04;
end
function av = InP_av()
av = 1.27;
end

function delta = InP_delta()
delta = 0.10;
end

function Evav = InP_Evav()
Evav = -7.04;
end
function [Ec, Ev, Eso]=InP_band()
Evav = InP_Evav();
delta = InP_delta();
Eg = InP_Eg();

Ev = Evav + delta/3;
Eso = Ev - delta;
Ec = Ev + Eg;
end

function [exx,eyy,ezz]=InP_strain()
a0_InP = 5.8688;%A
a0 = InP_a0();
exx = (a0_InP-a0)/a0;
eyy = exx;
c11 = InP_c11();
c12 = InP_c12();
ezz = -2*c12/c11*exx;
end

function [Ec,Ev, EvH, EvL,Eso]=InP_band_strain()
[Ec0, Ev0, Eso0]=InP_band();
Evav0 = InP_Evav();
[exx,eyy,ezz]=InP_strain();
delta = InP_delta();

Eg = InP_Eg();
av = InP_av();
b = InP_b();
ac = InP_ac();

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

function me = InP_me()
me= 0.077;
end
function mHH = InP_mHH()
mHH = 0.61;
end
function mLH = InP_mLH()
mLH = 0.12;
end

function mHHz = InP_mHHz()
gamma1= InP_gamma1();
gamma2= InP_gamma2();

mHHz = 1/(gamma1 - 2* gamma2);
end

function mHHt = InP_mHHt()
gamma1= InP_gamma1();
gamma2= InP_gamma2();

mHHt = 1/(gamma1 + gamma2);
end

function mLHz = InP_mLHz()
gamma1= InP_gamma1();
gamma2= InP_gamma2();

mLHz = 1/(gamma1 + 2* gamma2);
end

function mLHt = InP_mLHt()
gamma1= InP_gamma1();
gamma2= InP_gamma2();

mLHt = 1/(gamma1 - gamma2);
end

function epsr = InP_epsr()
epsr = 12.56;
end
% function [taun, taup] = InP_tau()
% taup = w*taup_InAs + (1-w-v)*taup_GaAs + v*taup_AlAs;
% end
% 
% function [nst, nlm] = InP_n(wl)
% %%%%(* S. Adachi, JAP, vol .53, no .8 (1982) pp .5863-5869. *)
% %%%%(* Adachi's Parameter*)
% %%%% wl in nm
% 
% %%%%(* P.K. Bhattacharya et al., Solid State Electron., vol .29, no .2 (1986) pp. 261-267.*)
% end
