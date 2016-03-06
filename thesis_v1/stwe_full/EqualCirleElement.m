function [Rc,Lm,Rsm,Cim,Ce] = EqualCirleElement...
         (freq,Zc,Gamma,Wa,Wg,delta_InGaAs,delta_p,delta_n,delta_sub_n,...
         hInGaAs,hp,hn,d)
% freq is the frequency
% Zc is the characteristic WG impendence
% Gamma is the propagation constant
% Wa is the width of waveguide
% Wg is the distance between the ground
% delta_p is the doping of the p_InP
% delta_n is the doping of the mesa_n_InP
% delta_sub_n is the doping of the sub_n_InP
% hp is the thickness of the p_InP
% hn is the thickness of the n_InP
% d is the thickness of the sub_n_InP
% reference:PhD thesis Stefan Irmscher, Appendix A
%ohmn contact resistor
Rohmn = 0.8*1e-3;
Rsm = Rohmn+hp./(Wa.*delta_p)+hn./(Wa.*delta_n)+Wg./(4*d.*delta_sub_n)+hInGaAs./(Wa.*delta_InGaAs);
% Rsm = 3e-3;
% reference; 2003,Robert Lewen,
% "High-Speed Electroabsorption Modulators and p-i-n Photodiodes for 
%  Fiber-Optic Communications"

omega = 2*pi*freq;
Zs = Zc.*Gamma;
Rc = real(Zs);
Lm = imag(Zs)./omega;
Yp = Gamma./Zc;
% Cim = sqrt(real(Yp)./((1-real(Yp).*Rsm).*omega.^2.*Rsm));
% Ce = imag(Yp)./omega- real(Yp)./(omega.^2.*Rsm.*Cim);
Yp = Zc./Gamma;
Cim = -1./(omega.*imag(Yp));
Ce = imag(1./Yp)./omega- Cim./(1+omega.^2.*Rsm.^2.*Cim.^2);

end