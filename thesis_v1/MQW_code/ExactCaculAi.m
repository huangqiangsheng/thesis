function [Eout waveFun]=ExactCaculAi(V0,Lw,Lb,mw,mb,F)
% reference "Exact calculation of quasibound states of an isolated QM with uniform electric field:Quatum-well Stark resonace"
% HQS 2011_11_28    
% V0(eV) Lw,Lb(m) mw,mb(m0),F(V/m)
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    V0 = V0*q;
    mw = mw*m0;
    mb = mb*m0;
    E0_w = (hbar*pi/Lw)^2/(2*mw);
    E0_b = (hbar*pi/Lw)^2/(2*mb);
    Fbar_w = q*F*Lw/E0_w;
    Fbar_b = q*F*Lw/E0_b;
    V0bar_b = V0/E0_b;
    E = linspace(-F*(Lw/2+Lb)*q-abs(V0)/q,abs(V0)/q+F*(Lw/2+Lb)*q,1e3);
    fval = zeros(size(E));
    for iter = 1:length(E)
        fval(iter) = BoundAi([E(iter) 0],Fbar_w,E0_w,Fbar_b,V0bar_b,E0_b);
    end
    figure(100);
    plot(E,log(fval)); xlabel('E(eV)');
    hold on;
    [gx,gy] = ginput(1);
    plot([gx(1) gx(1)],[0 gy(1)],'r');
    [Eout,fout] = fminsearch(@(E)BoundAi(E,Fbar_w,E0_w,Fbar_b,V0bar_b,E0_b),[gx(1) 0],optimset('TolX',1e-32));
    plot(Eout(1),log(fout),'r*');hold off
%     fprintf('E = %.6e - i%.3e/2\nfval=%.1e + i%.1e\n',Eout(1),Eout(2),real(fout),imag(fout));
    Eout = Eout(1) - 1i*Eout(2)/2;
    Ebar_w = Eout * q / E0_w;
    Ebar_b = Eout * q / E0_b;
    waveFun = DrawAi(Ebar_w,Fbar_w,Ebar_b,Fbar_b,V0bar_b,E0_b,mw,mb,Lw,Lb);
end

function fval = BoundAi(Ein,Fbar_w,E0_w,Fbar_b,V0bar_b,E0_b)
    q = 1.6021765e-19;  %C
    Ebar_w = Ein*q/E0_w;
    Ebar_b = Ein*q/E0_b;
    E_w = Ebar_w(1) - 1i * Ebar_w(2)/2;
    E_b = Ebar_b(1) - 1i * Ebar_b(2)/2;
    n1_p = -1*(pi/Fbar_w)^(2/3)*(E_w - 1/2 * Fbar_w);
    n1_n = -1*(pi/Fbar_w)^(2/3)*(E_w + 1/2 * Fbar_w);
    n2_p = -1*(pi/Fbar_b)^(2/3)*(E_b - V0bar_b - 1/2 * Fbar_b);
    n2_n = -1*(pi/Fbar_b)^(2/3)*(E_b - V0bar_b + 1/2 * Fbar_b);
    AiMatrix = [ airy(0,n1_p) airy(2,n1_p) -airy(0,n2_p) 0;...
                 airy(1,n1_p) airy(3,n1_p) -airy(1,n2_p) 0;...
                 airy(0,n1_n) airy(2,n1_n) 0 -(airy(2,n2_n)+1i*airy(0,n2_n));...
                 airy(1,n1_n) airy(3,n1_n) 0 -(airy(3,n2_n)+1i*airy(1,n2_n))];

    fval = abs(det(AiMatrix));
end

function [waveFun] = DrawAi(Ebar_w,Fbar_w,Ebar_b,Fbar_b,V0bar_b,E0_b,mw,mb,Lw,Lb)
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    n1_p = -1*(pi/Fbar_w)^(2/3)*(Ebar_w - 1/2 * Fbar_w);
    n1_n = -1*(pi/Fbar_w)^(2/3)*(Ebar_w + 1/2 * Fbar_w);
    n2_p = -1*(pi/Fbar_b)^(2/3)*(Ebar_b - V0bar_b - 1/2 * Fbar_b);
    n2_n = -1*(pi/Fbar_b)^(2/3)*(Ebar_b - V0bar_b + 1/2 * Fbar_b);
    fprintf('Figure = %.3e\n',real(n2_n));
    AiMatrix = [ airy(0,n1_p) airy(2,n1_p) -airy(0,n2_p) 0;...
                 airy(1,n1_p) airy(3,n1_p) -airy(1,n2_p) 0;...
                 airy(0,n1_n) airy(2,n1_n) 0 -(airy(2,n2_n)+1i*airy(0,n2_n));...
                 airy(1,n1_n) airy(3,n1_n) 0 -(airy(3,n2_n)+1i*airy(1,n2_n))];
    V = null(AiMatrix);
    a0 = V(1);
    b0 = V(2);
    a2 = V(3);
    a1 = V(4);
    E = Ebar_b*E0_b;
    F = Fbar_b*E0_b/(Lw*q);
    V0 = V0bar_b*E0_b;
    sample_num = 1e3;
    dz = Lw/sample_num;
    z1 = -Lw/2-Lb : dz: -Lw/2;
    z2 = -Lw/2 : dz : Lw/2;
    z3 = Lw/2 : dz : Lw/2 + Lb;
    n2_n = -(2*mb/(q*hbar*F)^2)^(1/3)*(E-V0-q*F*z1);
    n1 = -(2*mw/(q*hbar*F)^2)^(1/3)*(E-q*F*z2);
    n2_p = -(2*mb/(q*hbar*F)^2)^(1/3)*(E-V0-q*F*z3);
    waveFun = zeros(2,length(z1)+length(z2)+length(z3));
    waveFun(2,1:length(z1)) = a1*(airy(2,n2_n) + 1i*airy(0,n2_n));
    waveFun(2,length(z1)+1:length(z2)+length(z1)) = a0*airy(0,n1)+b0*airy(2,n1);
    waveFun(2,1+length(z1)+length(z2):end) = a2*airy(0,n2_p);
    waveFun(2,:) = waveFun(2,:)./sqrt(sum(abs(waveFun(2,:)).^2) * dz);
    waveFun(1,:) = [z1,z2,z3];
end
