function [waveFun, Energy,flag] = TMM_EXCAL(x, v, m, pic,tol,V0,Lw,Lb,mw,mb,F,Flag)
%% 20111124 by HQS
%% Flag: 0 electron 1 hole
%% modified by HQS 2012_1_16 correct Ai Function  
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    [waveFun, Energy] = QW_FD1d(x, v, m);
    Energy = Energy *q;
    
    flag = 3;
  
    E = linspace(min(v)*1.2,max(v)*0.8,500);
 
    f22 = ones(1,length(E));
    for iter = 1:length(E)
        f22(iter) = F_E(E(iter),x, v, m);
    end
    figure(pic);plot(E,f22);
    xlabel('E(eV)');hold on;
    turnp = turnpoint(f22);
    dE = E(2)-E(1);
    gx = [E(turnp(1)-10), E(turnp(1)+10)]; 
    gy = [f22(turnp(1)-10), f22(turnp(1)+10)];  

%     [gx,gy] = ginput(2);
%     gx =0.0;
    plot([gx(1) gx(1)],[0 gy(1)],'r');
    plot([gx(2) gx(2)],[0 gy(2)],'r');
%     [E0,fval,exitflag] = fminsearch(@(E)F_E(E,x, v, m),gx(1), optimset('TolX',tol));
    [E0,fval,exitflag] = fminbnd(@(E)F_E(E,x, v, m),gx(1),gx(2),optimset('TolX',tol));
    plot(E0(1),log(abs(fval)+1),'r*');
    if fval > 5
    % TMM cannot find a reasonable solution then adopt Ai function
        fprintf('Ai Function\n');
        V0 = V0*q;
        mw = mw*m0;
        mb = mb*m0;
        E0_w = (hbar*pi/Lw)^2/(2*mw);
        E0_b = (hbar*pi/Lw)^2/(2*mb);
        Fbar_w = q*F*Lw/E0_w;
        Fbar_b = q*F*Lw/E0_b;
        V0bar_b = V0/E0_b;
        [Eout,fout] = fminsearch(@(E)BoundAi(E,Fbar_w,E0_w,Fbar_b,V0bar_b,E0_b,mw,mb),[gx(1) 0],optimset('TolX',tol));
        plot(Eout(1),log(abs(fout)+1),'y*');
        Eout = Eout(1) - 1i*Eout(2)/2;
        Ebar_w = Eout * q / E0_w;
        Ebar_b = Eout * q / E0_b;
        waveFun = DrawAi(Ebar_w,Fbar_w,Ebar_b,Fbar_b,V0bar_b,E0_b,mw,mb,Lw,Lb,x);
        if Flag == 1
            waveFun = waveFun(end:-1:1);
        end
        Energy = Eout;
    else
        fprintf('TMM Method\n');
        Energy = E0;
        waveFun = DrawWav(Energy,x,v,m);
    end
    hold off;
end

function waveFun = DrawWav(Ein,x,v,m)
%% 20111124 by HQS
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    v = v * q;
    m = m * m0;
    E = Ein(1);
%     E = Ein(1) - 1i*Ein(2)/2;
    E = E * q;
    dh = x(2) - x(1);
    k = sqrt(2*m./(hbar^2).*(E - v));
    P = m(2:end).*k(1:end-1)./(m(1:end-1).*k(2:end));
    A = ones(1,length(x));
    B = ones(1,length(x));
    A(1) = 0;
    B(1) = 1;
    F = [1 0; 0 1];
    for iter = 1 : (length(x) -1)
        F = 1/2*[(1+P(iter))*exp(1i*k(iter+1)*dh) (1-P(iter))*exp(1i*k(iter+1)*dh);...
                 (1-P(iter))*exp(-1i*k(iter+1)*dh) (1+P(iter))*exp(-1i*k(iter+1)*dh)]*...
                 F;
        C = F * [A(1);B(1)];
        A(iter+1) = C(1);
        B(iter+1) = C(2);
    end
    waveFun = A + B;
    waveFun = waveFun/sqrt((sum(abs(waveFun).^2))*dh); %归一化
end
function f22 = F_E(Ein,x, v, m)
%% 20111124 by HQS
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    v = v * q;
    m = m * m0;
    E = Ein;
%     E = Ein(1) - 1i*Ein(2)/2;
    E = E * q;
    dh = x(2) - x(1);
    
    k= sqrt(2*m./(hbar^2).*(E - v));
    
    P = m(2:end).*k(1:end-1)./(m(1:end-1).*k(2:end));
    
    F = [1,0;0,1];
    
    for iter = 1:(length(x)-1)
        F = 1/2*[(1+P(iter))*exp(1i*k(iter+1)*dh) (1-P(iter))*exp(1i*k(iter+1)*dh);...
             (1-P(iter))*exp(-1i*k(iter+1)*dh) (1+P(iter))*exp(-1i*k(iter+1)*dh)]*...
              F;
    end
    f22 = log(abs(F(2,2))+1);
end
function fval = BoundAi(Ein,Fbar_w,E0_w,Fbar_b,V0bar_b,E0_b,mw,mb)
%% modified by HQS 2012_1_16 correct AiMatrix
   %寻找最小值
    q = 1.6021765e-19;  %C
    Ebar_w = Ein*q/E0_w;
    Ebar_b = Ein*q/E0_b;
    E_w = Ebar_w(1) - 1i * Ebar_w(2)/2;
    E_b = Ebar_b(1) - 1i * Ebar_b(2)/2;
    n1_p = -1*(pi/Fbar_w)^(2/3)*(E_w - 1/2 * Fbar_w);
    n1_n = -1*(pi/Fbar_w)^(2/3)*(E_w + 1/2 * Fbar_w);
    n2_p = -1*(pi/Fbar_b)^(2/3)*(E_b - V0bar_b - 1/2 * Fbar_b);
    n2_n = -1*(pi/Fbar_b)^(2/3)*(E_b - V0bar_b + 1/2 * Fbar_b);
    AiMatrix = [ airy(0,n1_p), airy(2,n1_p), -airy(0,n2_p), 0;...
                 airy(1,n1_p)*(mb^(2/3)), airy(3,n1_p)*(mb^(2/3)), -airy(1,n2_p)*(mw^(2/3)), 0;...
                 airy(0,n1_n), airy(2,n1_n), 0, -(airy(2,n2_n)+1i*airy(0,n2_n));...
                 airy(1,n1_n)*(mb^(2/3)), airy(3,n1_n)*(mb^(2/3)), 0, -(airy(3,n2_n)+1i*airy(1,n2_n))*(mw^(2/3))];
% the old type is error because the bondary condition is wrong    
% AiMatrix = [ airy(0,n1_p) airy(2,n1_p) -airy(0,n2_p) 0;...
%                  airy(1,n1_p) airy(3,n1_p) -airy(1,n2_p) 0;...
%                  airy(0,n1_n) airy(2,n1_n) 0 -(airy(2,n2_n)+1i*airy(0,n2_n));...
%                  airy(1,n1_n) airy(3,n1_n) 0 -(airy(3,n2_n)+1i*airy(1,n2_n))];
    fval = abs(det(AiMatrix));
end

function [waveFun] = DrawAi(Ebar_w,Fbar_w,Ebar_b,Fbar_b,V0bar_b,E0_b,mw,mb,Lw,Lb,z)
   %draw the wavefunciton
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    n1_p = -1*(pi/Fbar_w)^(2/3)*(Ebar_w - 1/2 * Fbar_w);
    n1_n = -1*(pi/Fbar_w)^(2/3)*(Ebar_w + 1/2 * Fbar_w);
    n2_p = -1*(pi/Fbar_b)^(2/3)*(Ebar_b - V0bar_b - 1/2 * Fbar_b);
    n2_n = -1*(pi/Fbar_b)^(2/3)*(Ebar_b - V0bar_b + 1/2 * Fbar_b);
%     AiMatrix = [ airy(0,n1_p) airy(2,n1_p) -airy(0,n2_p) 0;...
%                  airy(1,n1_p) airy(3,n1_p) -airy(1,n2_p) 0;...
%                  airy(0,n1_n) airy(2,n1_n) 0 -(airy(2,n2_n)+1i*airy(0,n2_n));...
%                  airy(1,n1_n) airy(3,n1_n) 0 -(airy(3,n2_n)+1i*airy(1,n2_n))];
AiMatrix = [ airy(0,n1_p), airy(2,n1_p), -airy(0,n2_p), 0;...
                 airy(1,n1_p)*(mb^(2/3)), airy(3,n1_p)*(mb^(2/3)), -airy(1,n2_p)*(mw^(2/3)), 0;...
                 airy(0,n1_n), airy(2,n1_n), 0, -(airy(2,n2_n)+1i*airy(0,n2_n));...
                 airy(1,n1_n)*(mb^(2/3)), airy(3,n1_n)*(mb^(2/3)), 0, -(airy(3,n2_n)+1i*airy(1,n2_n))*(mw^(2/3))];
    V = null(AiMatrix);
    a0 = V(1);
    b0 = V(2);
    a2 = V(3);
    a1 = V(4);
    E = Ebar_b*E0_b;
    F = Fbar_b*E0_b/(Lw*q);
    V0 = V0bar_b*E0_b;
    dz = abs(z(2)-z(1));
    z1 = z(find(z < -Lw/2));
    z2 = z(find(z >= -Lw/2 & z <= Lw/2));
    z3 = z(find(z > Lw/2));
    n2_n = -(2*mb/(q*hbar*F)^2)^(1/3)*(E-V0-q*F*z1);
    n1 = -(2*mw/(q*hbar*F)^2)^(1/3)*(E-q*F*z2);
    n2_p = -(2*mb/(q*hbar*F)^2)^(1/3)*(E-V0-q*F*z3);
    waveFun = zeros(1,length(z1)+length(z2)+length(z3));
    waveFun(1:length(z1)) = a1*(airy(2,n2_n) + 1i*airy(0,n2_n));
    waveFun(length(z1)+1:length(z2)+length(z1)) = a0*airy(0,n1)+b0*airy(2,n1);
    waveFun(1+length(z1)+length(z2):end) = a2*airy(0,n2_p);
    waveFun = abs(waveFun./sqrt(sum(abs(waveFun).^2) * dz));

end

% function [waveFun, Energy] = QW_FD1d(x, v, m, eigNum, flagPeriodic,tol)
% %% Solve 1D Schrodinger equation using Finite difference method
% %% 20110122 by TYB
% %% x, xgrid point
% %% v, potential
% %% m, effective mass
% %% eigNum, eigen number
% %% flagPeriodic: 1(default), periodic boundary, 0, zero boundary
% %% tol, 1e-16(default), tolerance for the function "eigs"
% %% waveFun, wave function
% %% Energy, bounded state energy.
% 
% if nargin <4
%     eigNum = 1;
% end
% if nargin <5
%     flagPeriodic = 0;
% end
% if nargin<6
%     tol = 1e-16;
% end
% options.tol = tol;
% 
% %%constant
% h = 6.62606896e-34; % J s       Planck constant
% hbar = h/2/pi;      % J s       reduced Planck constant
% q = 1.6021765e-19;  %C
% m0 = 9.10938215e-31; %kg
% 
% dx = abs(x(2)-x(1));
% v = v*q;
% m = m*m0;
% 
% %% get spdiags
% main_diag=(hbar/dx)^2./m + v; % main diagonal
% sub_diag=-(hbar/dx)^2./(2*m); % -1 and +1 diagonal
% N = length(main_diag); % matrix scale  N x N
% A_diag = zeros(N, 3);
% A_diag(:,3) = [0,sub_diag(1:N-1)];
% A_diag(:,2) = main_diag(:);
% A_diag(:,1) = [sub_diag(2:N),0];
% A = spdiags(A_diag,[-1,0,1],N,N);%asemble the matrix
% 
% if flagPeriodic
%     A(1,N) = -(hbar/dx)^2./(2*m(1));
%     A(N,1) = -(hbar/dx)^2./(2*m(end));
% end
% %%  caculate En and Fi
% NSTM = eigNum; % eigen root number 
% [Fic,En,flag] = eigs(A,NSTM,'sr',options); %solve the equation. why "sr"??
% %  [Fic,En,flag] = eigs(A,NSTM,'si',options);
% % [Fic,En,flag] = eigs(A,NSTM);
% waveFun = -Fic';
% waveFun = abs(waveFun);
% waveFun = waveFun/sqrt((sum(abs(waveFun).^2))*dx);
% Energy = diag(En)/q;
% end


