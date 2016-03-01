function [Phix,Phi]=Draw_DepletionWidth(x,epsr,N,Eband,interp_num)
%xp 为p的耗尽层的厚度
%x(n)~x(n-1)为epsr(n-1),N(n-1)的材料
%V0为外加电压
%Phi(x) = Ax^2 + Bx + C
    
q = 1.6021765e-19;  %C
epsr0 = 8.85418782e-12;
N = N*1e6;
A = q*N./(2*epsr*epsr0);
B = zeros(size(A));
C = zeros(size(A));
for iter = 2:length(A)
B(iter) = epsr(iter-1)/epsr(iter)*...
          (2*A(iter-1)*(x(iter)-x(iter-1)) + B(iter-1));
C(iter) = A(iter-1)*(x(iter)-x(iter-1))^2+...
          B(iter-1)*(x(iter)-x(iter-1))+C(iter-1);
end
Phix = zeros(1,length(A)*interp_num);
PhiE = zeros(1,length(A)*interp_num);
Phi = zeros(1,length(A)*interp_num);
for iter = 1:(length(x)-1)
   dx = linspace(0,x(iter+1)-x(iter),interp_num);
   Phix((iter-1)*interp_num+1 : iter*interp_num) = dx + x(iter);
   PhiE((iter-1)*interp_num+1 : iter*interp_num) = Eband(iter);
   Phi((iter-1)*interp_num+1 : iter*interp_num) = A(iter)*dx.^2 +...
                                                  B(iter)*dx + C(iter); 
end
Phi = -Phi;%相对于电子的势能
Phi = Phi + PhiE;
% plot(Phix,Phi);
end