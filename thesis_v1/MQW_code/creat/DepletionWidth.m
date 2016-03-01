function DeltaV = DepletionWidth(xp,x,epsr,N,V0)
%xp 为p的耗尽层的厚度
%x(n)~x(n-1)为epsr(n-1),N(n-1)的材料
%V0为外加电压
%Phi(x) = Ax^2 + Bx + C

q = 1.6021765e-19;  %C
epsr0 = 8.85418782e-12;
xp = xp * 1e-6;
x = [-xp x];
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
V = C(end) - B(end)^2/(4*A(end));
DeltaV = abs(V-V0);
end