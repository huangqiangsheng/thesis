function [waveFun, Energy] = QW_FD1d(x, v, m, eigNum, flagPeriodic,tol)
%% Solve 1D Schrodinger equation using Finite difference method
%% 20110122 by TYB
%% x, xgrid point
%% v, potential
%% m, effective mass
%% eigNum, eigen number
%% flagPeriodic: 1(default), periodic boundary, 0, zero boundary
%% tol, 1e-16(default), tolerance for the function "eigs"
%% waveFun, wave function
%% Energy, bounded state energy.

if nargin <4
    eigNum = 1;
end
if nargin <5
    flagPeriodic = 0;
end
if nargin<6
    tol = 1e-16;
end
options.tol = tol;

%%constant
h = 6.62606896e-34; % J s       Planck constant
hbar = h/2/pi;      % J s       reduced Planck constant
q = 1.6021765e-19;  %C
m0 = 9.10938215e-31; %kg

dx = abs(x(2)-x(1));
v = v*q;
m = m*m0;

%% get spdiags
main_diag=(hbar/dx)^2./m + v; % main diagonal
sub_diag=-(hbar/dx)^2./(2*m); % -1 and +1 diagonal
N = length(main_diag); % matrix scale  N x N
A_diag = zeros(N, 3);
A_diag(:,3) = [0,sub_diag(1:N-1)];
A_diag(:,2) = main_diag(:);
A_diag(:,1) = [sub_diag(2:N),0];
A = spdiags(A_diag,[-1,0,1],N,N);%asemble the matrix

if flagPeriodic
    A(1,N) = -(hbar/dx)^2./(2*m(1));
    A(N,1) = -(hbar/dx)^2./(2*m(end));
end
%%  caculate En and Fi
NSTM = eigNum; % eigen root number 
[Fic,En,flag] = eigs(A,NSTM,'sr',options); %solve the equation. why "sr"??
%[Fic,En,flag] = eigs(A,NSTM);
waveFun = -Fic';
Energy = diag(En)/q;
end