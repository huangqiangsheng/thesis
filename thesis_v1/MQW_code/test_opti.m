function val = test_opti(x)
%% 优化量子阱的参数
% well_w,well_v,ww分别为掺杂比例和well的宽度
% barrier_w,barrier_v,bw分别掺杂比例和barrier的宽度
% 2012-2-27,HQS 修正了strain balance的评价标准.
well_w = x(1)/100;
well_v = x(2)/100;
barrier_w = x(3)/100;
barrier_v = x(4)/100;
ww = 110;
goal_lambda_CHH = 1550;
goal_lambda_CLH= 1485;
CTE = 20;
CTM = 0;
CD = 0;
Coef_width = 0;
Coef_ezz = 0;%500
Coef_DEc = 0;
Coef_DEvH = 0;
Coef_DEvL = 0;
Coef_WEg = 0;

ww = ww*1e-10;
bw = ww;%假设量子阱的宽度相等
well = InAlGaAs_params(well_w,well_v,1550);
barrier = InAlGaAs_params(barrier_w,barrier_v,1550);
% well = InGaAsP_params(well_w,well_v,1550);
% barrier = InGaAsP_params(barrier_w,barrier_v,1550);
MQW.barrier = barrier;
MQW.well = well;
MQW.tb = bw;
MQW.tw = ww;
MQW.num_b = 2;
MQW.num_w = 1;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

flag1 = 0;
if well.Ec>=barrier.Ec
    flag1 = 1;
elseif well.EvH <= barrier.EvH
    flag1 = 1;
elseif well.EvL <= barrier.EvL
    flag1 = 1;
end

if flag1 ~= 1
    tol=1e-32; 
    sample_num = 200;
    x = linspace(-MQW.t/2, MQW.t/2, sample_num); %sample point
    dx = x(2) - x(1);
    MQW.grid = x; %grid point
    Ebias = 0; % applied electric field
    Vb = (1:sample_num)*Ebias*dx - fix(sample_num/2)*Ebias*dx; %applied bias potential V.
    %%
    v0 = MQW.barrier.Ec -MQW.well.Ec;
    v = ones(1,sample_num)*v0;
    m = ones(1,sample_num)*MQW.barrier.me;
    for iter = 1:MQW.num_w;
        v(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = 0;
        m(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = MQW.well.me;
    end
    v = v + Vb;
    [Ec,exitflag] = TMM_opt(x, v, m,tol);
    if exitflag ~= 1
        flag1 = 1;
    end
    %%
    v0 = MQW.barrier.EvH -MQW.well.EvH;
    v = ones(1,sample_num)*v0;
    m = ones(1,sample_num)*MQW.barrier.mHHz;
    for iter = 1:MQW.num_w;
        v(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = 0;
        m(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) =  MQW.well.mHHz;
    end
    v = v + Vb;
    [EvH,exitflag] = TMM_opt(x, -v, m,tol);
    if exitflag ~= 1
        flag1 = 1;
    end
    %%
    v0 = MQW.barrier.EvL -MQW.well.EvL;
    v = ones(1,sample_num)*v0;
    m = ones(1,sample_num)*MQW.barrier.mLHz;
    for iter = 1:MQW.num_w;
        v(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) = 0;
        m(find(x>(-MQW.t/2+MQW.tb+ (iter-1)* MQW.period) & x < (-MQW.t/2+(iter)*MQW.period))) =  MQW.well.mLHz;
    end
    v = v + Vb;
    [EvL,exitflag] = TMM_opt(x, -v, m,tol);
    if exitflag ~= 1
        flag1 = 1;
    end
    %%
    MQW.Ec = Ec + MQW.well.Ec;     % conduct band
    MQW.EvH = -EvH + MQW.well.EvH;  % for heavy hole
    MQW.EvL = -EvL + MQW.well.EvL;  % for light hole
    evnm =1239.842;% 1 eV = evnm/wl; evnm = hc/q*1e9
    MQW.PL_cHH = evnm./(MQW.Ec-MQW.EvH);
    MQW.PL_cLH = evnm./(MQW.Ec-MQW.EvL);
    if flag1 ~= 1
        val = CTE*((MQW.PL_cHH - goal_lambda_CHH))^2+...
               CTM*((MQW.PL_cLH - goal_lambda_CLH))^2+...
               CD*((MQW.PL_cHH - MQW.PL_cLH))^2+...
               Coef_ezz * abs(MQW.well.ezz-0.009)*1e3+...
               Coef_ezz * abs(MQW.barrier.ezz+0.009)*1e3+...
               (Coef_width *((ww-110e-10)*1e9))^2 + ...
               (Coef_WEg*(evnm/MQW.well.Eg-1550)^2)+...
               (Coef_DEvH *(-(MQW.barrier.EvH -MQW.well.EvH) -0.12))^2+...
               (Coef_DEvL*(-(MQW.barrier.EvL -MQW.well.EvL) - 0.12))^2+...
               (Coef_DEc*((MQW.barrier.Ec -MQW.well.Ec) - 0.14))^2;
    else
        val = 1e10;
    end
else
    val = 1e10;
end

end
function [Energy,exitflag] = TMM_opt(x, v, m, tol)

    [E0,fval] = fminsearch(@(E)F_E(E,x, v, m),0, optimset('TolX',tol));
    exitflag = 1;
    if fval > 0.1
        exitflag = 0;
    end
    Energy = E0;
end

function f22 = F_E(E,x, v, m)
%% 20111124 by HQS
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    v = v * q;
    m = m * m0;
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

