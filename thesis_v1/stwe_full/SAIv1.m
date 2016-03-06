function output = SAIv1(x)
    Len = 100;ZS=50;ZL=50;
    %[Lens,Cps,Cpl,FlagRtn] = f1s_4(x,Len);
    %f1s_4 15 f2s_4 39 f3s_4 61
    [Lens,Cps,Cpl,FlagRtn] = f4s_1(x,Len);
    
    %fp 16 fp_2 32
    if (FlagRtn == 1)
        output = -1;
        return;
    end

    Lens = Lens*1e-6;    
    filter=0;% apply different filter mechanics
%     Cps=15e-15;Cpl=15e-15;% take pad capacitance into account.
    Cps=0;Cpl=0;
    BD3GSVm = get_f3dB_stwe(Lens,Cps,Cpl,ZS,ZL,filter,0);
    output = -BD3GSVm;

end

function [Lens,Cps,Cpl,FlagRtn] = f2s_1(x,Len) %bits = 15
% without pad capacitance and extensional segment
    if length(x) ~= 2
        error ('输入参数bit数不对');
    end
    
    lm1 = Len;
    lp1 = x(1);lp2 = x(2);    
    Cps = 0;
    Cpl = 0;
    
    Lens = [lp1;lm1;lp2];  
    FlagRtn =0;
end

function [Lens,Cps,Cpl,FlagRtn] = f3s_1(x,Len) %bits = 15
% without pad capacitance and extensional segment
    if length(x) ~= 3
        error ('输入参数bit数不对');
    end
    
    lp1 = x(1);    
    lm1 = Len;
    lo1 = x(2);
    lp2 = x(2);
    lm2 = Len;    
    lp3 = x(3);    
    Cps = 0;
    Cpl = 0;
    Lens = [lp1;lm1;lo1;lp2;lm2;lp3]; 
    FlagRtn = 0;
    
end

function [Lens,Cps,Cpl,FlagRtn] = f4s_1(x,Len) 
% without pad capacitance and extensional segment
    if length(x) ~= 4
        error ('输入参数bit数不对');
    end
    
    lp1 = x(1);    
    lm1 = x(2);
    lo1 = x(3);
    lp2 = x(3);
    lm2 = Len-lm1;    
    lp3 = x(4);    
    Cps = 0;
    Cpl = 0;
    Lens = [lp1;lm1;lo1;lp2;lm2;lp3]; 
    FlagRtn = 0;
    
end