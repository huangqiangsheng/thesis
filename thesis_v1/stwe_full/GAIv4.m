function output = GAIv4(x)
% 20100111 by TYB, modified for new version get_f3dB_stwe.m
% interface script for GA toolbox. 
% 2007-07-09  by Yongbo Tang. 
% 2009-04-29  by Y.B Tang add function f1s_2 f2s_2 for extensional segment
% 2009-04-30  by Y.B. T   modify the function f**, modify decode (adjust the bit position)
% 2010-01-20  cowork with get_f3dB_stwe.m and related, add comments.
    Len = 200;ZS=50;ZL=50;
    %[Lens,Cps,Cpl,FlagRtn] = f1s_4(x,Len);
    %f1s_4 15 f2s_4 39 f3s_4 61
    [Lens,Cps,Cpl,FlagRtn] = f1s_4(x,Len);
    %fp 16 fp_2 32
    if (FlagRtn == 1)
        output = -1;
        return;
    end

    Lens = Lens*1e-6;    
    filter=0;% apply different filter mechanics
    Cps=15e-15;Cpl=15e-15;% take pad capacitance into account.
    BD3GSVm = get_f3dB_stwe(Lens,Cps,Cpl,ZS,ZL,filter,0);
    
    output = -BD3GSVm;
end

% perpare the input parameters for get_f3dB_stwe. 
function [Lens,Cps,Cpl,FlagRtn] = f3s_1(x,Len) %bits = 73
% 3 active segments, with pad capacitance
    global para;
    para= [%bits, start, end; 
            7,0,635;   %lp1     1
            8,0,255;   %lm1     2
            7,0,635;   %lo1     3
            8,0,1275;   %lp2     4
            8,0,255;    %lm2    5
            7,0,635;    %lo2    6
            8,0,1275;    %lp3    7
            %lm3 = Len - lm1
            8,0,1275;    %lp4    8
            7,0,127;    %Cps    9
            7,0,127     %Cpl    10
                    ];

    out = decode(x,para);
    lp1 = out(1);lm1 = out(2);
    lo1 = out(3);lp2 = out(4);lm2 = out(5);
    lo2 = out(6);lp3 = out(7);lm3 = round(Len-lm1-lm2);
    lp4 = out(8);
    Cps = out(9)*1e-15;
    Cpl = out(10)*1e-15;
    
    Lens = [lp1;lm1;lo1;lp2;lm2;lo2;lp3;lm3;lp4];
    FlagRtn = 0;
    if (lm3<0-1e-8 || lm2<1-1e-8 || lm1<45 || lm2<45 || lm3 <45 ) 
        FlagRtn = 1;
    end 

    if (lo1>(lp2-20+1e-8) ||lo2>(lp3-20+1e-8) || lp1<15-1e-8 || lp3<15-1e-8) 
        FlagRtn = 1;
    end  
    
    if (Cps<-9.9e-15 || Cpl<-9.9e-15)
        FlagRtn = 1;
    end
    
end

function [Lens,Cps,Cpl,FlagRtn] = f3s_4(x,Len) %bits = 59
% 3 active segmetns, without pad capacitance and extensional segment
    global para;
    %Len = 165;
    para= [%bits, start, end;  
            7,0,635;   %lp1     1
            8,0,255;   %lm1     2
            7,0,635;   %lo1     3
            8,0,1275;   %lp2     4
            8,0,255;    %lm2    5
            7,0,635;    %lo2    6
            8,0,1275;    %lp3    7
            %lm3 = Len - lm1
            8,0,1275;    %lp4    8
%            7,0,127;    %Cps    9
%            7,0,127     %Cpl    10
                    ];

    out = decode(x,para);
    lp1 = out(1);lm1 = out(2);
    lo1 = out(3);lp2 = out(4);lm2 = out(5);
    lo2 = out(6);lp3 = out(7);lm3 = round(Len-lm1-lm2);
    lp4 = out(8);
    Cps = 0*1e-15;
    Cpl = 0*1e-15;
    
    Lens = [lp1;lm1;lo1;lp2;lm2;lo2;lp3;lm3;lp4];
    FlagRtn = 0;
    if (lm3<0-1e-8 || lm2<1-1e-8 || lm1<30 || lm2<30 || lm3 <30 ) 
        FlagRtn = 1;
    end 

    if (lo1>(lp2-40+1e-8) ||lo2>(lp3-40+1e-8) || lp1<15-1e-8 || lp3<15-1e-8) 
        FlagRtn = 1;
    end      
end

function [Lens,Cps,Cpl,FlagRtn] = f2s_1(x,Len)% bits= 51  
    global para;
    %Len = 165;
    para= [%bits, start, end;  
            7,0,635;   %lp1
            8,0,255;   %lm1
            7,0,635;   %lo1
            7,0,635;   %lp2
            %lm2 = Len - lm1
            8,0,1275;    %lp3
           7,0,127;    %Cps
           7,0,127     %Cpl
                    ];
    out = decode(x,para);
    lp1 = out(1);lm1 = out(2);
    lo1 = out(3);lm2 = Len-lm1;
    lp2 = out(4);lp3 = out(5);
    Cps = out(6)*1e-15;
    Cpl = out(7)*1e-15;
    Lens = [lp1;lm1;lo1;lp2;lm2;lp3];
    FlagRtn = 0;
    if (lm2<0-1e-8 || lm1<70 || lm2<70) 
        FlagRtn = 1;
    end 

    if (lo1>(lp2-20+1e-8) || lp1<15-1e-8 || lp3<15-1e-8) 
        FlagRtn = 1;
    end   
    if (Cps<-9.9e-15 || Cpl<-9.9e-15)
        FlagRtn = 1;
    end    
end

function [Lens,Cps,Cpl,FlagRtn] = f2s_4(x,Len)% bits= 39
% without pad capacitance and extensional segment
    global para;
    %Len = 165;
    para= [%bits, start, end; 
         	7,0,635;   %lp1
            8,0,255;   %lm1
            8,0,1275;   %lo1
            8,0,1275;   %lp2
 %           8,0,255;%lm2 = Len - lm1
            8,0,1275    %lp3
                    ];
    out = decode(x,para);
    lp1 = out(1);lm1 = out(2);
    lo1 = out(3);lp2 = out(4);
    lm2 = Len-lm1;lp3 = out(5);
    Cps = 0*1e-15;
    Cpl = 0*1e-15;
    Lens = [lp1;lm1;lo1;lp2;lm2;lp3];
    FlagRtn = 0;
    if (lm1<20 || lm2<20) 
        FlagRtn = 1;
    end 

    if (lo1>(lp2-40+1e-8) || lp1<15-1e-8 || lp3<15-1e-8) 
        FlagRtn = 1;
    end   
end

function [Lens,Cps,Cpl,FlagRtn] = f1s_1(x,Len) %bits = 29
    global para;

    %Len = 165;
    para= [ %bits, start, end; 
	        7,0,635;   %lp1
%             8,0,255;   %lm1
            8,0,1275;    %lp3
            7,0,127;     %Cps
            7,0,127     %Cpl
                    ];
    out = decode(x,para);
    lp1 = out(1);lp2 = out(2);
    lm1 = Len;
    Cps = out(3)*1e-15;
    Cpl = out(4)*1e-15;
    
    Lens = [lp1;lm1;lp2];  
    FlagRtn =0;
    if (lp1<15-1e-8 || lp2<15-1e-8 ) 
        FlagRtn =1;
    end 
 
end

function [Lens,Cps,Cpl,FlagRtn] = f1s_4(x,Len) %bits = 15
% without pad capacitance and extensional segment
    global para;

    %Len = 165;
    para= [ 
    %       6,0,63;   %ls1
	%                  %ls2 = Len - ls1-lm1
			7,0,635;   %lp1
    %        7,Len-63.5,Len;   %lm1
            8,0,1275    %lp2
                    ];
    out = decode(x,para);
    lm1 = Len;
    lp1 = out(1);lp2 = out(2);    
    Cps = 0;
    Cpl = 0;
    Lens = [lp1;lm1;lp2];  
    FlagRtn =0;
    if (lp1<15-1e-8 || lp2<15-1e-8 ) 
        FlagRtn =1;
    end
	if (lm1 < 20 )
		FlagRtn =1;
    end 
end

function [Lens,Cps,Cpl,FlagRtn] = fp(x,Len,ns)% bits= 16
% periodic electrode, ns means the segment number
    global para;
    para= [ 8,1,256;% dtp/dm
            8,1,256;% dOp/dm
                    ];
    out = decode(x,para);
    ratio1 = out(1)/para(1,3)*16;
    ratio2 = out(2)/para(1,3)*16;
    lm= Len/ns;
    lo= lm*ratio1;
    lp= lm*ratio2;
    Lens = zeros(3*ns,1);
    Lens(2:3:end-1)=lm;
    Lens(3:3:end-1)=lo;
    Lens(4:3:end-1)=lp;
    Lens(1)=lp/2;
    Lens(end)=lp/2;
    Cps = 0*1e-15;
    Cpl = 0*1e-15;
    FlagRtn = 0;

    if (lo>(lp-20+1e-8)) 
        FlagRtn = 1;
    end   
end

function [Lens,Cps,Cpl,FlagRtn] = fp_2(x,Len,ns)% bits= 32
% take the first and the last TML lengths as variables. 
% periodic electrode
    global para;
    para= [ 8,1,256;% dtp/dm
            8,1,256;% dOp/dm
            8,1,256;% 
            8,1,256;%
                    ];
    out = decode(x,para);
    ratio1 = out(1)/para(1,3)*8;
    ratio2 = out(2)/para(1,3)*8;
    rt3 = out(3)/para(1,3)*16;
    rt4 = out(4)/para(1,3)*16;
    lm= Len/ns;
    lo= lm*ratio1;
    lp= lm*ratio2;
    Lens = zeros(3*ns,1);
    Lens(2:3:end-1)=lm;
    Lens(3:3:end-1)=lo;
    Lens(4:3:end-1)=lp;
    Lens(1)=lm*rt3;
    Lens(end)=lm*rt4;
    Cps = 0*1e-15;
    Cpl = 0*1e-15;
    FlagRtn = 0;

    if (lo>(lp-20+1e-8)) 
        FlagRtn = 1;
    end   
end



    

