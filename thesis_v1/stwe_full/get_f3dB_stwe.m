function BD3GSVm = get_f3dB_stwe(Lens,CS,CL,ZS,ZL,filter,plotFlag,reload_data);
% stwe's factory to get 3dB bandwidth
if nargin <= 7
    reload_data=0;
end
%run stwe_setup;  % setup file including ZC, ZS, data.
global passive active optical stwe_setup_flag;
if reload_data
    stwe_setup_flag = []
end

if isempty(stwe_setup_flag)|| stwe_setup_flag==0
    stwe_setup();
end
% calculate
[GSVm,GSVm0,Ta,ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL);

% post process
dfsp = active.freq(2)-active.freq(1);
BD3GSVm = getBW3dB(GSVm,dfsp,2);
BD3i = round(BD3GSVm/dfsp);

if filter == 1 
%limitation
	exceed = GSVm0-min(GSVm(1:round(0.4*BD3i)))-1.0;
	if (exceed>0) 
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0~0.4×BD3i 内，有值小于-1dB的情形
    exceed = GSVm0-min(GSVm(round(0.4*BD3i):round(0.6*BD3i)))-1.5;
	if (exceed>0) 
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0.4~0.6×BD3i 内，有值小于-1.5dB的情形
    exceed = GSVm0-min(GSVm(round(0.6*BD3i):round(0.8*BD3i)))-2.0;
	if (exceed>0) % apply penalty
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0.6~0.8×BD3i 内，有值小于-2dB的情形
elseif filter == 2
	exceed = GSVm0-min(GSVm(1:round(0.6*BD3i)))-1.5; 
	if (exceed>0) % apply penalty
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0~0.6×BD3i 内，有值小于-1.5dB的情形
elseif filter == 3  % corresponding to the response of a TWL EAM.  
%限制条件，对频率响应曲线的两个限制 
	exceed = GSVm0-min(GSVm(1:round(0.5*BD3i)))-1.0;
	if (exceed>0) % apply penalty
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0~0.5×BD3i 内，有值小于-1dB的情形
    exceed = GSVm0-min(GSVm(round(0.5*BD3i):round(0.75*BD3i)))-2.0;
	if (exceed>0) % apply penalty
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0.5~0.75×BD3i 内，有值小于-2dB的情形   
elseif filter == 4 % 40G >-1dB
    exceed = GSVm0-min(GSVm(1:30))-1.5;
	if (exceed>0) % apply penalty
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0~0.5×BD3i 内，有值小于-1dB的情形
elseif filter == 5
	exceed = GSVm0-min(GSVm(1:round(0.7*BD3i)))-2; 
	if (exceed>0) % apply penalty
		BD3GSVm=BD3GSVm-10*(exp((exceed)*3)-1)*dfsp;
	end  %剔除 0~0.6×BD3i 内，有值小于-1.5dB的情形    
end

if plotFlag==1
outputstring=sprintf('GSVm(1)=%f, GSVm0=%f, detal=%f, f3=%e',GSVm(1),GSVm0,GSVm(1)-GSVm0,BD3GSVm);
disp(outputstring);
figure(10);hold on
[AX,H1,H2]=plotyy(active.freq, Ta,active.freq, GSVm);
set(AX(1),'YLim',[-30,0]);
set(AX(2),'YLim',[-5,0.3]);
set(AX(1),'XTickLabel',[]);
set(AX(1),'YTickMode','auto');
set(AX(2),'YTickMode', 'auto');
set(H1,'LineWidth',2.5);
set(H2,'LineWidth',2.5);
grid(AX(1));
grid(AX(2));
end
end

function  BD3 = getBW3dB(data,df,type)
% data，frequency response， df，deta frequency， 
% type for the ref point: 1- maximal point，2- the first point）
if type ==1 
    db3p = max(data(:))-3;%get the 3dB point in the curve
else if type ==2
        db3p = data(1)-3;
    else
        error('Wrong input： type');
    end
end
        
B3p = find(data<db3p,1,'first'); %seach the 3dB frequency
if isempty(B3p)||B3p==length(data)
    BD3 = length(data)*df;
%    disp('Increase the scan range!');
    return;
else
    BD3 = B3p-(data(B3p)-db3p)/(data(B3p)-data(B3p-1));% interpolate  
end

B3p = find(data>db3p,1,'first'); %seach the 3dB frequency
if B3p >1
    BDs = B3p+(db3p-data(B3p))/(data(B3p)-data(B3p-1));% interpolate 
    BD3 = BD3-BDs;
    %B3p == 1  %if B3p == 1，ignore，default is 0。
end
BD3 = BD3*df;% convert to Hz。
end