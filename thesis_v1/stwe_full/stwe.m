function [GSVm,GSVm0, Ta, ABCD] = stwe(passive,active,optical,Lens,ZS,ZL,CS,CL)
%20100111 by TYB updated from old version.
%calculate the frequency response/see Y. Tang, et al, OC 281 (2008) 5177-5182.
%active, passive, optical:  structure, see stwe_setup
%Lens, array-type, segmented electrode length list, according to p,a,o,p,a,o,p,...a,p
%ZS,ZL, source resistance and load resistance
%CS,CL, pad capacitance for the inport and output pads.

nLen = length(Lens);
if nLen==1  % Conventional TW type
    oLen = [];
    mLen = Lens;          %普通的行波电极
    tLen = [0,Lens,0];
elseif mod(nLen,3)==0 % meandering TW type
    oLen = Lens;
    oLen(end)=[];
    oLen(1:3:end)=[];      %光波导的长度
    mLen = Lens(2:3:end-1);%有源混合波导
    tLen = Lens;
    tLen(3:3:end-1)=[];    %电极长度
else
    error ('wrong size of the Lens');
end

freq=active.freq;     %扫描频率
nElec  = length(tLen);% electrode section number
nFreq  = length(freq);% frequency point for calculation

% calculate the input resistance Zin  at every node.
Zin = zeros(nElec,nFreq); % input resistance looking at the load
Aa =1;Ba=0;Ca=0;Da=1;
Ca = j*2*pi*freq*CL;  %% capacitance
for ii = nElec:-1:1% calculate Zin from right to the left.
    if mod(ii,2) == 1
        [Aa,Ba,Ca,Da] = preABCD(Aa,Ba,Ca,Da,passive.Zc,passive.Gamma,tLen(ii));
        Zin(ii,:) = (Aa.*ZL+Ba)./(Ca.*ZL+Da);
    else
        [Aa,Ba,Ca,Da] = preABCD(Aa,Ba,Ca,Da,active.Zc,active.Gamma,tLen(ii)); 
        Zin(ii,:) = (Aa.*ZL+Ba)./(Ca.*ZL+Da);
    end
end
[Aa,Ba,Ca,Da] = multipyMatrix(1,0,j*2*pi*freq*CS,1,Aa,Ba,Ca,Da);% get the ABCD for the whole device.
Zin(1,:) = (Aa.*ZL+Ba)./(Ca.*ZL+Da);% Zin at source port
Ta = 20*log10(abs((Zin(1,:)-ZS)./(Zin(1,:)+ZS)));% reflection to the source. 

%calculate the voltage at every discontinuity.
Vtemp = Zin(1,:)./(Zin(1,:)+ZS);% V1 at input port 
Vdis = Vtemp; % V1, at input port
for ii = 1:nElec-1% V, from left to right
    if mod(ii,2) == 1
        Vtemp = VotageSection(tLen(ii),Vtemp,passive.Zc,Zin(ii+1,:),passive.Gamma,tLen(ii));
    else
        Vtemp = VotageSection(tLen(ii),Vtemp,active.Zc,Zin(ii+1,:),active.Gamma,tLen(ii));
    end
    Vdis = [Vdis;Vtemp];
end
Vtemp = VotageSection(tLen(nElec),Vtemp,passive.Zc,ZL,passive.Gamma,tLen(nElec));
% Votage at the output port/over ZL
Vdis = [Vdis;Vtemp];

c = 2.997956376932163e+008 ;% light speed
betao = optical.ng*2*pi*freq/c;% the wavenumber of the optical wave packet

SVm  = 0;% sum of the small signal RF link gain. 
oPos = cumsum(oLen);% node position along optical wg, starting with 0
oPos = [0;oPos];
for ii = 2:2:nElec
    SVmtemp = IntVotageSection(Vdis(ii,:).*exp( j*betao*oPos(ii-1) ),active.Zc,Zin(ii+1,:),active.Gamma,betao,tLen(ii));
    SVm = SVm + SVmtemp./(1+j*2*pi*freq*active.tao);
end
Response = SVm*2/sum(mLen);
GSVm = 20*log10( abs(Response) ); 
%normalize it to ideal lossless electrode with perfectly matched impedance.

%%calculate DC response. Ignore all capacitances and inductances. 
SVm0 = 0;
ZDC  = sum(tLen)*passive.Rdc+ZS+ZL;%% ignore the active part's influence.
for ii = 2:2:nElec
	SVm0 = SVm0+(sum(tLen(ii:nElec))*passive.Rdc+ZL-0.5*tLen(ii)*passive.Rdc)/ZDC*tLen(ii);
end    
GSVm0 = 20*log10( abs(SVm0*2/sum(mLen)) ); 
%normalize it to ideal lossless electrode with perfectly matched impedance.
ABCD = [Aa;Ba;Ca;Da];

%% additional function may be needed.
% [s11,s12,s21,s22] = ABCD2S(Aa,Ba,Ca,Da,ZS,ZL);
% S11 = 20*log10(abs(s11));
% S21 = 20*log10(abs(s21));%convert Np/m to dB/m
% %calculat the chacteristic impedance.
% delta = (Aa-Da)./Ca/2;
% ZZ = sqrt((Aa+Da).*(Aa+Da)-4)./2./Ca;
% ZZ = ((real(ZZ)>0)*2-1).*ZZ;
% Zc1 = delta+ZZ;
% Zc2 = -(delta-ZZ);
% % calculate the microwave loss, including the first and the last p sections
% gam = acosh(Aa/2+Da/2);
% gam_i = find(imag(gam)<0);
% gam(gam_i:end) = gam(gam_i:end)+j*2*pi;
% gam = gam./sum(Len);
% alpha = 8.6859*real(gam);
% 
% % calculate microwave group index normalized to optical length, excluding the first and the last p sections 
% Aa =1;Ba=0;Ca=0;Da=1;
% for ii = Ns-1:-1:2
%     if mod(ii,2) == 1
%         [Aa,Ba,Ca,Da] = preABCD(Aa,Ba,Ca,Da,Zcp,gammap,Len(ii));
%     else
%         [Aa,Ba,Ca,Da] = preABCD(Aa,Ba,Ca,Da,Zcm,gammam,Len(ii)); 
%     end
% end  
% gam2 = acosh(Aa/2+Da/2);
% gam2_i = find(imag(gam2)<0);
% gam2(gam2_i:end) = gam2(gam2_i:end)+j*2*pi;
% gam2 = gam2./(OLen(end)+Len(end-1)); % index nomarlized to optical waveguide length
% nm = imag(gam2)*c./freq;% effective microwave index excluding the two leading electrode sections.
end

function Vz = VotageSection(z,Vpre,ZC,ZL,gamma,length)
%calculate the voltage at z point.

TaoL = (ZL-ZC)./(ZL+ZC);
Vz  = Vpre.*( exp(-gamma*(z-length))+TaoL.*exp(gamma*(z-length)) )./(exp(gamma*length)+TaoL.*exp(-gamma*length));

end

function IntVz = IntVotageSection(Vpre,ZC,ZL,gamma,betao,length)
%calculate the voltage integal along an electrode segment。
% Vpre: V at beginning， length is the length of the Integral range
% ZC:characteristic impedance, gamma: the propagation constant
% ZL is the input impedance at the end point
% betao is the wavenumber assiciated with optical wave packet.

TaoL = (ZL-ZC)./(ZL+ZC);
IntVz  = Vpre./(1+TaoL.*exp(-gamma*2*length)).*( (exp((j*betao-gamma)*length)-1)./(j*betao-gamma)...
    +TaoL.*exp(-2*gamma*length).*(exp((j*betao+gamma)*length)-1)./(j*betao+gamma) ) ;

end

function [Ao,Bo,Co,Do] = multipyMatrix(A1,B1,C1,D1,A2,B2,C2,D2)
% A1*A2

Ao = A1.*A2+B1.*C2;
Bo = A1.*B2+B1.*D2;
Co = C1.*A2+D1.*C2;
Do = C1.*B2+D1.*D2;
end

function [Ao,Bo,Co,Do] = preABCD(Ai,Bi,Ci,Di,Zc,gamma,Len)
% 从右往左计算级联传输线的总的等效ABCD矩阵。 

A = cosh(gamma*Len);
B = Zc.*sinh(gamma*Len);
C = sinh(gamma*Len)./Zc;
D = A;

[Ao,Bo,Co,Do] = multipyMatrix(A,B,C,D,Ai,Bi,Ci,Di);
% Ao = A.*Ai+B.*Ci;
% Bo = A.*Bi+B.*Di;
% Co = C.*Ai+D.*Ci;
% Do = C.*Bi+D.*Di;
end

function [S11,S12,S21,S22]=ABCD2S(A,B,C,D,ZS,ZL)
% convert the ABCD matrix to S-Parameters matrix;
% reference: David M. Polar micrwave engineering 2rd P211

den = A.*ZL+B+C.*ZS.*ZL+D.*ZS;
S11 = (A.*ZL+B-C.*ZS.*ZL-D.*ZS)./den;
S12 = 2*sqrt(ZL.*ZS).*(A.*D-B.*C)./den;
S21 = 2*sqrt(ZL.*ZS)./den;
S22 = (-A.*ZL+B-C.*ZS.*ZL+D.*ZS)./den;

end

