function [Eb,wavlen] = VarMethod(ze,zh,phie,phih,mu,epsr)
%% "Modeling of self-electro-optic-effect devices" 
%   P. J. Mares and S. L. Chuang,1993
%% created by HQS 2011_12_4
%ze,zh对应phie,phih的格点
    h = 6.62606896e-34; % J s       Planck constant
    hbar = h/2/pi;      % J s       reduced Planck constant
    q = 1.6021765e-19;  %C
    m0 = 9.10938215e-31; %kg
    epsr0 = 8.85418782e-12;
    mu = mu*m0;
    C1 = hbar^2/(2 * mu*q);
    C2 = q/(4*pi*epsr*epsr0);
    startel = 50;
    endel = 150;
    sample_num = 5e1;
    lambda = linspace(startel,endel,sample_num);
%     fval = zeros(size(lambda));
%     for iter = 1:length(lambda)
%         fval(iter) = CalculEb(lambda(iter),C1,C2,ze,zh,phie,phih);
%     end
%     figure(100);
%     plot(lambda,fval); xlabel('lambda(A)');
%     hold on;
%     [gx,gy] = ginput(1);
%     plot([gx(1) gx(1)],[0 gy(1)],'r');
    [wavlen,Eb] = fminsearch(@(lambda)CalculEb(lambda,C1,C2,ze,zh,phie,phih),100,optimset('TolX',1e-32));
%     plot(wavlen,Eb,'r*');hold off
    wavlen = wavlen * 1e-10;
end

function Eb = CalculEb(lambda,C1,C2,ze,zh,phie,phih)
    Eb = 0;
    Eb1 = 0;
    f = 0;
    lambda = abs(lambda)*1e-10;
    dze = abs(ze(2)-ze(1));
    dzh = abs(zh(2)-zh(1));
    zh = repmat(zh,length(ze),1);
    ze = ze';
    ze = repmat(ze,1,length(zh));
    f = G(2*abs(zh - ze)/lambda);
    phih = abs(phih)';
    phih = repmat(phih,1,length(zh));
    Eb1 = sum(f.*phih.^2,1);
    Eb =  sum(Eb1.*abs(phie).^2) * dzh;

    Eb = C1/(lambda^2) - C2*2/lambda*Eb*dze;    
end

function f = G(x)
    f= ones(size(x));
    iter = find((x > 0 & x <= 6.8));
    f(iter) = 1 + (-8.9707e-1)*x(iter) + (-2.5262e-1)*x(iter).^2.*log(x(iter)/2)...
              +(2.2576e-1)*x(iter).^2+...
              (3.2373e-2)*x(iter).^3+(-4.1396e-4)*x(iter).^4.*log(x(iter)/2);
    iter = find(x > 6.8);
    f(iter) = 1./x(iter) - 3./x(iter).^3+45./x(iter).^5 -1575./x(iter).^7;
end

