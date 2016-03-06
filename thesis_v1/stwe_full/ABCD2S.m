function [S11,S12,S21,S22]=ABCD2S(A,B,C,D,ZS,ZL)
% convert the ABCD matrix to S-Parameters matrix;
% reference: David M. Polar micrwave engineering 2rd P211

den = A.*ZL+B+C.*ZS.*ZL+D.*ZS;
S11 = (A.*ZL+B-C.*ZS.*ZL-D.*ZS)./den;
S12 = 2*sqrt(ZL.*ZS).*(A.*D-B.*C)./den;
S21 = 2*sqrt(ZL.*ZS)./den;
S22 = (-A.*ZL+B-C.*ZS.*ZL+D.*ZS)./den;

end

