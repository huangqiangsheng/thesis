function [mseq]= PRBSgenerator(GeneratorPolynomial, allBits)
%PRBSgenerator function: 根据输入多项式系数，输出PRBS码，psedurandom bit sequences
% GeneratorPolynomial: PN 的多项式系数，由高到低， 如 [1 0 0 0 0 0 1 0 1] 或者 [8 2 0],或者直接输入8.
% 如果 GeneratorPolynomial为一个整数，则表示为选择随机码长度，自动选择相应的多项式（支持范围为2～53）。
% 常用的有 [21 19 0]； 参考 matlab 帮助 PN Sequence Generator
%
% fbconnection =[5 3 0];
% mseq = PRBSgenerator(fbconnection)
% rGenerator PolynomialrGenerator Polynomial

if nargin<2 
	N = 2^n-1;
else
	N = allBits;
end

len = length(GeneratorPolynomial);
if len == 1
	if GeneratorPolynomial>53 || GeneratorPolynomial<2
	 error('GeneratorPolynomial is set to be in [2,53]');
	end
GP={  %PolynomialrGenerator Polynomial cell
[2 1 0]				
[3 2 0]				
[4 3 0]				
[5 3 0]				
[6 5 0]				
[7 6 0]				
[8 6 5 4 0]			
[9 5 0]				
[10 7 0]			
[11 9 0]			
[12 11 8 6 0]		
[13 12 10 9 0]		
[14 13 8 4 0]		
[15 14 0]			
[16 15 13 4 0]		
[17 14 0]			
[18 11 0]			
[19 18 17 14 0]	
[20 17 0]
[21 19 0]
[22 21 0]
[23 18 0]
[24 23 22 17 0]
[25 22 0]
[26 25 24 20 0]
[27 26 25 22 0]
[28 25 0]
[29 27 0]
[30 29 28 7 0]
[31 28 0]
[32 31 30 10 0]
[33 20 0]
[34 15 14 1 0]
[35 2 0]
[36 11 0]
[37 12 10 2 0] 
[38 6 5 1 0] 
[39 8 0] 
[40 5 4 3 0]		
[41 3 0]			
[42 23 22 1 0]		
[43 6 4 3 0]		
[44 6 5 2 0]		
[45 4 3 1 0]		
[46 21 10 1 0]		
[47 14 0] 
[48 28 27 1 0] 
[49 9 0] 
[50 4 3 2 0] 
[51 6 3 1 0] 
[52 3 0] 
[53 6 2 1 0]
};
GeneratorPolynomial = GP{GeneratorPolynomial-1};
end
[n,ni] = max(GeneratorPolynomial(:));
if n>1 
	coef = zeros(1,n);
	GeneratorPolynomial(ni) = [];
	coef(n-GeneratorPolynomial) = 1;
else 
	n = len-1;
	coef = GeneratorPolynomial(2:end);
end

mseq = zeros(1,N);
register = ones(1,n);%定义移位寄存器的初始状态
for i = 1:N
    newregister= mod(sum(coef.*register),2);
	register(2:end) = register(1:end-1);
    register(1) = newregister;
    mseq(i) = register(n);
end

end

	
		


