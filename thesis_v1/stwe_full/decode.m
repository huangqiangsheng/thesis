function output = decode(x,para)
%将二进制编码的基因翻译为物理参数 
%2007-07-09, 修改为通用型的decode函数

%x 码字
%para， n*3 matrix    para(:,1) bit数， para(:,2) 最小值， para(:,3) 最大值
%n为decode输出参数个数
%%
    
    %采用二进制编码,处理输入参数
    bits = length(x);%输入参数bit数
    if(sum(para(:,1))>bits) 
		sum(para(:,1))
		bits
        error ('输入参数bit数太少');
    end
    paraNum = size(para,1); %输出参数个数
    output = zeros(paraNum,1);

    argu = zeros(1,paraNum);%物理参数的十进制值(0~1)
%% 中间部分的代码使恢复时，产生随机效果也可以直接用xi = 1:bits，然后argu(ii) =  x(sum(para(1:(ii-1),1))+1:sum(para(1:ii,1)))*weight2'/(2^para(ii,1)-1);
	sqr = floor(sqrt(bits));
	tim = floor(bits/sqr);
	mi  = zeros(tim,sqr);
	mi(1:sqr,:)=magic(sqr);% try to separate the bits converted from the same parameter
    if tim>sqr
        mi(sqr+1:end,:)=reshape((sqr*sqr+1):sqr*tim,(tim-sqr),sqr);
    end
	xi = 1:bits;
	xi(1:tim*sqr)=mi(:)';
%%	
    for ii = 1:paraNum
        weight2 = 2.^(para(ii,1)-1:-1:0);  
        %argu(ii) = x(sum(para(1:(ii-1),1))+1:sum(para(1:ii,1)))*weight2'/(2^para(ii,1)-1); 
		argu(ii) = x(xi(sum(para(1:(ii-1),1))+1:sum(para(1:ii,1))))*weight2'/(2^para(ii,1)-1);
        btm = para(ii,2);
        top = para(ii,3);
        output(ii) = btm + (top-btm)*argu(ii);  
    end
end