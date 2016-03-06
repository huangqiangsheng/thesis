function output = decode(x,para)
%�������Ʊ���Ļ�����Ϊ������� 
%2007-07-09, �޸�Ϊͨ���͵�decode����

%x ����
%para�� n*3 matrix    para(:,1) bit���� para(:,2) ��Сֵ�� para(:,3) ���ֵ
%nΪdecode�����������
%%
    
    %���ö����Ʊ���,�����������
    bits = length(x);%�������bit��
    if(sum(para(:,1))>bits) 
		sum(para(:,1))
		bits
        error ('�������bit��̫��');
    end
    paraNum = size(para,1); %�����������
    output = zeros(paraNum,1);

    argu = zeros(1,paraNum);%���������ʮ����ֵ(0~1)
%% �м䲿�ֵĴ���ʹ�ָ�ʱ���������Ч��Ҳ����ֱ����xi = 1:bits��Ȼ��argu(ii) =  x(sum(para(1:(ii-1),1))+1:sum(para(1:ii,1)))*weight2'/(2^para(ii,1)-1);
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