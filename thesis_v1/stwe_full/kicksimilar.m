function state = kicksimilar(state,simiRange,Ntemp)
% ����GA ����gaoutputfcn�������������������޳���Ⱥ�еĲ����������ӡ�
%state GA״̬�ṹ
%simiRange, �ųⷶΧbits
%Ntemp, ������ģ����, ����ֵǰ��λ��Ϊģ��

% size(state.Score)%%% np*1
% size(state.Population(1,:))%%%1*bits
nP = length(state.Score);
Bits = length(state.Population(1,:));

[unused,order] = sort(state.Score);
mm=0;
for jj = 1:Ntemp
tempPopu  = state.Population(order(jj+mm),:);%%%�����һ��ģ������Ƹ���
point     = order(jj+mm);
tempScore = state.Score(point);
if tempScore ==0  %%%��ֹ�໥ɾ�����, ��1��2����,��1��ɾ2,���޴���,2��ɾ1,���¶�ʧ����
    mm=mm+1;
	continue;
end
simiFlag = zeros(size(tempPopu));
    for ii = 1:nP
        simiFlag = (tempPopu == state.Population(ii,:) );
        R = sum(simiFlag);
        if R>Bits-simiRange 
            state.Score(ii) = 0;
        end
    end
state.Score(point) = tempScore;
end
if mm>Ntemp/2
   str = sprintf('%-5.0f    %2d\n',...
                   state.Generation,mm);
  % disp(str);
   end
end
    