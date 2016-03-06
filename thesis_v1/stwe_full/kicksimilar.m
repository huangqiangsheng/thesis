function state = kicksimilar(state,simiRange,Ntemp)
% 用于GA 置于gaoutputfcn或者衍生函数，用于剔除种群中的部分相似因子。
%state GA状态结构
%simiRange, 排斥范围bits
%Ntemp, 操作的模板数, 最优值前几位作为模板

% size(state.Score)%%% np*1
% size(state.Population(1,:))%%%1*bits
nP = length(state.Score);
Bits = length(state.Population(1,:));

[unused,order] = sort(state.Score);
mm=0;
for jj = 1:Ntemp
tempPopu  = state.Population(order(jj+mm),:);%%%提出第一个模板的相似个体
point     = order(jj+mm);
tempScore = state.Score(point);
if tempScore ==0  %%%防止相互删的情况, 如1和2相似,则1会删2,若无处理,2会删1,导致丢失最优
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
    