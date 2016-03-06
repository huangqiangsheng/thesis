function EyeDiagramFilter(Response)
% ��������������ת��Ϊ������Ӧ��ͨ��filter����������PRBS���ε������ͨ���ص���������ͼ
% ��Ӧ��������Ƶ����Ӧ��������Ӧ�Ͳ�������������ݵĶԱȡ�
% ���ṩ�����Ƶ����Ӧ���ݡ�

% clear all;
% 
% % load response_data_165.mat
% load response_data_210_lewen.mat

freqNum = length(Response); % Ƶ����Ӧ ����

df = 1;
dfsp = df*1e9;
Tp = 1/dfsp;  % repetition rate  TpӦ�õ���TotalTime���ɴ˿���ȷ��Response����

taos = 0.5e-12; % the signal source's RC time constant
taod = 0.5e-12; % the detector's RC time constant
z_cn = find(abs(Response)<1.0e-10);
Response(z_cn)=0;% ǿ�Ƹ��������µ�С��Ϊ0��
Response = Response./(1+j*2*pi*(1:freqNum)*dfsp*taos)./(1+j*2*pi*(1:freqNum)*dfsp*taod);
endDrawNum=120;
% tt= 1./(1+j*2*pi*(1:freqNum)*dfsp*taos)./(1+j*2*pi*(1:freqNum)*dfsp*taos);
% 
% figure; 
% plot(df*(1:endDrawNum),20*log10(real(abs(tt(1:endDrawNum))))); grid on;
% x1=round(5/df);x2=round(49/df);
% Response(x1:x2)=Response(x1:x2)+(Response(x1)+(Response(x2)-Response(x1))/(x2-x1)*((x1:x2)-x1)-Response(x1:x2))*0.8;
% % x1=round(49/df);x2=round(85/df);
% % Response(x1:x2)=Response(x1:x2)+(Response(x1)+(Response(x2)-Response(x1))/(x2-x1)*((x1:x2)-x1)-Response(x1:x2))*1;


% figure; subplot(2,1,1);plot(df*(1:endDrawNum),20*log10(abs(Response(1:endDrawNum)))); grid on;
% subplot(2,1,2);plot(df*(1:endDrawNum),angle(Response(1:endDrawNum))*180/pi); grid on;
 
N=2^(nextpow2(freqNum));
ts = Tp/N;
pulseResponse = real(ifft(Response,N));
pulseResponse = fftshift(pulseResponse);
temp=pulseResponse(1:N/2-200);
pulseResponse(1:N/2-200)=[];
pulseResponse=[pulseResponse,temp];
p_cn = find(abs(pulseResponse)<2.0e-4);
smalldc = sum(pulseResponse(p_cn))/length(p_cn);% ǿ�Ƹ��������µ�С��Ϊ0��
pulseResponse = pulseResponse - smalldc;
% x1=round(58.84/ts/1e12);x2=round(66.65/ts/1e12); 
% pulseResponse(x1:x2)=pulseResponse(x1:x2)+(pulseResponse(x1)+(pulseResponse(x2)-pulseResponse(x1))/(x2-x1)*((x1:x2)-x1)-pulseResponse(x1:x2))*1;
% x1=round(61.77/ts/1e12);x2=round(71.29/ts/1e12); 
% pulseResponse(x1:x2)=pulseResponse(x1:x2)+(pulseResponse(x1)+(pulseResponse(x2)-pulseResponse(x1))/(x2-x1)*((x1:x2)-x1)-pulseResponse(x1:x2))*(-1);
% x1=round(66.65/ts/1e12);x2=round(79.59/ts/1e12); 
% pulseResponse(x1:x2)=pulseResponse(x1:x2)+(pulseResponse(x1)+(pulseResponse(x2)-pulseResponse(x1))/(x2-x1)*((x1:x2)-x1)-pulseResponse(x1:x2))*1.0;
% x1=round(58.84/ts/1e12);x2=round(79.59/ts/1e12); 
% pulseResponse(x1:x2)=pulseResponse(x1:x2)+(pulseResponse(x1)+(pulseResponse(x2)-pulseResponse(x1))/(x2-x1)*((x1:x2)-x1)-pulseResponse(x1:x2))*0.65;

%  test = fft(pulseResponse);
% test = fftshift(test); test(1:2048)=[];
%  figure;plot(20*log10(abs(test(1:120)))-20*log10(abs(test(1))));grid on
 
% figure; plot(ts*(1:1201)*1e12,real(pulseResponse(1:1201))); grid on;
figure; plot(ts*(1:1201)*1e12-40,real(pulseResponse(1:1201))/max(real(pulseResponse(1:1201)))); grid on;
xlabel('Time (ps)'); ylabel('Normalized impulse response');
%legend('IFFT Result'); 

IFFTImpulseResp = pulseResponse/ts;

fbconnection =13; %2^fbconnection-1
%׼��������ɢ��
BitNum = 2^fbconnection+3*fbconnection;

OverSampleNum = 60;% ÿbit�ظ��������100->40Gbps, 40->102.5Gbps;25->~160Gbps;34->120Gbps; 82->50Gbps
SampleNum = OverSampleNum * BitNum;% �ܲ���������
SampleFreq=1/ts; % ʵ�ʲ���Ƶ�ʡ�
LEPfreq = 1/(ts*OverSampleNum)% ��ͼƵ��
NyquestFre=SampleFreq/2; % NyquestƵ��

mseq = PRBSgenerator(fbconnection,BitNum);%ʹ��α����롣�����������������ͷ��ʼ�ظ������㣬��ǰ����ȡ��
InputBit = 2*mseq(1: BitNum)-1;% ȥ��ֱ������
InputSeq = repmat(InputBit, [OverSampleNum, 1]); % ÿbit�ظ�չ������ɢΪ���������С�
clear InputBit
InputSeq = InputSeq(:);
% ������ɢֵ׼����ϡ�
clear mseq 


OutputSeq = ts * filter(IFFTImpulseResp',1,InputSeq);


% figure
% subplot(2,1,1)
% plotNum = OverSampleNum*200;
% plot(1:plotNum,InputSeq(1:plotNum),'r');
% subplot(2,1,2)
% plot(1:plotNum,real(OutputSeq(1:plotNum)),'g');

%%% draw eye diagram
OutputSignal = real(OutputSeq);
PlotNum = floor(1.3* OverSampleNum);
startPlotNum = OverSampleNum*10+50;% shift the diagram
endPlotNum = startPlotNum+2^fbconnection+fbconnection+PlotNum;
%k = 0; 
k = startPlotNum;
m = linspace(-0.5*PlotNum*ts*1e12, 0.5*PlotNum*ts*1e12, PlotNum);
figure
subplot(1,1,1);hold on;
axis([ -0.55 0.55 min(m) max(m)]);
while k+PlotNum <= endPlotNum %%length(OutputSignal)
	plot(m,OutputSignal(k+1:k+PlotNum));
    k = k+OverSampleNum; 
end
hold off; 
%title('Eye Diagram of the 100 Gbps Output Signal', 'fonts', 12);
ylabel('Amplitude'); xlabel('Time (ps)'); axis([-inf, inf, -inf, inf]); hold off
