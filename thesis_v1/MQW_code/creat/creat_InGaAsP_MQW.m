%UCSB tw 12.8nm, tb 8nm
% well = InGaAsP_params(0.485,0.979,1550);
% barrier = InGaAsP_params(0.923,0.325,1550);
% v2 tw 12.8nm, tb 8nm
% well = InGaAsP_params(0.458,0.989,1550);
% barrier = InGaAsP_params(0.996,0.251,1550);
% v3 tw 117nm, tb 72nm
% well = InGaAsP_params(0.457,0.993,1550);
% barrier = InGaAsP_params(0.996,0.252,1550);
%v5 125.244 50.153 95.414 96.070 25.000 
% well = InGaAsP_params(0.502,0.954,1550);
% barrier = InGaAsP_params(0.961,0.25,1550);
%v6 101.126 78.044 73.768 74.393 34.694 
well = InGaAsP_params(0.780,0.737,1550);
barrier = InGaAsP_params(0.780,0.347,1550);
MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 80e-10/2;
MQW.tw = 101e-10;
MQW.num_w = 10;
MQW.num_b = 11;
MQW.t = MQW.num_b*MQW.tb*2+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

save InGaAsP_MQW MQW;