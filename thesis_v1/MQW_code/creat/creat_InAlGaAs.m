barrier = InAlGaAs_params(0.47,0.2,1550);
well = InAlGaAs_params(0.59,0.08,1550);
SCH = InAlGaAs_params(0.52,0.16,1550);
pcontact = InAlGaAs_params(0.53,0,1550);
MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 70e-10;
MQW.tw = 110e-10;
MQW.num_w = 1;
MQW.num_b = 2;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

% save InAlGaAs MQW;