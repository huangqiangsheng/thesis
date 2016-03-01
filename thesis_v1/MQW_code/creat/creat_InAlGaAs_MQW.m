well = InAlGaAs_params(0.59,0.08,1550);
barrier = InAlGaAs_params(0.47,0.20,1550);

MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 70e-10;
MQW.tw = 110e-10;
MQW.num_w = 10;
MQW.num_b = 11;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

save InAlGaAs_MQW MQW;