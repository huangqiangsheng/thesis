well = InAlGaAs_params(0.574,0.111,1550);
barrier = InAlGaAs_params(0.468,0.217,1550);

MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 50e-10;
MQW.tw = 80e-10;
MQW.num_w = 1;
MQW.num_b = 2;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

save InAlGaAs_MZM MQW;