well = InAlGaAs_params(0.65,0.09,1550);
barrier = InAlGaAs_params(0.42,0.17,1550);

MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 70e-10;
MQW.tw = 110e-10;
MQW.num_w = 1;
MQW.num_b = 2;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

save InAlGaAs_MQW_new_2 MQW;