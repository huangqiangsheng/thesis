clear
well = InGaAsP_params(0.8,0.8,1550);
barrier = InGaAsP_params(0.786,0.466,1550);

MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 100e-10;
MQW.tw = 55e-10;
MQW.num_w = 1;
MQW.num_b = 2;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

save CC MQW;