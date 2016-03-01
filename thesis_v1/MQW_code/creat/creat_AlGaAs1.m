well.Ec = 0;
well.EvH = 0;
well.EvL = 0;
well.me = 0.067;
well.mHH = 0.51;
well.mLH = 0.085;

barrier.Ec = 0.4;
barrier.EvH = -0.071;
barrier.EvL = -0.071;
barrier.me = 0.0985;
barrier.mHH = 0.601;
barrier.mLH = 0.109;

MQW.barrier = barrier;
MQW.well = well;
MQW.tb = 35e-10;
MQW.tw = 30e-10;
MQW.num_w = 1;
MQW.num_b = 2;
MQW.t = MQW.num_b*MQW.tb+MQW.num_w*MQW.tw;
MQW.period = MQW.tb + MQW.tw;

save AlGaAs MQW;