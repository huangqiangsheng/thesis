const = init_const();
well = InAlGaAs_params(0.59,0.08,1550);
well.tw = 110e-10;
dJ = (1.24/1500-1.24/(1550))*const.q;
N1 = const.k*const.q*300*log(2)/(well.tw*(pi*const.hbar^2/(well.me*const.m0)));
N2 = (dJ+const.k*const.q*300*log(2)*2)/(2*well.tw*(pi*const.hbar^2/(well.me*const.m0)));