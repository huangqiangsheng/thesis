 function ConstPara = init_const()
    ConstPara.mu0 = 12.5663706144e-7;
    ConstPara.epsr0 = 8.854187818e-12;
    ConstPara.c0 = 1/sqrt(ConstPara.mu0*ConstPara.epsr0);
    ConstPara.h = 6.62607e-34;
    ConstPara.hbar =  ConstPara.h/(2*pi);
    ConstPara.q = 1.6021892e-19;
    ConstPara.m0 = 9.109534e-31;
    ConstPara.k = 8.617347e-5; %eV/K
 end