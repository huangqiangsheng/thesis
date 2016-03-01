function Loss = Absorption2dB(alpha,Length,To,C)
    Loss = To*alpha*Length*10*log(exp(1))+C;
end