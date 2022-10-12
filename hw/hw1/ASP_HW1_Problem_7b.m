clear;
R = [3 0.7 0.3j; 0.7 3 0.7; -0.3j 0.7 3]; p = [6 -3j 1+3j].'; sd2 = 25;
wopt = R\p;
jmin = ASP_Wiener_MSE(R,wopt,p,sd2);