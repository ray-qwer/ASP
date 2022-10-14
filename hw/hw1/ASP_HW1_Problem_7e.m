clear;
R = [3 0.7 0.3j; 0.7 3 0.7; -0.3j 0.7 3]; p = [6 -3j 1+3j].'; sd2 = 25;
w0_r = linspace(-3,3,201);
w1_i = linspace(-3,3,201);
[w0_R,w1_I] = meshgrid(w0_r,w1_i);
j_surf = zeros(size(w0_R));
for i = 1:size(j_surf(:))
    w = [w0_R(i)+0.3j; -0.6+w1_I(i)*j; 0.5+1.6j ];
    j_surf(i) = ASP_Wiener_MSE(R,w,p,sd2);
end

contour(w0_R,w1_I,real(j_surf),[2,3,4,5,10,20,50,100],"ShowText","on");
xlabel("Re\{W0\}");
ylabel("Im\{W1\}",'rotation',0,'HorizontalAlignment','right')
title("ASP\_HW1\_Problem\_7e")
hold on;

syms W0_R W1_I;
w_syms = [W0_R+0.3j; -0.6+W1_I*j; 0.5+1.6j];
funcJ = ASP_Wiener_MSE(R,w_syms,p,sd2);
g = gradient(funcJ,[W0_R,W1_I]);
g_s = solve(g == sym(0));
W0R = g_s.W0_R; W1I = g_s.W1_I;
j_ans = subs(funcJ,{W0_R,W1_I},{W0R,W1I});
plot(W0R,W1I,'*r');
text(W0R,W1I,sprintf("Re w0:%.2f\nIm w1: %.2f\nJ: %.2f",W0R,W1I,j_ans),'VerticalAlignment','top','HorizontalAlignment','left');
hold off;
