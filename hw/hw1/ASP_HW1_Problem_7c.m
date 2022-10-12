clear;
R = [3 0.7 0.3j; 0.7 3 0.7; -0.3j 0.7 3]; p = [6 -3j 1+3j].'; sd2 = 25;
reWo = linspace(-4,4,201);
J = zeros(size(reWo));
for i = 1:length(reWo)
    J(i) = ASP_Wiener_MSE(R,[reWo(i);-1;1],p,sd2);
end

syms w0;
assume(w0,'real');
wopt = [w0;-1;1];
funcJ = wopt'*R*wopt - p'*wopt - wopt'*p + sd2;
dif = diff(funcJ, w0);
w_ans = (solve(dif == sym(0), w0));
j_ans = (subs(funcJ,w0,w_ans));

plot(reWo,J);
ylabel("J",'rotation',0,'HorizontalAlignment','right');
xlabel("Re\{w0\}");
title("ASP HW1 problem 7c");
hold on;
plot(w_ans,j_ans,'*r');
text(w_ans,j_ans,sprintf("w_0:%f\nj_{min}:%f",w_ans,j_ans),'VerticalAlignment','top','HorizontalAlignment','left');
hold off;