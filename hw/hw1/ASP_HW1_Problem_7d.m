clear;
R = [3 0.7 0.3j; 0.7 3 0.7; -0.3j 0.7 3]; p = [6 -3j 1+3j].'; sd2 = 25;
w_x = linspace(-4,4,201);
w_y = linspace(-4,4,201);
[w0_x, w0_y] = meshgrid(w_x,w_y);
j_surf = zeros(size(w0_y));
for i = 1:size(j_surf(:))
    j_surf(i) = ASP_Wiener_MSE(R,[w0_x(i)+w0_y(i)*j ; 3+j ; 2+4j],p , sd2);
end

surf(w0_x,w0_y,real(j_surf),abs(j_surf));
colorbar;
xlabel("w0\{re\}");
ylabel("w0\{im\}");
zlabel("J");
title("ASP\_HW1\_Problem\_7d");
hold on;

syms w0R w0I;
wopt = [w0R+w0I*j;3+j;2+4j];
funcJ = ASP_Wiener_MSE(R, wopt, p, sd2);
g = gradient(funcJ, [w0R,w0I]);
[W0I,W0R] = solve(g==sym(0));
j_ans = subs(funcJ, {w0R,w0I},{W0R,W0I});

plot3(W0R,W0I,j_ans,'*r');
text(W0R,W0I,j_ans,sprintf("Re:%.2f\nIm:%.2f\nJ:%.2f",W0R,W0I,j_ans),'VerticalAlignment','top','HorizontalAlignment','left');
hold off;