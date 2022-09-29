k = [-100:100];
fset = [-1:0.001:1];
rx = (-1/3).^abs(k);
s = zeros(size(fset));
for i = 1:length(s)
    tmp = exp(-j*2*pi*fset(i).*k);
    s(i) = 1/81*sum(tmp.*rx);
end
subplot(2,1,1);
plot(fset,abs(s));
subplot(2,1,2);
plot(fset,angle(s));