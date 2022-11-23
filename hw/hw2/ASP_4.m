syms z n;
H(z) = (1-z^(-1)/5)/((1-z^(-1)/2)*(1+z^(-1)/3));
% h(z) = 1/(1+3/4*z^-1);
Sx = H(z)*conj(H(1/conj(z)));
rx(n) = iztrans(Sx);
R = zeros(10);
for i = 1:10
%     for j = 1:10
%         R(i,j) = rx(abs(i-j));
%     end
    R(i,:) = rx(-i+1:-i+10);
end

e = eig(R);