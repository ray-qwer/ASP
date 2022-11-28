function [e_LMS, w] = ASP_NLMS_1(xn_input,dn_input, mu)
M = 10;
L = length(xn_input);

% function input: x(n) 0~499 d(n) 0~499
% output: w 10 0~499 y 1
w = zeros(M,L); %10 0~499
e_LMS = zeros(1,L);

xn = zeros(M,1);

% n = 0
xn(1) = xn_input(1);
e_LMS(1) = dn_input(1);
for i = 2:L
    
    w(:,i) = w(:,i-1) + mu/(norm(xn,2)^2)*xn*conj(e_LMS(i-1));
    % next xn
    xn = zeros(M,1);
    for m = i:-1:max(1,i-M+1)
        xn(i-m+1) = xn_input(m); 
    end
    %next e_LMS
    e_LMS(i) = dn_input(i) - w(:,i)'*xn;
end
end