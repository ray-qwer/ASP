function [e_RLS, w] = ASP_RLS_1(xn_input,dn_input, delta, lamda)
M = 10;
L = length(xn_input);

w = zeros(M,L); %10 0~499
e_RLS = zeros(1,L);
e_RLS(1) = dn_input(1);

P_before = eye(M)/delta;
for n = 2:L
    xn = zeros(M,1);
    for m = n:-1:max(1,n-M+1)
        xn(n-m+1) = xn_input(m); 
    end
    k = (1/lamda*P_before*xn)/(1+1/lamda*xn'*P_before*xn);
    ee = dn_input(n) - w(:,n-1)'*xn;
    w(:,n) = w(:,n-1) + k*conj(ee);
    P_before = P_before/lamda - 1/lamda*k*xn'*P_before;

    e_RLS(n) = dn_input(n) - w(:,n)'*xn;
end
end