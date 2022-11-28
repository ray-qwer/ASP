% % input
% clear; % M = 10
% data = load('ASP_HW3_Problem_4.mat');
% V = data.matV;
% V_R = size(V,1);
% V_L = size(V,2);
% 
% % determine x(n)
% % x(n)-1/4*x(n-1) - 1/8*x(n-2) = v(n) + 1/8*v(n-1)
% x = zeros(V_R, V_L);
% x(:,1) = V(:,1);
% x(:,2) = V(:,2) + V(:,1)/8 + x(:,1)/4;
% for i = 3:V_L
%     x(:,i) = V(:,i) + V(:,i-1)/8 + x(:,i-1)/4 + x(:,i-2)/8;
% end
% 
% mu = 0.02; %0.1
% [w, e_LMS] = ASP_LMS1(x(1,:), V(1,:), mu);
% figure(1)
% plot(1:500, abs(e_LMS));

% LMS filter
function [e_LMS,w ] = ASP_LMS_1(xn_input,dn_input, mu)
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
    
    w(:,i) = w(:,i-1) + mu*xn*conj(e_LMS(i-1));
    % next xn
    xn = zeros(M,1);
    for m = i:-1:max(1,i-M+1)
        xn(i-m+1) = xn_input(m); 
    end
    %next e_LMS
    e_LMS(i) = dn_input(i) - w(:,i)'*xn;
end
end