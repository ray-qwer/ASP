clear;
H = @(f) (1-exp(-1j*2*pi*f)/5)./((1-exp(-1j*2*pi*f)/2).*(1+exp(-1j*2*pi*f)/3));
Sx = @(f) abs(H(f)).^2;
M = 10;
Rx = zeros(M,M);
for k = 0:M-1
    Sxk =@(f) Sx(f).*exp(1j*2*pi*f*k);
    Rx(1,k+1) = integral(Sxk, -1/2, 1/2);
end
% autocorrelation R
for m = 2:M
    for n = 1:M
        k = m-n;
        if k<0 
            Rx(m,n) = conj(Rx(1,abs(k)+1));
        else 
            Rx(m,n) = Rx(1,k+1);
        end
    end
end

% wopt
p = zeros(M,1);
for k = 0:-1:-M+1
    Hk = @(f) H(f).*exp(1j*2*pi*f*k);
    p(abs(k)+1) = integral(Hk, -1/2, 1/2);
end
% for k = 0:1:M-1
%     Hk = @(f) H(f).*exp(1j*2*pi*f*k);
%     p(k+1) = integral(Hk, -1/2, 1/2);
% end
% p(2:end) = conj(p(2:end));

wopt = Rx\p;

% v, & construct x
v =load("ASP_HW2_Problem_5.mat","matV");
v = v.matV;
[R, L] = size(v);
x = zeros(R, L);
x(:,1) = v(:,1);
x(:,2) = v(:,2) -v(:,1)/5+ x(:,1)/6;
for k = 3:L
    x(:,k) = v(:,k) - v(:,k-1)/5 + x(:,k-1)/6 + x(:,k-2)/6;
end

% a, R = 5 
R = 5;

% LMS mu = 0.05
mu = 0.05;
[e_LMS_1, d_LMS_1] = getEnW(x(1:R,:), v(1:R,:), mu, M, R, 0, 0, wopt, "LMS");
e_LMS_1 = e_LMS_1/ R;
d_LMS_1 = d_LMS_1/ R;

% LMS mu = 0.1
mu = 0.1;
[e_LMS_2, d_LMS_2] = getEnW(x(1:R,:), v(1:R,:), mu, M, R, 0, 0, wopt, "LMS");
e_LMS_2 = e_LMS_2/ R;
d_LMS_2 = d_LMS_2/ R;

% NLMS mu_hat = 0.8
mu_hat = 0.8;
[e_NLMS_1, d_NLMS_1] = getEnW(x(1:R,:), v(1:R,:), mu_hat, M, R, 0, 0, wopt, "NLMS");
e_NLMS_1 = e_NLMS_1/ R;
d_NLMS_1 = d_NLMS_1/ R;

% NLMS mu_hat = 0.8
mu_hat = 0.9;
[e_NLMS_2, d_NLMS_2] = getEnW(x(1:R,:), v(1:R,:), mu_hat, M, R, 0, 0, wopt, "NLMS");
e_NLMS_2 = e_NLMS_2/ R;
d_NLMS_2 = d_NLMS_2/ R;

% RLS delta = 0.01 lambda = 0.5
delta = 0.01; lambda = 0.5;
[e_RLS_1, d_RLS_1] = getEnW(x(1:R,:), v(1:R,:), 0, M, R, lambda, delta, wopt, "RLS");
e_RLS_1 = e_RLS_1/ R;
d_RLS_1 = d_RLS_1/ R;

% RLS delta = 0.01 lambda = 0.8
delta = 0.01; lambda = 0.8;
[e_RLS_2, d_RLS_2] = getEnW(x(1:R,:), v(1:R,:), 0, M, R, lambda, delta, wopt, "RLS");
e_RLS_2 = e_RLS_2/ R;
d_RLS_2 = d_RLS_2/ R;

% RLS delta = 0.01 lambda = 0.9
delta = 0.01; lambda = 0.9;
[e_RLS_3, d_RLS_3] = getEnW(x(1:R,:), v(1:R,:), 0, M, R, lambda, delta, wopt, "RLS");
e_RLS_3 = e_RLS_3/ R;
d_RLS_3 = d_RLS_3/ R;

t = 1:L;
figure(1);
sgtitle("Problem 5(a)");
subplot(3,2,1);
plot(t, e_LMS_1);
hold on;
title("MSE of LMS");
plot(t, e_LMS_2);
legend("\mu = 0.05", "\mu = 0.1");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
hold off;
subplot(3,2,2);
plot(t, d_LMS_1);
hold on;
plot(t, d_LMS_2);
title("MSD of LMS");
legend("\mu = 0.05", "\mu = 0.1");
xlabel("n");
ylabel("$\hat{D}(n)$","Interpreter","latex");
hold off;
subplot(3,2,3);
plot(t, e_NLMS_1);
hold on;
plot(t, e_NLMS_2);
title("MSE of NLMS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\tilde{\mu} = 0.8$$', '$$\tilde{\mu} = 0.9$$','Interpreter','latex');
hold off;
subplot(3,2,4);
plot(t, d_NLMS_1);
hold on;
plot(t, d_NLMS_2);
title("MSD of NLMS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\tilde{\mu} = 0.8$$', '$$\tilde{\mu} = 0.9$$','Interpreter','latex');
hold off;
subplot(3,2,5);
plot(t, e_RLS_1);
hold on;
plot(t, e_RLS_2);
plot(t, e_RLS_3);
title("MSE of RLS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\delta = 0.01, \lambda = 0.5$$', '$$\delta = 0.01, \lambda = 0.8$$','$$\delta = 0.01, \lambda = 0.9$$','Interpreter','latex');
hold off;
subplot(3,2,6);
plot(t, d_RLS_1);
hold on;
plot(t, d_RLS_2);
plot(t, d_RLS_2);
title("MSD of RLS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\delta = 0.01, \lambda = 0.5$$', '$$\delta = 0.01, \lambda = 0.8$$','$$\delta = 0.01, \lambda = 0.9$$','Interpreter','latex');
hold off;
% ===================================================================================================
% b, R = 1000
R = 1000;

% LMS mu = 0.05
mu = 0.05;
[e_LMS_1, d_LMS_1] = getEnW(x(1:R,:), v(1:R,:), mu, M, R, 0, 0, wopt, "LMS");
e_LMS_1 = e_LMS_1/ R;
d_LMS_1 = d_LMS_1/ R;

% LMS mu = 0.1
mu = 0.1;
[e_LMS_2, d_LMS_2] = getEnW(x(1:R,:), v(1:R,:), mu, M, R, 0, 0, wopt, "LMS");
e_LMS_2 = e_LMS_2/ R;
d_LMS_2 = d_LMS_2/ R;

% NLMS mu_hat = 0.8
mu_hat = 0.8;
[e_NLMS_1, d_NLMS_1] = getEnW(x(1:R,:), v(1:R,:), mu_hat, M, R, 0, 0, wopt, "NLMS");
e_NLMS_1 = e_NLMS_1/ R;
d_NLMS_1 = d_NLMS_1/ R;

% NLMS mu_hat = 0.8
mu_hat = 0.9;
[e_NLMS_2, d_NLMS_2] = getEnW(x(1:R,:), v(1:R,:), mu_hat, M, R, 0, 0, wopt, "NLMS");
e_NLMS_2 = e_NLMS_2/ R;
d_NLMS_2 = d_NLMS_2/ R;

% RLS delta = 0.01 lambda = 0.5
delta = 0.01; lambda = 0.5;
[e_RLS_1, d_RLS_1] = getEnW(x(1:R,:), v(1:R,:), 0, M, R, lambda, delta, wopt, "RLS");
e_RLS_1 = e_RLS_1/ R;
d_RLS_1 = d_RLS_1/ R;

% RLS delta = 0.01 lambda = 0.8
delta = 0.01; lambda = 0.8;
[e_RLS_2, d_RLS_2] = getEnW(x(1:R,:), v(1:R,:), 0, M, R, lambda, delta, wopt, "RLS");
e_RLS_2 = e_RLS_2/ R;
d_RLS_2 = d_RLS_2/ R;

% RLS delta = 0.01 lambda = 0.9
delta = 0.01; lambda = 0.9;
[e_RLS_3, d_RLS_3] = getEnW(x(1:R,:), v(1:R,:), 0, M, R, lambda, delta, wopt, "RLS");
e_RLS_3 = e_RLS_3/ R;
d_RLS_3 = d_RLS_3/ R;

t = 1:L;
figure(2);
sgtitle("Problem 5(b)");
subplot(3,2,1);
plot(t, e_LMS_1);
hold on;
title("MSE of LMS");
plot(t, e_LMS_2);
legend("\mu = 0.05", "\mu = 0.1");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
hold off;
subplot(3,2,2);
plot(t, d_LMS_1);
hold on;
plot(t, d_LMS_2);
title("MSD of LMS");
legend("\mu = 0.05", "\mu = 0.1");
xlabel("n");
ylabel("$\hat{D}(n)$","Interpreter","latex");
hold off;
subplot(3,2,3);
plot(t, e_NLMS_1);
hold on;
plot(t, e_NLMS_2);
title("MSE of NLMS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\tilde{\mu} = 0.8$$', '$$\tilde{\mu} = 0.9$$','Interpreter','latex');
hold off;
subplot(3,2,4);
plot(t, d_NLMS_1);
hold on;
plot(t, d_NLMS_2);
title("MSD of NLMS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\tilde{\mu} = 0.8$$', '$$\tilde{\mu} = 0.9$$','Interpreter','latex');
hold off;
subplot(3,2,5);
plot(t, e_RLS_1);
hold on;
plot(t, e_RLS_2);
plot(t, e_RLS_3);
title("MSE of RLS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\delta = 0.01, \lambda = 0.5$$', '$$\delta = 0.01, \lambda = 0.8$$','$$\delta = 0.01, \lambda = 0.9$$','Interpreter','latex');
hold off;
subplot(3,2,6);
plot(t, d_RLS_1);
hold on;
plot(t, d_RLS_2);
plot(t, d_RLS_2);
title("MSD of RLS");
xlabel("n");
ylabel("$\hat{J}(n)$","Interpreter","latex");
legend('$$\delta = 0.01, \lambda = 0.5$$', '$$\delta = 0.01, \lambda = 0.8$$','$$\delta = 0.01, \lambda = 0.9$$','Interpreter','latex');
hold off;
% ===================================================================================================
% c, misadjustment
% Jex(n) = J(n) - Jmin
% M = 1/Jmin lim{n->inf} Jex
Mis = zeros(1,7);
Mis(1) = e_LMS_1(end) / min(e_LMS_1);
Mis(2) = e_LMS_2(end) / min(e_LMS_2);
Mis(3) = e_NLMS_1(end) / min(e_NLMS_1);
Mis(4) = e_NLMS_2(end) / min(e_NLMS_2);
Mis(5) = e_RLS_1(end) / min(e_RLS_1);
Mis(6) = e_RLS_2(end) / min(e_RLS_2);
Mis(7) = e_RLS_3(end) / min(e_RLS_3);
name_list = ["LMS,u= 0.05", "LMS, u= 0.1","NLMS,u= 0.8", "NLMS, u= 0.9", ...
    "RLS,lambda= 0.5, delta= 0.01", "RLS,lambda= 0.8, delta= 0.01","RLS,lambda= 0.9, delta= 0.01"];

for i = 1:length(name_list)
    disp(name_list(i)+": "+Mis(i));
end



function [e, d] = getEnW(x, v, mu, M, R, lambda, delta, wopt, filter)
    L = size(x,2);
    e = zeros(1, L);
    d = zeros(1, L);
    for i = 1:R
        if filter == "LMS"
            [e_i, w_i] = ASP_LMS(x(i,:), v(i,:), mu, M);
        elseif filter == "NLMS"
            [e_i, w_i] = ASP_NLMS(x(i,:), v(i,:), mu, M);
        elseif filter == "RLS"
            [e_i, w_i] = ASP_RLS(x(i,:), v(i,:), delta, lambda, M);
        else
            return;
        end
        e = e + abs(e_i).^2;
        d = d + sum(abs(wopt- w_i).^2, 1);
    end
    e = e/R; d = d/R;
end