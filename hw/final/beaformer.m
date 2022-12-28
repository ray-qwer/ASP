% LCMV
matX = load("ASP_Final_Data.mat","matX");
theta_s_noisy = load("ASP_Final_Data.mat","theta_s_noisy");
theta_i_noisy = load("ASP_Final_Data.mat","theta_i_noisy");
matX = matX.matX;
theta_i_noisy = theta_i_noisy.theta_i_noisy;
theta_s_noisy = theta_s_noisy.theta_s_noisy;
t = [1:length(theta_s_noisy)];

threshold = 0.006;
IMF_set_i = hht(theta_i_noisy, t, 4);
mean_set_i = mean(IMF_set_i,2);
start_index_i = find(abs(mean_set_i) > threshold,1);
theta_i_denoise = sum(IMF_set_i(start_index_i:end,:),1);

IMF_set_s = hht(theta_s_noisy, t, 3);
mean_set_s = mean(IMF_set_s,2);
start_index_s = find(abs(mean_set_s) > threshold,1);
theta_s_denoise = sum(IMF_set_s(start_index_s:end,:),1);
figure(1);
plot(t,theta_s_denoise);
hold on;
plot(t,theta_i_denoise);
hold off;
title("denoised DOA");
legend('$$\hat{\theta _s}(t)$$','$$\hat{\theta _i}(t)$$','Interpreter','latex');
xlabel("time");
ylabel("angle(degree)");

g = [1 0.002].';

delta = 0.01; lambda = 0.995;
w = zeros(size(matX));
P = 1/delta*eye(size(matX,1));
for i = 1:length(theta_s_noisy)
    C = [exp(1j*pi*(([0:9].'))*sin(theta_s_denoise(i)/180*pi)), exp(1j*pi*(([0:9].'))*sin(theta_i_denoise(i)/180*pi))];
    x_hat = matX(:,i);
    k = (1/lambda*P*x_hat)/(1+1/lambda*((x_hat')*P*x_hat));
    P = 1/lambda.*P - 1/lambda*k*x_hat'*P;
    w(:,i) = P*C*((C'*P*C)\g);
end

y = zeros(size(theta_s_noisy));
for i = 1:length(theta_s_noisy)
    y(i) = w(:,i)' * matX(:,i);
end

figure(2);
% title("estimated source signal");
subplot(211);
plot(t, real(y));
title('real part of $$\hat{s}(t)$$','Interpreter','latex');
xlabel("time");
ylabel("magnitude");
subplot(212);
plot(t, imag(y));
title('imag part of $$\hat{s}(t)$$','Interpreter','latex');
xlabel("time");
ylabel("magnitude");
sgtitle('estimated source signal $$\hat{s}(t)$$','Interpreter','Latex');