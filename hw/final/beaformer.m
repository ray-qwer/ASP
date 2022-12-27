matX = load("ASP_Final_Data.mat","matX");
theta_s_noisy = load("ASP_Final_Data.mat","theta_s_noisy");
theta_i_noisy = load("ASP_Final_Data.mat","theta_i_noisy");
matX = matX.matX;
theta_i_noisy = theta_i_noisy.theta_i_noisy;
theta_s_noisy = theta_s_noisy.theta_s_noisy;
t=[1:length(theta_i_noisy)];
figure(1);
subplot(2,1,1);
plot(t,theta_s_noisy);
subplot(2,1,2);
plot(t,theta_i_noisy);

[~, w] = ASP_RLS(theta_s_noisy, theta_s_noisy, 0.01, 0.5, 10);
w = w(:,end);
theta_s_tmp = ([zeros(1,9) theta_s_noisy]);
theta_s_hat = zeros(size(theta_s_noisy));
for i = 1:size(theta_s_noisy,2)
    x = flip(theta_s_tmp(i:i+9));
    theta_s_hat(i) = x*conj(w);
end
figure(2);
subplot(2,1,1);
plot(t, theta_s_noisy);
subplot(2,1,2);
plot(t, theta_s_hat);