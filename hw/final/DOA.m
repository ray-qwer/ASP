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
% 
% [~, w] = ASP_RLS(theta_s_noisy, theta_s_noisy, 0.01, 0.5, 10);
% w = w(:,end);
% theta_s_tmp = ([zeros(1,9) theta_s_noisy]);
% theta_s_hat = zeros(size(theta_s_noisy));
% for i = 1:size(theta_s_noisy,2)
%     x = flip(theta_s_tmp(i:i+9));
%     theta_s_hat(i) = x*conj(w);
% end


IMF_set_i = hht(theta_i_noisy, t, 4);
figure(2)
for i = 1:size(IMF_set_i,1)
    subplot(size(IMF_set_i,1),1,i);
    plot(t, IMF_set_i(i,:));
end
figure(3);
recon_i = sum(IMF_set_i(2:end,:),1);
plot(t, theta_i_noisy);
hold on;
plot(t,recon_i);
hold off;

IMF_set_s = hht(theta_s_noisy, t, 3);
figure(4)
for i = 1:size(IMF_set_s,1)
    subplot(size(IMF_set_s,1),1,i);
    plot(t, IMF_set_s(i,:));
end
figure(5);
recon_s = sum(IMF_set_s(3:end,:),1);
plot(t, theta_s_noisy);
hold on;
plot(t,recon_s);
hold off;

figure(6);
plot(t,recon_i);
hold on;
plot(t,recon_s);
hold off;

figure(7);
[~,w]=ASP_RLS(theta_i_noisy, recon_i,0.01, 0.9,10);
theta_i_hat = zeros(size(theta_i_noisy));
tmp = [ones(1,9)*theta_i_noisy(1) theta_i_noisy];
for i = 1:length(theta_i_hat)
    theta_i_hat(i) = w(:,i)'* flip(tmp(i:i+9)).';
end

plot(t, theta_i_noisy);
hold on;
plot(t, theta_i_hat);
hold off;

figure(8);
[~,w]=ASP_RLS(theta_s_noisy, recon_s,0.01, 0.9,10);
theta_s_hat = zeros(size(theta_s_noisy));
tmp = [ones(1,9)*theta_s_noisy(1) theta_s_noisy];
for i = 1:length(theta_s_hat)
    theta_s_hat(i) = w(:,i)'* flip(tmp(i:i+9)).';
end
plot(t, theta_s_noisy);
hold on;
plot(t, theta_s_hat);
hold off;

figure(9);
plot(t, theta_s_hat);
hold on;
plot(t, recon_s);
plot(t, recon_i);
plot(t, theta_i_hat);
hold off;
legend("theta_s_hat","recon_s","recon_i","theta_i_hat");