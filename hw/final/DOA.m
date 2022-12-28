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

figure(2);
fft_s = fft(theta_s_noisy);
fft_s( abs(fft_s) < max(abs(fft_s))/70)=0;
% fft_s = fftshift(fft_s);
% fft_s = fftshift(fft_s);
subplot(2,1,1);
plot(1:length(fft_s), fft_s);

theta_s_ifft = ifft(fft_s);
subplot(2,1,2);
plot(t, theta_s_ifft);
hold on;
% subplot(3,1,3);
plot(t, theta_s_noisy);
hold off;

figure(3);
fft_i = fft(theta_i_noisy);
fft_i( abs(fft_i) < max(abs(fft_i))/70)=0;
% fft_s = fftshift(fft_s);
% fft_s = fftshift(fft_s);
subplot(3,1,1);
plot(1:length(fft_i), fft_i);

theta_i_ifft = ifft(fft_i);
subplot(3,1,2);
plot(t, theta_i_ifft);
hold on;
plot(t, theta_i_noisy);
hold off;

figure(4);
[~,w]=ASP_RLS(theta_i_noisy, theta_i_ifft,0.01, 0.9,10);
theta_i_hat = zeros(size(theta_i_noisy));
tmp = [ones(1,9)*theta_i_noisy(1) theta_i_noisy];
for i = 1:length(theta_i_hat)
    theta_i_hat(i) = w(:,i)'* flip(tmp(i:i+9)).';
end
subplot(211);
plot(t, theta_i_hat);
hold on;
plot(t, theta_i_noisy);
plot(t, theta_i_ifft);
hold off;
subplot(212);
plot(t,theta_i_hat);
hold on;
plot(t, theta_i_ifft);
hold off;

IMF_set = hht(theta_i_noisy, t, 4);
figure(5)
for i = 1:size(IMF_set,1)
    subplot(size(IMF_set,1),1,i);
    plot(t, IMF_set(i,:));
end
figure(6);
recon_s = sum(IMF_set(2:end,:),1);
plot(t, theta_i_noisy);
hold on;
plot(t,recon_s);
hold off;