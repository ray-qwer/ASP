% LCMV
matX = load("ASP_Final_Data.mat","matX");
theta_s_noisy = load("ASP_Final_Data.mat","theta_s_noisy");
theta_i_noisy = load("ASP_Final_Data.mat","theta_i_noisy");
matX = matX.matX;
theta_i_noisy = theta_i_noisy.theta_i_noisy;
theta_s_noisy = theta_s_noisy.theta_s_noisy;

fft_s = fft(theta_s_noisy);
fft_s( abs(fft_s) < max(abs(fft_s))/70)=0;
theta_s_denoise = ifft(fft_s);
fft_i = fft(theta_i_noisy);
fft_i( abs(fft_i) < max(abs(fft_i))/70)=0;
theta_i_denoise = ifft(fft_i);

g = [1 1E-5].';

delta = 0.01; lambda = 0.999;
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
t = [1:length(theta_s_noisy)];
figure(2);
subplot(211);
plot(t, real(y));
subplot(212);
plot(t, imag(y));