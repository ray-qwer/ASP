data = load("ASP_HW3_Problem_4.mat");
C_n = data.C_n; F_n1_n = data.F_n1_n;
Q1_n = data.Q1_n; Q2_n = data.Q2_n;
Y_tilde = data.Y_tilde;

M = size(C_n, 2);
n_max  = size(Y_tilde, 2);
% initail state
x_n1_yn = zeros(M,1);
K_n1_n = eye(M);        % K(1,0)
x_n_yn_list = zeros(M,n_max);
t = 1:n_max;

% iteration
for i = 1:n_max
    % get Kalman Gain
    G = Kalman_Gain(C_n, Q2_n, F_n1_n, K_n1_n);
    % get new K(n+1, n)
    K_n1_n = Riccati(K_n1_n, F_n1_n, G, C_n, Q1_n);
    % get new x(n+1| yn)
    x_n1_yn = One_step_Prediction(Y_tilde(:,i), C_n, x_n1_yn, G, F_n1_n);
    x_n_yn = F_n1_n \ x_n1_yn;
    x_n_yn_list(:,i) = x_n_yn;
end

% plot
% real part
figure(1);

for i = 1:M
    subplot(M,1,i);
    plot(t,real(x_n_yn_list(i,:)));
    title(['$$\hat{x}_{',num2str(i),'}(n|\mathcal{Y}_n)$$'],'Interpreter','Latex');
    xlabel("n");
    ylabel("real part");
end
sgtitle("ASP\_HW3\_Problem\_4\_Real");

figure(2);
for i = 1:M
    subplot(M,1,i);
    plot(t,imag(x_n_yn_list(i,:)));
    title(['$$\hat{x}_{',num2str(i),'}(n|\mathcal{Y}_n)$$'],'Interpreter','Latex');
    xlabel("n");
    ylabel("imag part");
end
sgtitle("ASP\_HW3\_Problem\_4\_Imag");

figure(3);
for i = 1:M
    subplot(M,1,i);
    plot(t,abs(x_n_yn_list(i,:)));
    title(['$$\hat{x}_{',num2str(i),'}(n|\mathcal{Y}_n)$$'],'Interpreter','Latex');
    xlabel("n");
    ylabel("magnitude");
end
sgtitle("ASP\_HW3\_Problem\_4\_Mag");

phase = (atan2(imag(x_n_yn_list),real(x_n_yn_list)));
figure(4);
for i = 1:M
    subplot(M,1,i);
    plot(t,(unwrap(phase(i,:))));
    title(['$$\hat{x}_{',num2str(i),'}(n|\mathcal{Y}_n)$$'],'Interpreter','Latex');
    xlabel("n");
    ylabel("phase");
end
sgtitle("ASP\_HW3\_Problem\_4\_Phase");

% b: regression
phase = phase.';
t = t.';
regress = zeros(2,M);
for i = 1:M
    b = [ones(n_max,1), t]\unwrap(phase(:,i));
    regress(:,i) = b;
    disp(['The slope of phase{',int2str(i),'} is ',num2str(b(2)),'.']);
end

figure(5);
for i = 1:M
    subplot(M,1,i);
    plot(t,unwrap(phase(:,i)));
    hold on;
    plot(t, [ones(n_max,1),t]*regress(:,i));
    title(['$$\hat{x}_{',num2str(i),'}(n|\mathcal{Y}_n)$$'],'Interpreter','Latex');
    xlabel("n");
    ylabel("phase");
    hold off;
end
sgtitle("ASP\_HW3\_Problem\_4\_Phase & regression line");

function x_n1_yn = One_step_Prediction(y, C, x_n_y_1, G, F_n1_n)
    alpha = y - C*x_n_y_1;
    x_n1_yn = F_n1_n*x_n_y_1 + G* alpha;
end

function G = Kalman_Gain(C, Q2, F_n1_n, K_n_n_1)
    R = C*K_n_n_1 * C' +Q2;
    G = F_n1_n*K_n_n_1*C'*inv(R);
end

function K_n1_n = Riccati(K_n_n_1, F_n1_n, G, C, Q1)
    K_n = K_n_n_1 - inv(F_n1_n)*G*C*K_n_n_1;
    K_n1_n = F_n1_n * K_n *(F_n1_n') +Q1;
end