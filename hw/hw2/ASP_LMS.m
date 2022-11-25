% M = 10; w = zeros(M,1);
% matV = load("ASP_HW2_Problem_5.mat","matV");
% v = matV.matV;
% [R, L] = size(v);
% x = zeros(R, L);
% x(:,1) = v(:,1);
% x(:,2) = v(:,2) -v(:,1)/5+ x(:,1)/6;
% for k = 3:L
%     x(:,k) = v(:,k) - v(:,k-1)/5 + x(:,k-1)/6 + x(:,k-2)/6;
% end
% 
% r = 5;
% e = zeros(1,L);
% mu = 0.02;
% data = [zeros(R, M-1), x];
% w = zeros(M, L);
% % t = 0
% e(1) = v(1, 1);
% for m = 2:L
% %     e(m) = v(1,m) - flip(data(1,m:m+M-1))*conj(w(:,m));
%     w(:,m) = w(:,m-1) + mu*conj(e(m-1)).*(flip(data(1,m-1:m+M-2)).');
%     e(m) = v(1, m) - flip(data(1,m:m+M-1))*conj(w(:,m));
% end
function [e_LMS, w_LMS] = ASP_LMS(x_input, d_input, mu, M)
    L = length(x_input);
    w_LMS = zeros(M, L);
    e_LMS = zeros(1, L);

    % init
    e_LMS(1) = d_input(1);
    data = [zeros(1, M-1), x_input];
    x = flip(data(1:M));
    % iterate
    for m = 2:L
        w_LMS(:, m) = w_LMS(:, m-1) + mu* conj(e_LMS(m-1)).* (x.');
        % update for next loop
        x = flip(data(m:m+M-1));
        e_LMS(m) = d_input(m) - x*conj(w_LMS(:,m));
    end
end