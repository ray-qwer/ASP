function [e_NLMS, w_NLMS] = ASP_NLMS(x_input, d_input, mu_hat, M)
    L = length(x_input);
    w_NLMS = zeros(M, L);
    e_NLMS = zeros(1, L);

    % init
    e_NLMS(1) = d_input(1);
    data = [zeros(1, M-1), x_input];
    x = flip(data(1:M));
    % iterate
    for m = 2:L
        mu = mu_hat/ (x*(x'));
        w_NLMS(:, m) = w_NLMS(:, m-1) + mu* conj(e_NLMS(m-1)).* (x.');
        % update for next loop
        x = flip(data(m:m+M-1));
        e_NLMS(m) = d_input(m) - x*conj(w_NLMS(:,m));
    end
end