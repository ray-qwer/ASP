function [e_RLS, w_RLS] = ASP_RLS(x_input, d_input, delta, lambda, M)
    L = length(x_input);
    w_RLS = zeros(M, L);
    e_RLS = zeros(1, L);

    % initj
    e_RLS(1) = d_input(1);
    data = [zeros(1, M-1), x_input];
    P = eye(M)/delta;
    % iterate
    for m = 2:L
        x = flip(data(m:m+M-1));
        k = (1/lambda .* P * (x.'))/(1 + 1/lambda .* conj(x) * P * (x.')); % size k = (M, 1)
        e_RLS(m) = d_input(m)- x*conj(w_RLS(:,m-1)); % size e_priori: 1*1
        w_RLS(:, m) = w_RLS(:, m-1) + k* conj(e_RLS(m));
        % update for next loop
        P = 1/lambda .* P - 1/lambda.*(k*conj(x)*P);
%         e_RLS = d_input(m) - x*conj(w_RLS(:,m));
    end
end