r = [3;2;1];


% check dimension R

% initial
% M = length(r);
% P = zeros(M, 1);    kappa = zeros(M-1, 1); delta = zeros(M, 1);
% a = cell(M,1); 
% a{1} = 1; P(1) = r(1); delta(1) = conj(r(2));
% for i = 2:M
%     kappa(i-1) = - delta(i-1)/ P(i-1);
%     P(i) = P(i-1)* (1-abs(kappa(i-1)^2));
%     a{i} = [a{i-1}; 0]+ kappa(i-1).*[0; conj(flip(a{i-1}))];
%     if i < M
%         delta(i) = flip(r(2:i+1)).' * a{i};
%     end
% end
% a = a(2:end);
[a, P, kappa] = ASP_Levinson_Durbin(r);
