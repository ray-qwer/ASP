% R = [3 0.7 0.3j; 0.7 3 0.7; -0.3j 0.7 3]; p = [6 -3j 1+3j]'; sd2 = 25;
% wopt = R\p;
% jmin = ASP_Wiener_MSE(R,wopt,p,sd2);

function J = ASP_Wiener_MSE(R, w, p ,sd2)
    
    if ~isreal(sd2) || sd2 < 0
        error("sd2 is negative or complex-valued.")
    end
    if size(R,1) ~= size(R,2) || size(R,1) ~= length(w) || length(w) ~= length(p)
        error("The dimensions of these input arguments are not suitable.")
    end
    [V,D] = eig(R);
    if ~isequal(round(abs(V*V')), diag(ones(1,length(R)))) || any(D<0,"all")
        error("R is not positive semidefinite.")
    end

    J = w'*R*w - p'*w - w'*p + sd2;
end