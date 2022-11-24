clear;
H = @(f) (1-exp(-1j*2*pi*f)/5)./((1-exp(-1j*2*pi*f)/2).*(1+exp(-1j*2*pi*f)/3));
Sx = @(f) abs(H(f)).^2;
M = 10;
u = 0.02;
Rx = zeros(M,M);
for k = 0:M-1
    Sxk =@(f) Sx(f).*exp(1j*2*pi*f*k);
    Rx(1,k+1) = integral(Sxk, -1/2, 1/2);
end
% autocorrelation R
for m = 2:M
    for n = 1:M
        k = m-n;
        if k<0 
            Rx(m,n) = conj(Rx(1,abs(k)+1));
        else 
            Rx(m,n) = Rx(1,k+1);
        end
    end
end

% v, x 
v = load("ASP_HW2_Problem_5.mat","matV");
[R, L] = size(v);
x = zeros(R, L);
x(:,1) = v(:,1);
x(:,2) = v(:,2) -v(:,1)/5+ x(:,1)/6;
for k = 3:L
    x(:,k) = v(:,k) - v(:,k-1)/5 + x(:,k-1)/6 + x(:,k-2)/6;
end

[MSD, MSE] = ASP_LMS();