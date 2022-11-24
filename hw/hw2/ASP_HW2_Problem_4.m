clear;
H = @(f) (1-exp(-1j*2*pi*f)/5)./((1-exp(-1j*2*pi*f)/2).*(1+exp(-1j*2*pi*f)/3));
Sx = @(f) abs(H(f)).^2;
M = 10;
u = 0.02;
R = zeros(M,M);
for k = 0:M-1
    Sxk =@(f) Sx(f).*exp(1j*2*pi*f*k);
    R(1,k+1) = integral(Sxk, -1/2, 1/2);
end

for m = 2:M
    for n = 1:M
        k = m-n;
        if k<0 
            R(m,n) = conj(R(1,abs(k)+1));
        else 
            R(m,n) = R(1,k+1);
        end
    end
end

sd2 = 1;

p = zeros(M,1);
for k = 0:-1:-M+1
    Hk = @(f) H(f).*exp(1j*2*pi*f*k);
    p(abs(k)+1) = integral(Hk, -1/2, 1/2);
end

wopt = R\p;
disp("i:");
disp("wopt= [");
for k = 1:length(wopt)
    disp(wopt(k));
end
disp("]");
e = eig(R);
disp("ii:");
disp(max(e)-min(e));

disp("iii:");
disp(length(e));
disp("iv");
for k = 1:length(e)
    words = "";
    t_c = 0;
    vk = 1-u*e(k);
    if  vk< -1 || vk > 1
        words = "unstable";
    elseif vk < 0
        words = "underdamped";
    else
        words = "overdamped";
        t_c = -1/log(vk);
    end
    if t_c~=0
        disp("v"+k+": "+vk+" "+words+", time constant: "+t_c);
    else
        disp("v"+k+": "+vk+" "+words);
    end
       
end
