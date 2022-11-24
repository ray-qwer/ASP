M = 10; w = zeros(M,1);
matV = load("ASP_HW2_Problem_5.mat","matV");
v = matV.matV;
[R, L] = size(v);
x = zeros(R, L);
x(:,1) = v(:,1);
x(:,2) = v(:,2) -v(:,1)/5+ x(:,1)/6;
for k = 3:L
    x(:,k) = v(:,k) - v(:,k-1)/5 + x(:,k-1)/6 + x(:,k-2)/6;
end

r = 5;
e = zeros(1,5);
u = 0.02;
for k = 1:r
    w = zeros(M,1);
    data = [zeros(1,M-1), x(k,:)];
    for i = 1:L
        en = v(k,i) - flip(data(i:i+M-1))*conj(w); 
        w= w+u.*flip(data(i:i+M-1))'.*en;
    end
    e(k) = en;
end
