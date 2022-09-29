k = [-10:10];
fset = [-1:0.001:1];
rx = 1/81*(-1/3).^abs(k);
fs = fftshift(fft(rx));
