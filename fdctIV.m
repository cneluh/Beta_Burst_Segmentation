function dx=fdctIV(x)
% Fast DCT-IV 
% Refer: A Wavelet tour of Signal Processing, Mallat,1998
% Chapter VIII 

x=x(:);
N=length(x);
b=1:2:N;
c=N-b+1;
n=(0:1:N/2-1)';

g=(x(b)+i*x(c));
g=g.*exp(-i.*(n+0.25).*pi/N);
fg=fft(g);
G=fg(n+1).*exp(-i*pi*n/N);

dx=zeros(1,N);
dx(b)=real(G);
dx(c)=-imag(G);
dx=dx.*sqrt(2/N);

% Written by Nuri Firat INCE
% 07.06.2004
% firat@umn.edu