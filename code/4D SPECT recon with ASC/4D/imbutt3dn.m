function x=imbutt3dn(x,cf,order)
%3D Butterworth no padding. Aliasing happens
%x: 3D image. cf: cut-off freq. (0~1).
%Mingwu Jin, MIRC IIT
%May 26, 2006

[m,n,k]=size(x);%m=n
step=[2/m 2/n 2/k];
[I,J,K]=meshgrid(-1:step(1):1-step(1),-1:step(2):1-step(2),-1:step(3):1-step(3));
D=sqrt(I.^2+J.^2+K.^2);
H=1./(1+(D/cf).^(2*order));
H=circshift(H,[m/2,n/2,k/2]);
fx=fftn(x);
fx=fx.*H;
x=real(ifftn(fx));