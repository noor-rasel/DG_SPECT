function x=imbutt3d(x,cf,order)
%3D Butterworth
%x: 3D image. cf: cut-off freq. (0~1).
%Mingwu Jin, MIRC IIT
%May 26, 2006

[m,n,k]=size(x);%m=n
cf=cf*2*m;
px=zeros(2*m,2*n,2*k);
px(1:m,1:n,1:k)=x;
[I,J,K]=meshgrid(1:2*m,1:2*n,1:2*k);
% s=(-1).^(I+J+K);
% px=px.*s;
fx=fftn(px);
D=sqrt((I-(m+1)).^2+(J-(n+1)).^2+(K-(k+1)).^2);
H=1./(1+(D/cf).^(2*order));
H=circshift(H,[m,n,k]);
fx=fx.*H;
px=real(ifftn(fx));
% px=px.*s;
x=px(1:m,1:n,1:k);