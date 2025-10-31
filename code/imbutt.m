function x=imbutt(x,cf,order)

[m,n]=size(x);%m=n
cf=cf*2*m;
px=zeros(2*m,2*n);
px(1:m,1:n)=x;
[I,J]=meshgrid(1:2*m,1:2*n);
s=(-1).^(I+J);
px=px.*s;
fx=fft2(px);
D=sqrt((I-(m+1)).^2+(J-(n+1)).^2);
H=1./(1+(D/cf).^(2*order));
fx=fx.*H;
px=real(ifft2(fx));
px=px.*s;
x=px(1:m,1:n);