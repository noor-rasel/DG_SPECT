function M=motionmatrix3d(vx,vy,vz)
%Generate compensation matrix given motion field by 
%using bilinear interpolation.
%Mingwu Jin, MIRC IIT
%Apr 26, 2006
[m,n,k]=size(vx);
safe_buff=500;
M=spalloc(m*n*k,m*n*k,m*n*k+(nnz(vx(:))+safe_buff)*7);
[x,y,z]=meshgrid(1:n,1:m,1:k);%n is column(x). m is row(y)
x=x(:)-vx(:);
y=y(:)-vy(:);
z=z(:)-vz(:);
%decide outside region
xl1=find(x<1);x(xl1)=1;
xln=find(x>n);x(xln)=n;
yl1=find(y<1);y(yl1)=1;
ylm=find(y>m);y(ylm)=m;
zl1=find(z<1);z(zl1)=1;
zlk=find(z>k);z(zlk)=k;
xf=floor(x);xc=ceil(x);%eight apices:(xf,yf) (xc,yf) (xf,yc) (xc,yc)
yf=floor(y);yc=ceil(y);
zf=floor(z);zc=ceil(z);
dx=x-xf;
dy=y-yf;
dz=z-zf;
I=zeros(m*n*k,8);
I(:,1)=sub2ind([m n k],yf,xf,zf);
I(:,2)=sub2ind([m n k],yf,xc,zf);
I(:,3)=sub2ind([m n k],yc,xf,zf);
I(:,4)=sub2ind([m n k],yc,xc,zf);
I(:,5)=sub2ind([m n k],yf,xf,zc);
I(:,6)=sub2ind([m n k],yf,xc,zc);
I(:,7)=sub2ind([m n k],yc,xf,zc);
I(:,8)=sub2ind([m n k],yc,xc,zc);
D=zeros(m*n*k,8);
D(:,1)=(1-dx).*(1-dy).*(1-dz);
D(:,2)=dx.*(1-dy).*(1-dz);
D(:,3)=(1-dx).*dy.*(1-dz);
D(:,4)=dx.*dy.*(1-dz);
D(:,5)=(1-dx).*(1-dy).*dz;
D(:,6)=dx.*(1-dy).*dz;
D(:,7)=(1-dx).*dy.*dz;
D(:,8)=dx.*dy.*dz;
for j=1:m*n*k
        M(I(j,D(j,:)>0),j)=D(j,D(j,:)>0)';
end
%M=M';
% M([xl1;xln;yl1;ylm;zl1;zlk],:)=0;%reset outside region equal to zero
%using (x'*M)', otherwise M=M' ==> M*x; 