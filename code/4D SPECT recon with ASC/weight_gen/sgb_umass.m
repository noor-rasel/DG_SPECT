fp=fopen('.rec_p64d64s64np64_3D','r');
x=fread(fp,'single');
fclose(fp);
pxl_wd=x(1);
plane_wd=x(4);
ray_wd=x(2);
ror=x(3);
wnplane=64;

%Bing's input
slope=0.0189;
slope_z=0.0195;
intercept=0.1266;
intercept_z=0.1271;

%/* blur(sigma)=SLOPE*distance+INTERCEPT cm */
i=0:wnplane;
ftmp=(ror-wnplane/2*plane_wd)+plane_wd*i+plane_wd/2;
sigma=(slope*ftmp+intercept);%/ray_wd; no need for nomilizing to pixel index
sigma_z=(slope_z*ftmp+intercept_z);

[gb_x,gb_y]=meshgrid(-31*pxl_wd:pxl_wd:32*pxl_wd,-31*pxl_wd:pxl_wd:32*pxl_wd);
gb_temp=cell(64,1);
for n=1:64
    x=exp(-(gb_x.^2/2/sigma(n)^2+gb_y.^2/2/sigma_z(n)^2))/sqrt(2*pi)/sigma(n)/sigma_z(n);
%     xf(:,:,n)=fft2(x);
    mask=x>.1e-4;x=x.*mask;x=x/sum(x(:));
    x_ind=find(x(32,:)>0);
    x_min=x_ind(1);x_max=x_ind(end);
    y_ind=find(x(:,32)>0);
    y_min=y_ind(1);y_max=y_ind(end);
    gb_temp{n}=x(x_min:x_max,y_min:y_max);%z along the column(y)
end
save gbk64_umass gb_temp

load gbk64
%speed test
x=zeros(64,64,64);temp=ones(64,64,64);
tic
for n=1:64
%     x(:,:,n)=filter2(gb_temp{n},temp(:,:,n));
    x(:,:,n)=conv2(temp(:,:,n),gb_temp{n},'same');
end
toc%filt
%0.125(gbk64_umass) vs 0.11(gbk64)