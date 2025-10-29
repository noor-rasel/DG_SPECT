%spatial domain Gaussian blurring
fp=fopen('.rec_p128d128s64np128_3D','r');
x=fread(fp,'single');
fclose(fp);
pxl_wd=x(1);
plane_wd=x(4);
ray_wd=x(2);
ror=x(3);
wnplane=128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters can be changed based on real imaging system! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slope=0.0137772;
%slope_z=slope;
intercept=0.134814;
%intercept_z=intercept;

%/* blur(sigma)=SLOPE*distance+INTERCEPT cm */
i=0:wnplane;
ftmp=(ror-wnplane/2*plane_wd)+plane_wd*i+plane_wd/2;
sigma=(slope*ftmp+intercept);%/ray_wd; no need for nomilizing to pixel index
%because the real coordinates are used below!!! fixed Mar. 31, 2006
%sigma_z=sigma;

%Gaussian blob
% d=0:10;d=repmat(d',1,65);d=d*pxl_wd;
% wd=exp(-d.^2./(2*repmat(sigma,11,1).^2));plot(wd')
% for n=1:64
%     t_ind=find(wd(:,n)<.001);
%     wdd(n)=t_ind(1);
% end
%Gaussian spatial template with right coordinates
[gb_x,gb_y]=meshgrid(-63*pxl_wd:pxl_wd:64*pxl_wd,-63*pxl_wd:pxl_wd:64*pxl_wd);
gb_temp=cell(128,1);
for n=1:128
    x=exp(-(gb_x.^2+gb_y.^2)/2/sigma(n)^2)/sqrt(2*pi)/sigma(n);
    mask=x>.1e-4;x=x.*mask;x=x/sum(x(:));
    x_ind=find(x(64,:)>0);
    x_min=x_ind(1);x_max=x_ind(end);
    gb_temp{n}=x(x_min:x_max,x_min:x_max);
end
save gbk128 gb_temp
% x=zeros(64,64,64);tic
% for n=1:64
% %     x(:,:,n)=filter2(gb_temp{n},temp(:,:,n));
%     x(:,:,n)=conv2(temp(:,:,n),gb_temp{n},'same');
% end
% toc%filter2 0.078; conv2 0.062
%distance>4 very small!
%plane 0~17: 3*3
%plane 18~37: 5*5
%plane 38~56: 7*7
%plane 57~63: 9*9
% load test_temp
%First way to do blurring(wrong coordinates)
% gb_temp1=zeros(3,3,18);
% gb_temp2=zeros(5,5,20);
% gb_temp3=zeros(7,7,19);
% gb_temp4=zeros(9,9,7);
% [gb_x,gb_y]=meshgrid(1:3,1:3);
% for n=1:18    
%     gb_temp1(:,:,n)=exp(-((gb_x-2).^2+(gb_y-2).^2)/2/sigma(n)^2);
% end
% [gb_x,gb_y]=meshgrid(1:5,1:5);
% for n=19:38    
%     gb_temp2(:,:,n-18)=exp(-((gb_x-3).^2+(gb_y-3).^2)/2/sigma(n)^2);
% end
% [gb_x,gb_y]=meshgrid(1:7,1:7);
% for n=39:57    
%     gb_temp3(:,:,n-38)=exp(-((gb_x-4).^2+(gb_y-4).^2)/2/sigma(n)^2);
% end
% [gb_x,gb_y]=meshgrid(1:9,1:9);
% for n=58:64    
%     gb_temp4(:,:,n-57)=exp(-((gb_x-5).^2+(gb_y-5).^2)/2/sigma(n)^2);
% end
% x=zeros(64,64,64);tic
% for n=1:18
%     x(:,:,n)=filter2(gb_temp1(:,:,n),temp(:,:,n));
% end
% for n=19:38
%     x(:,:,n)=filter2(gb_temp2(:,:,n-18),temp(:,:,n));
% end
% for n=39:57
%     x(:,:,n)=filter2(gb_temp3(:,:,n-38),temp(:,:,n));
% end
% for n=58:64
%     x(:,:,n)=filter2(gb_temp4(:,:,n-57),temp(:,:,n));
% end
% toc%0.172~0.079sec
%Another way(wrong coordinates)
% gb_temp=zeros(9,9,64);
% [gb_x,gb_y]=meshgrid(1:9,1:9);
% for n=1:64    
%      gb_temp(:,:,n)=exp(-((gb_x-5).^2+(gb_y-5).^2)/2/sigma(n)^2);
% end
% gb_temp(gb_temp<0.001)=0;tic
% for n=1:64
%     x(:,:,n)=filter2(gb_temp(:,:,n),temp(:,:,n));
% end
% toc%0.14~0.094
%proj3d.c way(wrong coordinates): slowest!!!
% tic
% for n=1:64
%     x(:,:,n)=fft2(temp(:,:,n));
%     x(:,:,n)=x(:,:,n).*gauss(1:64,1:64,n);
%     x(:,:,n)=real(ifft2(x(:,:,n)));
% end
% toc%0.7sec

% fp=fopen('roi_p64d64s64np64_3D','r');
% roi=fread(fp,'uint16');
% fclose(fp);