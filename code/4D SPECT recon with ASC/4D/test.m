load weight_planes
wp1=wp(:,:,1);
wp1=sparse(wp1);
load test_temp
temp=reshape(temp,64^2,64);
v_slice=temp.*repmat(wp1(:,1),1,64);


load weight64_mn
load gbk64

temp=permute(temp,[2 1 3]);
dsp(temp)
%validate projection
k=33;
if k>32
    ang=k-32;
else
    ang=k;
end
tic,y=proj3d_sa(temp,k,wp_vray((ang-1)*64+1:ang*64),...
    wp_ipxl((ang-1)*64+1:ang*64),wp_wgt((ang-1)*64+1:ang*64),0,gb_temp);toc
if 17-k+1<=0
    kc=17-k+1+64;
else
    kc=17-k+1;
end
tic,Sin=proj3d(temp,1,1,'-p',64,'-d',64,'-s',64,...
    '-firstAngle',kc,'-lastAngle',kc,'-ror',28.5,'-np',64);toc
dsp(y),dsp(Sin(:,:,kc))
%proj3d_sa, Angle 1 <=> proj3d, Angle 17! Same for PET projection!
%                 2 <==>              16
%                 ...
%
%relationship: ang <==> 17-ang+1 (if <=0, +64)


%validate backprojection
tic,xm=back3d_sa(y,1,wp_vray(1:64),wp_ipxl(1:64),wp_wgt(1:64),0,gb_temp);toc
tic,xc=back3d(Sin,1,1,'-p',64,'-d',64,'-s',64,'-firstAngle',17,'-lastAngle',17,'-ror',28.5,'-np',64);toc
%same!!!

nx=64;
rsc=zeros(64,64,64);tic
for n=1:64
    if n<33
        rsc=rsc+back3d_sa(ones(nx),n,wp_vray((n-1)*64+1:n*64),...
            wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),0,gb_temp);
    else
        m=n-32;
        rsc=rsc+back3d_sa(ones(nx),n,wp_vray((m-1)*64+1:m*64),...
            wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),0,gb_temp);
    end
end;toc
rscc=ones(64,64,64);
rscc=back3d(rscc,1,1,'-p',64,'-d',64,'-s',64,'-ror',28.5,'-np',64);
rscc=rscc*NC;

dsp(sum(permute(temp,[3,1,2])))
dsp(sum(permute(temp,[3,1,2]),3))%<=> dsp(squeeze(sum(temp,2))')

temp=zeros(64,64,64);
temp(12:20,20,25:35)=1;
%first projection plane is perpendicular with 'column axis' at '0'

tic,y=zeros(64,64,64);
for n=1:64
    if n<33
        y(:,:,n)=proj3d_sa(temp,n,wp_vray((n-1)*64+1:n*64),...
            wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
    else
        m=n-32;
        y(:,:,n)=proj3d_sa(temp,n,wp_vray((m-1)*64+1:m*64),...
            wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),1,gb_temp);
    end
end,toc

x=back3d_sa(y(:,:,1),1,wp_vray(1:64),wp_ipxl(1:64),wp_wgt(1:64),1,gb_temp);
%speed: C-routine 26.625sec; MATLAB 13.265sec Local machine

%test static reconstruction
load test_temp
load weight64_mn
load gbk64
%C
Sin=proj3d(temp,1,1,'-p',64,'-d',64,'-s',64,'-ror',28.5,'-np',64);
Sin=Sin*sum(temp(:))/sum(Sin(:));
xc=em3d(Sin,ones(64,64,64),1);
%MATLAB
y=zeros(64,64,64);
for n=1:64
    if n<33
        y(:,:,n)=proj3d_sa(temp,n,wp_vray((n-1)*64+1:n*64),...
            wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
    else
        m=n-32;
        y(:,:,n)=proj3d_sa(temp,n,wp_vray((m-1)*64+1:m*64),...
            wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),1,gb_temp);
    end
end
blur=1;
tic,xm=em3dM(y,ones(64,64,64),20,blur);toc
tic,xr=em3dM(y(25:48,:,:),ones(64,64,24),20,blur);toc

%NCAT phantom
load ncat8gate Sin
nsin=random('poiss',Sin);
y=zeros(64,64,64,8);
for g=1:8
    y(:,:,:,g)=fbp_3d(nsin(:,:,:,g));
end
save fbp_ncat8g y
y=reshape(y,64^2,64,8);
sy=buttlpf(y,0.4,5);
sy=reshape(sy,64,64,64,8);
sy=flipdim(sy,2);
sy=flipdim(sy,3);%for right views in AMIDE!!!
fid=fopen('fbpbutt_ncat8g.bin','wb');
fwrite(fid,sy,'float');
fclose(fid);

x=flipdim(x,2);
x=flipdim(x,3);
fid=fopen('ncat8g.bin','wb');
fwrite(fid,x,'float');
fclose(fid);