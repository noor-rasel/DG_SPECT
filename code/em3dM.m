function x=em3dM(d,x,iter,blur,attn,scat)%,li,sbeta,gbeta,motion)
%EM with attenuation and scatter
%Mingwu Jin, Oct. 05, 2006

if nargin<4
    blur=1;
    attn=1;
    scat=0;
end%default SPECT with attenuation and without scatter
[nx,ny,nz,T,G]=size(x);
[ns,nr,na]=size(d);%nr=nx=ny; ns=nz!
if attn==1
    load weight64_mn
    load weight64_attn1.mat
    wp_wgt=wp_attnwgt;clear wp_attnwgt;
else
    load weight64_mn %no attenuation
end
load gbk64

% fp=fopen('roi_p64d64s64np64_3D','r');
% roi=fread(fp,'uint16');
% fclose(fp);
% roi=repmat(roi,1,64);
% ind=find(roi>0);clear roi

rsc=zeros(nx,ny,nz);
for n=1:64
    if n<33
        rsc=rsc+back3d_sa(ones(nz,nx),n,wp_vray((n-1)*64+1:n*64),...
            wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
    else
        m=n-32;
        if attn==1%attenuation weights
            rsc=rsc+back3d_sa(ones(nz,nx),n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
        else
            rsc=rsc+back3d_sa(ones(nz,nx),n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
        end
    end
end
ind=find(rsc>0);

for k=1:iter
    y=zeros(ns,nr,na);
    %projection
    for n=1:64
        if n<33
            y(:,:,n)=proj3d_sa(x,n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
        else
            m=n-32;
            if attn==1%attenuation weights
                y(:,:,n)=proj3d_sa(x,n,wp_vray((m-1)*64+1:m*64),...
                    wp_ipxl((m-1)*64+1:m*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
            else
                y(:,:,n)=proj3d_sa(x,n,wp_vray((m-1)*64+1:m*64),...
                    wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
            end
        end
    end
    %scatter correction
    y=y+scat;
    pro_nz=find(y>0);
    y(pro_nz)=d(pro_nz)./y(pro_nz);
    %backprojection
    t=zeros(nx,ny,nz);
    for n=1:64
        if n<33
            t=t+back3d_sa(y(:,:,n),n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
        else
            m=n-32;
            if attn==1%attenuation weights
                t=t+back3d_sa(y(:,:,n),n,wp_vray((m-1)*64+1:m*64),...
                    wp_ipxl((m-1)*64+1:m*64),wp_wgt((n-1)*64+1:n*64),blur,gb_temp);
            else
                t=t+back3d_sa(y(:,:,n),n,wp_vray((m-1)*64+1:m*64),...
                    wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),blur,gb_temp);
            end
        end
    end
    t(ind)=t(ind)./rsc(ind);
    x=x.*t;
end

