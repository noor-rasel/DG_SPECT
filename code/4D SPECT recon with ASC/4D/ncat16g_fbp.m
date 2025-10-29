%16 gates NCAT phantom, sinogram and FBP
load weight64_mn
load gbk64
load roi
ncat_phantom=zeros(64,64,64,16);
sino=zeros(64,64,64,16);
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'r');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    ncat_phantom(:,:,:,g)=temp;
    for n=1:64
        if n<33
            sino(:,:,n,g)=proj3d_sa(temp,n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
        else
            m=n-32;
            sino(:,:,n,g)=proj3d_sa(temp,n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),1,gb_temp);
        end
    end
end
ncat_phantom=ncat_phantom*4e6/sum(sino(:));%4e6/sum(ncat_phantom(:))=0.2192
sino=sino*4e6/sum(sino(:));%4e6/sum(sino(:))=0.2188
save ncat4D16g_input ncat_phantom sino
nsino=random('poiss',sino);
Im_rawfbp=fbp_3dM(nsino);
save ncat16g_rawfbp Im_rawfbp
% temp=clinicfilt3d_yyjin(temp,2.4,.4);
% temp=imbutt3d(temp,.2,2.4);
% temp(temp<0)=0;
% temp=temp*sum(nsino(:))/sum(temp(:));
load ncat4D16g_input sino
for n=1:16
    count_gate(n)=sum(sum(sum(sino(:,:,:,n))));
end
load ncat16g_rawfbp
temp=Im_rawfbp(:,:,:,1);ncatrawfbp_g1=temp;
temp=imbutt3d(temp,.2,2.4);
temp(temp<0)=0;
temp=temp*2.5e5/sum(temp(:));
ncatfbp_g1=temp;
temp=Im_rawfbp(:,:,:,2);ncatrawfbp_g2=temp;
temp=imbutt3d(temp,.2,2.4);
temp(temp<0)=0;
temp=temp*2.5e5/sum(temp(:));
ncatfbp_g2=temp;