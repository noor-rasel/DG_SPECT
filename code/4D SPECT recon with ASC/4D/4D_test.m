%4D test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimal of the 3D projection matrix. To be used to find upper bound in
%BSREM II algorithm.
load weight64_mn
load gbk64
pj=zeros(64,64,64);
for n=1:64
    if n<33
        pj(:,:,n)=projMat_min(n,wp_vray((n-1)*64+1:n*64),...
            wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
    else
        m=n-32;
        pj(:,:,n)=projMat_min(n,wp_vray((m-1)*64+1:m*64),...
            wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64),1,gb_temp);
    end
end
save projMat_min pj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NCAT phantom and sinogram
ncat_phantom=zeros(64,64,64,8);
sino=zeros(64,64,64,8);
load weight64_mn
load gbk64
load roi
tic
for g=1:8
    filename=['D:\imagereconstruction\4D\ncat\wholencat8g_act_' num2str(g) '.bin'];
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
toc%115 seconds
ncat_phantom=ncat_phantom*4e6/sum(sino(:));
sino=sino*4e6/sum(sino(:));
save ncat4D8g_input ncat_phantom sino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test reconstruction on 1 gate with spatial prior first.
%noiseless
%FBP
tic,temp=fbp_3dM(sino(:,:,:,1));toc%12.5 seconds; 2 times slower than C.(fbp_3d) without blurring
templpf=clinicfilt3d_yyjin(temp,2.4,.4);
temp(temp<0)=0;
temp=temp*sum(sum(sum(sino(:,:,:,1))))/sum(temp(:));
%EM
tic,x=em3dM(sino(:,:,:,1),temp,10,1);toc%325 seconds
x=x*sum(sum(sum(sino(:,:,:,1))))/sum(x(:));
x10=em3dM(sino(:,:,:,1),ones(64,64,64),10,1);
x10=x10*sum(sum(sum(sino(:,:,:,1))))/sum(x10(:));
tic,x50=em3dM(sino(:,:,:,1),x10,40,1);toc%1214 seconds
x50=x50*sum(sum(sum(sino(:,:,:,1))))/sum(x50(:));
tic,x100=em3dM(sino,x50,50,1);toc
x100=x100*sum(sum(sum(sino(:,:,:,1))))/sum(x10(:));
%Noisy
%OS
load ncat4D16g_input sino
sino=sino(:,:,:,1);
nsino=random('poiss',sino);
temp=fbp_3dM(nsino);
templpf=imbutt3d(temp,.2,2.4);
templpf(templpf<0)=0;
templpf=templpf*sum(nsino(:))/sum(templpf(:));
im_fbp=templpf;
snr_fbp=10*log10(ncatg1(:)'*ncatg1(:)/sum((ncatg1(:)-im_fbp(:)).^2))
roi_snr=[25:49,16:40,34:38];%34:38
ncatg1roi=ncatg1(25:49,16:40,34:38);dsp(ncatg1roi)
Imroi=im_fbp(25:49,16:40,34:38);
snr_fbp_roi=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
load roi
tic,[ysit5u,cfsu]=mbsrem4d(nsino,repmat(roi,[1,1,64]),16,5,1,0.01,0,1);toc
%363.844 seconds.(calculation of the objective function included)
%283.5 seconds without prior
snr_su=10*log10(ncatg1(:)'*ncatg1(:)/sum((ncatg1(:)-ysit5u(:)).^2))
Imroi=ysit5u(25:49,16:40,34:38);
snr_roi=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))

templpf=imbutt3d(temp,.05,1);
templpf(templpf<0)=0;
templpf=templpf*sum(nsino(:))/sum(templpf(:));
tic,[ysit5,cfs]=mbsrem4d(nsino,templpf,16,5,1,0.01,0,1);toc
Imroi=ysit5(25:49,16:40,34:38);
snr_roi=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))

sbeta=0:0.001:0.01
Im_maps=zeros(64,64,64,11);cf_maps=zeros(11,5);
for n=1:11
    [Im_maps(:,:,:,n),cf_maps(n,:)]=mbsrem4d(nsino,repmat(roi,[1,1,64]),16,5,1,sbeta(n),0,1);%357 seconds
    Imroi=ysit5u(25:49,16:40,34:38,n);
    snr_maps_roi(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
end
save os_maps&fbp_g1 Im_maps cf_maps snr_maps_roi im_fbp snr_fbp_roi

sbeta=[0 0.001 0.005 0.01:0.01:0.1] 
Im_maps=zeros(64,64,64,13);cf_maps=zeros(13,20);snr_maps_roi=zeros(13,1);
for n=1:13
    [Im_maps(:,:,:,n),cf_maps{n}]=mbsrem4d(nsino,yit1roi,16,20,1,sbeta(n),0,1);%357 seconds
    Imroi=Im_maps(25:49,16:40,34:38,n);
    snr_maps_roi(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
end
save os_maps20newprior Im_maps cf_maps snr_maps_roi yit1roi nsino%cf_maps are not right & alpha=0.015/n^.33*2
%uniform start
sbeta=[0 0.001 0.005 0.01:0.01:0.1] 
Im_maps=zeros(64,64,64,13);cf_maps=cell(13);snr_maps_roi=zeros(13,1);
for n=1:13
    [Im_maps(:,:,:,n),cf_maps{n}]=mbsrem4d(nsino,repmat(roi,[1,1,64]),16,10,1,sbeta(n),0,1);%357 seconds
    Imroi=Im_maps(25:49,16:40,34:38,n);
    snr_maps_roi(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
end
save os_maps10newprior_uni Im_maps cf_maps snr_maps_roi
%fbp start
sbeta=[0 0.001 0.005 0.01:0.01:0.1] 
Im_maps=zeros(64,64,64,13);cf_maps=cell(13);snr_maps_roi=zeros(13,1);
for n=1:13
    [Im_maps(:,:,:,n),cf_maps{n}]=mbsrem4d(nsino,im_fbp,16,10,1,sbeta(n),0,1);%357 seconds
    Imroi=Im_maps(25:49,16:40,34:38,n);
    snr_maps_roi(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
end
save os_maps10newprior_fbp Im_maps cf_maps snr_maps_roi
%motion
%2 gate no motion
MM=kron([1;-1],speye(64^3));
save testMM1 MM
MM=kron([-1;1],speye(64^3));
save testMM2 MM
load ncat4D16g_input sino
sino=sino(:,:,:,1);
nsino=random('poiss',sino);
nsino2=random('poiss',sino);
nsino=cat(4,nsino,nsino2);
load roi
tic,[Im_mapt,cf_mapt]=mbsrem4d(nsino,repmat(roi,[1,1,64,2]),16,10,1,0,0.5,1);toc

filename='D:\imagereconstruction\4D\thomas_motion_ncat\ncat_mov16_vec_1_2.txt';
fid = fopen(filename);
fgetl(fid);
fgetl(fid);
C = textscan(fid,'%*s %*s %f %f %f %*s %f %f %f %*s %f %f %f');
fclose(fid);
Xref = C{1};
Yref = C{2};
Zref = C{3};
Xtar = C{4};
Ytar = C{5};
Ztar = C{6};
Vtx = C{7};
Vty = C{8};
Vtz = C{9};
[lgr,ind]=size(Xref);
[x,y,z]=meshgrid(1:64);
temp=interp3(x,y,z,ncat_phantom(:,:,:,1),Xref,Yref,Zref);
Xmax=ceil(max(Xref));Xmin=floor(min(Xref));
Ymax=ceil(max(Yref));Ymin=floor(min(Yref));
Zmax=ceil(max(Zref));Zmin=floor(min(Zref));
dsp(ncat_phantom(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%realistic projection: attenuation and scatter
fid=fopen('c:\simind\ncatdwp1s.a02','rb');
bimS=fread(fid,'single');
fclose(fid);
bimS=reshape(bimS,[64,64,64]);
fid=fopen('c:\simind\ncatdwp1t.a02','rb');
bimT=fread(fid,'single');
fclose(fid);
bimT=reshape(bimT,[64,64,64]);

fid=fopen('c:\simind\ncatdwp1t.a01','rb');
bimLT=fread(fid,'single');
fclose(fid);
bimLT=reshape(bimLT,[64,64,64]);
bimT=bimLT+bimT;
fid=fopen('c:\simind\ncatdwp1s.a01','rb');
bimLS=fread(fid,'single');
fclose(fid);
bimLS=reshape(bimLS,[64,64,64]);
bimS=bimLS+bimS;%2 channels: 126-133kev, 133-154kev.
%scale and align for reconstruction
scal=2.5e5/sum(bimT(:));
bimT=bimT*scal;
bimS=bimS*scal;
for n=1:64
    bimT(:,:,n)=rot90(bimT(:,:,n));
    bimS(:,:,n)=rot90(bimS(:,:,n));
end
%add noise for primary (=total-scatter) and scatter counts.
sino_p=bimT-bimS;
sino_p=random('poiss',sino_p);
sino_s=random('poiss',bimS);
save sino_attn_scat sino_p sino_s
load sino_attn_scat
%reconstruction: FBP
Im_fbp_s=fbp_3dM(sino_p+sino_s);%FBP:total counts
Im_fbp_sf=clinicfilt3d_yyjin(Im_fbp_s,2.4,0.4);
Im_fbp_sf=Im_fbp_sf*sum(sino_p(:)+sino_s(:))/sum(Im_fbp_sf(:));
Im_fbp_p=fbp_3dM(sino_p);%FBP:primary counts
Im_fbp_pf=clinicfilt3d_yyjin(Im_fbp_p,2.4,0.4);
Im_fbp_pf=Im_fbp_pf*sum(sino_p(:))/sum(Im_fbp_pf(:));
dsp(Im_fbp_sf(:,:,35:38));dsp(Im_fbp_pf(:,:,35:38));
%reconstruction: MAP
load roi
sub_num=16;it_num=10;OF_tag=1;
sbeta=0.05;gbeta=0;blur=1;
tic,Im_maps=mbsrem4d(sino_p+sino_s,repmat(roi,[1,1,64]),sub_num,it_num,OF_tag,sbeta,gbeta);toc
Im_maps=Im_maps*sum(sino_p(:)+sino_s(:))/sum(Im_maps(:));
dsp(Im_maps(:,:,35:38))
tic,Im_maps_ac=mbsrem4d(sino_p+sino_s,repmat(roi,[1,1,64]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,0);toc
Im_maps_ac=Im_maps_ac*sum(sino_p(:)+sino_s(:))/sum(Im_maps_ac(:));
dsp(Im_maps_ac(:,:,35:38))
tic,Im_maps_asc=mbsrem4d(sino_p+sino_s,repmat(roi,[1,1,64]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,sino_s);toc
Im_maps_asc=Im_maps_asc*sum(sino_p(:))/sum(Im_maps_asc(:));
dsp(Im_maps_asc(:,:,35:38))
%SNR
load ncat4D16g_input ncat_phantom
ncatg1=ncat_phantom(:,:,:,1);clear ncat_phantom
ncatg1roi=ncatg1(25:49,16:40,34:38);
%fbp: total
Imroi=Im_fbp_sf(25:49,16:40,34:38);
10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
%fbp: primary
Imroi=Im_fbp_pf(25:49,16:40,34:38)*sum(sino_p(:)+sino_s(:))/sum(sino_p(:));
10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))
%MAP: no AC and SC, sbeta=0.01
Imroi=Im_maps001(25:49,16:40,34:38)*sum(sino_p(:)+sino_s(:))/sum(Im_maps001(:));
10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2))