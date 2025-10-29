
load weight64_mn
load gbk64
load weight64_attn1.mat
wp_wgt=wp_attnwgt;clear wp_attnwgt;
load ncat4D16g_input ncat_phantom
%sinogram with attenuation
tic
sino=zeros(64,64,64,16);
for g=1:16
    temp=ncat_phantom(:,:,:,g);
    for n=1:64
        if n<33
            sino(:,:,n,g)=proj3d_sa(temp,n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
        else
            m=n-32;
            sino(:,:,n,g)=proj3d_sa(temp,n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((n-1)*64+1:n*64),1,gb_temp);
        end
    end
end
toc%238 seconds.
%sum(ncat_phantom(:))/4e6=0.9913
sino=sino*4e6/sum(sino(:));%4e6/sum(sino(:))=3.375
save ncat4D16gattn_sino sino
%FBP on attenuated sinogram
load ncat4D16gattn_sino sino
sino_attn=sino(:,:,:,1);
Im_fbp=fbp_3dM(sino_attn);
Im_fbp(Im_fbp<0)=0;Im_fbp=Im_fbp*sum(sino_attn(:))/sum(Im_fbp(:));
Im_fbp_attn=Im_fbp;
dsp(Im_fbp_attn(:,:,35:38));
%FBP on ideal sinogram
load ncat4D16g_input sino
sino=sino(:,:,:,1);
Im_fbp=fbp_3dM(sino);
Im_fbp(Im_fbp<0)=0;Im_fbp=Im_fbp*sum(sino(:))/sum(Im_fbp(:));
dsp(Im_fbp(:,:,35:38));

%EM
tic,Im_em=em3dM(sino_attn,ones(64,64,64),50);toc%1560 seconds
Im_em=Im_em*sum(sino_attn(:))/sum(Im_em(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scatter
%DPW method
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
bimS=bimLS+bimS;
%Ideal scatter estimation: noiseless
scal=2.5e5/sum(bimT(:));
bimT=bimT*scal;
bimS=bimS*scal;
for n=1:64
    bimT(:,:,n)=rot90(bimT(:,:,n));
    bimS(:,:,n)=rot90(bimS(:,:,n));
end

Im_fbp_s=fbp_3dM(bimT);%FBP:total counts
Im_fbp_s(Im_fbp_s<0)=0;Im_fbp_s=Im_fbp_s*2.5e5/sum(Im_fbp_s(:));
dsp(Im_fbp_s(:,:,35:38));
Im_fbp_p=fbp_3dM(bimT-bimS);%FBP:primary counts
Im_fbp_p(Im_fbp_p<0)=0;Im_fbp_p=Im_fbp_p*2.5e5*sum(bimT(:)-bimS(:))/sum(bimT(:))/sum(Im_fbp_p(:));
dsp(Im_fbp_p(:,:,35:38));
tic,Im_em=em3dM(bimT,ones(64,64,64),20,1,1,bimS);toc%EM 20: SC
Im_em=Im_em*2.5e5/sum(Im_em(:));
dsp(Im_em(:,:,35:38))
tic,Im_em50=em3dM(bimT,Im_em,30,1,1,bimS);toc%EM 50:AC + SC
Im_em50=Im_em50*2.5e5/sum(Im_em50(:));
dsp(Im_em50(:,:,35:38))
tic,Im_em_nsc=em3dM(bimT,ones(64,64,64),50);toc%EM 50:AC + No SC
Im_em_nsc=Im_em_nsc*2.5e5/sum(Im_em_nsc(:));
dsp(Im_em_nsc(:,:,35:38)
tic,Im_emn=em3dM(bimT,ones(64,64,64),50,1,0,0);toc%EM 50:no AC and SC
Im_emn=Im_emn*2.5e5/sum(Im_emn(:));
dsp(Im_emn(:,:,35:38))
%Ideal scatter estimation: noise
sino_p=bimT-bimS;
sino_p=random('poiss',sino_p);
sino_s=random('poiss',bimS);
Im_fbp_s=fbp_3dM(sino_p+sino_s);%FBP:total counts
Im_fbp_sf=clinicfilt3d_yyjin(Im_fbp_s,2.4,0.4);
Im_fbp_sf=Im_fbp_sf*sum(sino_p(:)+sino_s(:))/sum(Im_fbp_sf(:));
dsp(Im_fbp_s(:,:,35:38));dsp(Im_fbp_sf(:,:,35:38));
Im_fbp_p=fbp_3dM(sino_p);%FBP:primary counts
Im_fbp_pf=clinicfilt3d_yyjin(Im_fbp_p,2.4,0.4);
Im_fbp_pf=Im_fbp_pf*sum(sino_p(:))/sum(Im_fbp_pf(:));
dsp(Im_fbp_p(:,:,35:38));dsp(Im_fbp_pf(:,:,35:38));
tic,Im_em_sc=em3dM(sino_p+sino_s,ones(64,64,64),50,1,1,sino_s);toc%EM 50: AC + SC
% Im_em_sc=Im_em_sc*2.5e5/sum(Im_em_sc(:));
% dsp(Im_em_sc(:,:,35:38));
Im_em_scf=clinicfilt3d_yyjin(Im_em_sc,2.4,0.4);
Im_em_scf=Im_em_scf*sum(sino_p(:))/sum(Im_em_scf(:));
dsp(Im_em_scf(:,:,35:38));
tic,Im_em_nsc=em3dM(sino_p+sino_s,ones(64,64,64),50);toc%EM 50: AC + No SC
Im_em_nscf=clinicfilt3d_yyjin(Im_em_nsc,2.4,0.4);
Im_em_nscf=Im_em_nscf*sum(sino_p(:)+sino_s(:))/sum(Im_em_nscf(:));
dsp(Im_em_nscf(:,:,35:38));

tic,Im_em=em3dM(sino_p+sino_s,ones(64,64,64),50,1,0,0);toc%EM 50: AC + SC
Im_emf=clinicfilt3d_yyjin(Im_em,2.4,0.3);
Im_emf=Im_emf*sum(sino_p(:)+sino_s(:))/sum(Im_emf(:));
dsp(Im_emf(:,:,35:38));