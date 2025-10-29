%get noiseless sinogram: primary(sino_p), scatter(sino_s), lower window
%total(sino_L)
sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
for g=1:16
    fid=fopen(['c:\simind\ncatlestew' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncatlestew' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncatlestew' num2str(g) 't.a01'],'rb');
    bimLT=fread(fid,'single');
    fclose(fid);%use for real compensation
    bimLT=reshape(bimLT,[64,64,64]);
    %2 channels: 122.5-126kev, 126-154kev.
    %scale and align for reconstruction
    scal=2.5e5/sum(bimT(:))*2;%8 million
    bimT=bimT*scal;
    bimS=bimS*scal;
    bimLT=bimLT*scal;
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
        bimS(:,:,n)=rot90(bimS(:,:,n));
        bimLT(:,:,n)=rot90(bimLT(:,:,n));
    end
    sino_p(:,:,:,g)=bimT-bimS;
    sino_s(:,:,:,g)=bimS;
    sino_L(:,:,:,g)=bimLT;
end
%save ncatlestew_NL sino_p sino_s sino_L
%noise data
load ncatlestew_NL sino_p sino_s sino_L
sino_p=random('poiss',sino_p);
sino_s=random('poiss',sino_s);
sino_L=random('poiss',sino_L);
save ncatlestew_n1 sino_p sino_s sino_L

%ideal noiseless reconstruction
sino=zeros(64,64,64,16);
for g=1:16    
    fid=fopen(['c:\simind\ncatlestew' num2str(g) 'a.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
    end
    sino(:,:,:,g)=bimT;
end
sino=sino*8e6/sum(sino(:));
save ncatlestew_air sino
%see ncatlestew_ideal_re.m

%TEW scatter
load ncatlestew_n1 sino_L
sino_L=reshape(sino_L,[64^2 64 16])*4;
tew_scat=buttlpf(sino_L,0.4*.634,3);
tew_scat=reshape(tew_scat,[64 64 64 16]);
save ncatlestew_esti_scat tew_scat

