function ncat_4D_ini(dummy)
%ncat_4D_ini.m  generate noise sinogram for reconstruction

%get noiseless sinogram(lesion): primary(sino_p), scatter(sino_s), lower window
%total(sino_L)
% sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
% for g=1:16
%     fid=fopen(['c:\simind\ncatlestew' num2str(g) 's.a02'],'rb');
%     bimS=fread(fid,'single');
%     fclose(fid);
%     bimS=reshape(bimS,[64,64,64]);
%     fid=fopen(['c:\simind\ncatlestew' num2str(g) 't.a02'],'rb');
%     bimT=fread(fid,'single');
%     fclose(fid);
%     bimT=reshape(bimT,[64,64,64]);
% 
%     fid=fopen(['c:\simind\ncatlestew' num2str(g) 't.a01'],'rb');
%     bimLT=fread(fid,'single');
%     fclose(fid);%use for real compensation
%     bimLT=reshape(bimLT,[64,64,64]);
%     %2 channels: 122.5-126kev, 126-154kev.
%     %scale and align for reconstruction
%     scal=2.5e5/sum(bimT(:))*2;%8 million
%     bimT=bimT*scal;
%     bimS=bimS*scal;
%     bimLT=bimLT*scal;
%     for n=1:64
%         bimT(:,:,n)=rot90(bimT(:,:,n));
%         bimS(:,:,n)=rot90(bimS(:,:,n));
%         bimLT(:,:,n)=rot90(bimLT(:,:,n));
%     end
%     sino_p(:,:,:,g)=bimT-bimS;
%     sino_s(:,:,:,g)=bimS;
%     sino_L(:,:,:,g)=bimLT;
% end
%save ncatlestew_NL sino_p sino_s sino_L

%lesion 30 noise sinogram
load ncatlestew_NL sino_p sino_s sino_L
for n=1:30
    sinop=random('poiss',sino_p);
    sinos=random('poiss',sino_s);
    sinoT=sinop+sinos;
    sinoL=random('poiss',sino_L);
%     sino_L=reshape(sino_L,[64^2 64 16])*4;
%     tew_scat=buttlpf(sino_L,0.4*.634,3);
%     tew_scat=reshape(tew_scat,[64 64 64 16]);
    filename=['ncatles_n' num2str(n) '.mat'];
    save(filename,'sinoT','sinoL');
end

%sino for ideal noiseless reconstruction
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get noiseless sinogram(normal): primary(sino_p), scatter(sino_s), lower window
%total(sino_L)
sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
for g=1:16
    fid=fopen(['c:\simind\ncat16gtew' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncat16gtew' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncat16gtew' num2str(g) 't.a01'],'rb');
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
save ncat16gtew_NL sino_p sino_s sino_L

%normal 30 noise sinogram
load ncat16gtew_NL sino_p sino_s sino_L
%add some offset
random('poiss',sino_p(:,:,:,1:8));
for n=1:30
    sinop=random('poiss',sino_p);
    sinos=random('poiss',sino_s);
    sinoT=sinop+sinos;
    sinoL=random('poiss',sino_L);
%     sino_L=reshape(sino_L,[64^2 64 16])*4;
%     tew_scat=buttlpf(sino_L,0.4*.634,3);
%     tew_scat=reshape(tew_scat,[64 64 64 16]);
    filename=['ncat16g_n' num2str(n) '.mat'];
    save(filename,'sinoT','sinoL');
end

%sino for ideal noiseless reconstruction
sino=zeros(64,64,64,16);
for g=1:16    
    fid=fopen(['c:\simind\ncat16gtew' num2str(g) 'a.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
    end
    sino(:,:,:,g)=bimT;
end
sino=sino*8e6/sum(sino(:));
save ncat16gtew_air sino
%see ncat16gtew_ideal_re.m

%get noiseless sinogram(lesion2): primary(sino_p), scatter(sino_s), lower window
%total(sino_L)
sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
for g=1:16
    fid=fopen(['c:\simind\ncatles2tew' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncatles2tew' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncatles2tew' num2str(g) 't.a01'],'rb');
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
save ncatles2tew_NL sino_p sino_s sino_L

%lesion 30 noise sinogram
load ncatles2tew_NL sino_p sino_s sino_L
for n=1:30
    sinop=random('poiss',sino_p);
    sinos=random('poiss',sino_s);
    sinoT=sinop+sinos;
    sinoL=random('poiss',sino_L);
%     sino_L=reshape(sino_L,[64^2 64 16])*4;
%     tew_scat=buttlpf(sino_L,0.4*.634,3);
%     tew_scat=reshape(tew_scat,[64 64 64 16]);
    filename=['ncatles2_n' num2str(n) '.mat'];
    save(filename,'sinoT','sinoL');
end

%sino for ideal noiseless reconstruction
sino=zeros(64,64,64,16);
for g=1:16    
    fid=fopen(['c:\simind\ncatles2tew' num2str(g) 'a.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
    end
    sino(:,:,:,g)=bimT;
end
sino=sino*8e6/sum(sino(:));
save ncatles2tew_air sino


%get noiseless sinogram(lesion LAD king): primary(sino_p), scatter(sino_s), lower window
%total(sino_L)
sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
for g=1:16
    fid=fopen(['c:\simind\ncatlesktew' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncatlesktew' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncatlesktew' num2str(g) 't.a01'],'rb');
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
save ncatlesktew_NL sino_p sino_s sino_L

%lesion 30 noise sinogram
load ncatlesktew_NL sino_p sino_s sino_L
for n=1:30
    sinop=random('poiss',sino_p);
    sinos=random('poiss',sino_s);
    sinoT=sinop+sinos;
    sinoL=random('poiss',sino_L);
%     sino_L=reshape(sino_L,[64^2 64 16])*4;
%     tew_scat=buttlpf(sino_L,0.4*.634,3);
%     tew_scat=reshape(tew_scat,[64 64 64 16]);
    filename=['ncatlesk_n' num2str(n) '.mat'];
    save(filename,'sinoT','sinoL');
end

%sino for ideal noiseless reconstruction
sino=zeros(64,64,64,16);
for g=1:16    
    fid=fopen(['c:\simind\ncatlesktew' num2str(g) 'a.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
    end
    scal=2.5e5/sum(bimT(:))*2;
    sino(:,:,:,g)=bimT*scal;
end
%sino=sino*8e6/sum(sino(:));
save ncatlesktew_air sino

%% 256 down to 64: Septal wall low contrast lesion
sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
for g=1:16
    fid=fopen(['c:\simind\ncatd64l' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncatd64l' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncatd64l' num2str(g) 't.a01'],'rb');
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
save ncatd64l sino_p sino_s sino_L
%add noise
load ncatd64l sino_p sino_s sino_L
for n=1:30
    sinop=random('poiss',sino_p);
    sinos=random('poiss',sino_s);
    sinoT=sinop+sinos;
    sinoL=random('poiss',sino_L);
%     sino_L=reshape(sino_L,[64^2 64 16])*4;
%     tew_scat=buttlpf(sino_L,0.4*.634,3);
%     tew_scat=reshape(tew_scat,[64 64 64 16]);
    filename=['ncatd64l_n' num2str(n) '.mat'];
    save(filename,'sinoT','sinoL');
end
%air
sino=zeros(64,64,64,16);
for g=1:16    
    fid=fopen(['c:\simind\ncatd64l' num2str(g) 'a.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
    end
    sino(:,:,:,g)=bimT;
end
sino=sino*8e6/sum(sino(:));
save ncatd64l_air sino
%% 256 down to 64: Normal
sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
for g=1:16
    fid=fopen(['c:\simind\ncatd64' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncatd64' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncatd64' num2str(g) 't.a01'],'rb');
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
save ncatd64 sino_p sino_s sino_L
%add noise
load ncatd64 sino_p sino_s sino_L
for n=1:30
    sinop=random('poiss',sino_p);
    sinos=random('poiss',sino_s);
    sinoT=sinop+sinos;
    sinoL=random('poiss',sino_L);
%     sino_L=reshape(sino_L,[64^2 64 16])*4;
%     tew_scat=buttlpf(sino_L,0.4*.634,3);
%     tew_scat=reshape(tew_scat,[64 64 64 16]);
    filename=['ncatd64_n' num2str(n) '.mat'];
    save(filename,'sinoT','sinoL');
end
%air
sino=zeros(64,64,64,16);
for g=1:16    
    fid=fopen(['c:\simind\ncatd64' num2str(g) 'a.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    for n=1:64
        bimT(:,:,n)=rot90(bimT(:,:,n));
    end
    sino(:,:,:,g)=bimT;
end
sino=sino*8e6/sum(sino(:));
save ncatd64_air sino