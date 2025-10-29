sino_p=zeros(64,64,64,16);sino_s=sino_p;sino_L=sino_p;
g=1;
    fid=fopen(['c:\simind\ncat16gh' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncat16gh' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);

    fid=fopen(['c:\simind\ncat16gtew' num2str(g) 't.a02'],'rb');
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

sino=sino*2.5e5/sum(sino(:));
n_sino=random('poiss',sino);

g=1;
fid=fopen(['c:\simind\ncat16ghh' num2str(g) 's.a02'],'rb');
    bimS=fread(fid,'single');
    fclose(fid);
    bimS=reshape(bimS,[64,64,64]);
    fid=fopen(['c:\simind\ncat16ghh' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    
    %07/09/2007
    g=1;
   filename=['c:\simind\ncat16g_actOrg_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=temp*10;
    filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
    
    fid=fopen(['c:\simind\ncat16gp' num2str(g) 't.a02'],'rb');
    bimT=fread(fid,'single');
    fclose(fid);
    bimT=reshape(bimT,[64,64,64]);
    
    %using downsampled 256^3 phantom!!!
    g=2;
    filename=['C:\simind\ncatd64_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    dsp(temp)