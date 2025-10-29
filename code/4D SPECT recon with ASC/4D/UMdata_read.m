dir_name='D:\imagereconstruction\4D\real_dataUM\P53606\';
filename='R19_1_float';
fid=fopen([dir_name filename],'rb');
temp=fread(fid,128^2*120,'single');
fclose(fid);
temp=reshape(temp,[128 128 120]);
% temp=permute(temp,[2 3 1]);
dsp(temp)
for n=1:120
%     dsp(temp(:,:,n));title(num2str(n)),pause;close
max_t(n)=max(max(temp(:,:,n)));
end
ind=find(max_t==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sinogram
dir_name='.\real_dataUM\P53606\R19_';
filename='_float';
for g=1:8
    fid=fopen([dir_name num2str(g) filename],'rb');
    temp=fread(fid,128^2*120,'single');
    fclose(fid);
    temp=reshape(temp,[128 128 120]);
    temp=permute(temp,[2 1 3]);
    temp=temp(:,:,42:109);
    sino(:,:,:,g)=temp(23:86,:,:);
end
save real_16g_sino sino
%attn
map_name='map_ostr.R35';
fid=fopen([dir_name map_name],'rb');
attn=fread(fid,128^3,'single');
fclose(fid);
attn=reshape(attn,[128 128 128]);
dsp(attn)
for n=1:120
    dsp(attn(:,:,n));pause;close
end
%scatter
dir_name='.\real_dataUM\P53606\R20_';
filename='_float';
scal=0.5*147.5*.14/(133.2*.05);
scat=zeros(64,128,68,8);
cf=0.2*.467;
for g=1:8
    fid=fopen([dir_name num2str(g) filename],'rb');
    temp=fread(fid,128^2*120,'single');
    fclose(fid);
    temp=reshape(temp,[128 128 120]);
    temp=permute(temp,[2 1 3]);
    temp=temp(:,:,42:109)*scal;
    tot_c=sum(temp(:));
    for n=1:68
        temp(:,:,n)=imbutt(temp(:,:,n),cf,3);
    end
    temp(temp<0)=0;
    temp=temp*tot_c/sum(temp(:));
    scat(:,:,:,g)=temp(23:86,:,:);
end
save real_16g_scat scat

%%%%%%%%%%%%%%%%%%%%
%Mingwu's ST121
load fbp_realdata42 rim

imir=rim(:,:,23:86);
imir=clinicfilt3d_yyjin(imir,2.4,.4);

imir=reshape(imir,[128^2*64 8]);
w_st121=[0.5 0.25 zeros(1,5) 0.25];G=8;
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end
imir121=imir*wm_st121;
imir121=reshape(imir121,[128 128 16 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inverse radon transform
temp=permute(temp,[2 1 3]);
temp=temp(:,:,42:109);
for n=1:68,imagesc(temp(32:32+63,32:32+63,n),[0 max(temp(:))]);pause;end

%65~78 heart
temp=permute(temp,[2 1 3]);
temp=temp(:,:,42:109);

n=23:105;immax=0.6471;
for m=1:30
sino=squeeze(temp(n(m),:,1:60));

imir=iradon(sino,123:3:300,'linear','Ram-Lak',1,128);
imir=imbutt(imir,0.2,5);
imir(imir<0)=0;
% tmax=max(imir(:));
% if immax<tmax
%     immax=tmax;%0.6471
% end
imagesc(imir,[0 immax]);title(num2str(n(m)));
mov(m)=getframe;pause(0.1)
end

imir=iradon(sino,0:3:177,'linear','Ram-Lak',1,128);

%8 gate slice 45
%immax=0.7725;slice=45;
slice=40:55;
for n=1:16
    for g=1:8
        dir_name='D:\imagereconstruction\4D\real_dataUM\P53606\';
        filename=['R19_' num2str(g) '_float'];
        fid=fopen([dir_name filename],'rb');
        temp=fread(fid,128^2*120,'single');
        fclose(fid);
        temp=reshape(temp,[128 128 120]);
        temp=permute(temp,[2 1 3]);
        temp=temp(:,:,42:109);
        sino=squeeze(temp(slice(n),:,1:60));
        imir(:,:,n,g)=iradon(sino,45:3:222,'linear','Ram-Lak',1,128);%123:3:300
        imir(:,:,n,g)=flipud(rot90(imir121(:,:,n,g),-1));
        imir(:,:,n,g)=imbutt(imir(:,:,n,g),0.2,2.4);
        %imir(imir<0)=0;
        %
        % imagesc(imir,[0 immax]);axis image,clinicalcolor%title(num2str(n(m)));
        % mov(g)=getframe;pause(0.1)
    end
end
imir=reshape(imir,[128^2*16 8]);
w_st121=[0.5 0.25 zeros(1,5) 0.25];G=8;
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end
imir121=imir*wm_st121;
imir121=reshape(imir121,[128 128 16 8]);
%save real_data_s40_51 imir121
save real_data_s40_55_cf02 imir121

im_show=imir121(:,:,4:12,:);
im_comb=zeros(128*3,128*3,8);
im_comb=zeros(90*3,90*3,8);
for n=1:3
    for m=1:3
        im_comb((n-1)*90+1:n*90,(m-1)*90+1:m*90,:)=squeeze(im_show(21:110,21:110,(m-1)*3+n,:));
    end
end
immax=max(im_comb(:));
for n=1:8
    imagesc(im_comb(:,:,n),[0 immax]);clinicalcolor,axis image,axis off
    mov(n)=getframe;pause(0.1);
end
movie2avi(mov,'slice43_51');

for n=1:8
    imagesc(im_comb(:,:,n),[0 immax]);clinicalcolor,axis image,axis off
    mov(n)=getframe;pause(0.1);
end
movie2avi(mov,'slice43_51mag');


clear
slice=23:10:103;
for n=1:9
    for g=1:8
        dir_name='D:\imagereconstruction\4D\real_dataUM\P53606\';
        filename=['R19_' num2str(g) '_float'];
        fid=fopen([dir_name filename],'rb');
        temp=fread(fid,128^2*120,'single');
        fclose(fid);
        temp=reshape(temp,[128 128 120]);
        temp=permute(temp,[2 1 3]);
        temp=temp(:,:,42:109);
        sino=squeeze(temp(slice(n),:,1:60));
        imir(:,:,n,g)=iradon(sino,45:3:222,'linear','Ram-Lak',1,128);%123:3:300
        imir(:,:,n,g)=imbutt(imir(:,:,n,g),0.1,2.4);
        %imir(imir<0)=0;
        %
        % imagesc(imir,[0 immax]);axis image,clinicalcolor%title(num2str(n(m)));
        % mov(g)=getframe;pause(0.1)
    end
end
imir=reshape(imir,[128^2*9 8]);
w_st121=[0.5 0.25 zeros(1,5) 0.25];G=8;
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end
imir121=imir*wm_st121;
imir121=reshape(imir121,[128 128 9 8]);
save real_data_s23_103 imir121


imir(imir<0)=0;
immax=max(imir(:));
for g=1:8
    imagesc(cat(2,imir(:,:,g),imir121(:,:,g)),[0 immax]);axis image,clinicalcolor%title(num2str(n(m)));
mov(g)=getframe;pause(0.1)
end

imir(imir<0)=0;
immax=max(imir(:));
for g=1:8
    imagesc(imir121(:,:,g),[0 immax]);axis image,clinicalcolor%title(num2str(n(m)));
mov(g)=getframe;pause(0.1)
end


figure,movie(mov,10)
% movie2avi(mov,'realdata_fbp_slice45')
movie2avi(mov,'realdata_fbp_slice47')

%%%%%%%%%%%%%%%%%%
weigen3d('-p',128,'-d',128,'-s',120,'-np',128,'-rpx',0.467,'-apx',0.467,'-ror',32.31)
%,'-roi',53.2

fid=fopen('roi_p128d128s120np128_3D','rb');
temp=fread(fid,'single');
fclose(fid);
temp=reshape(temp,[128 64]);

