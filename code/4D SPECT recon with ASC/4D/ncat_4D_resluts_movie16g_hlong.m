%movies
%short axis images
b_slice=8;
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,:);

M=30;N=28;G=16;
%slice 18 for 'short'; 28*20*30
%slice 16 for 'vlong'; 20*30*28
%slice 8 for 'hlong'; 30*28*20

cview='hlong';p_slice=8;
for n=1:G
    temp=Im_part(:,:,:,n);
    temp=cardiac3drot(temp,cview);
    id_short(:,:,n)=temp(:,:,p_slice);
end
immax=max(id_short(:));
Im_comb=zeros(M,N,G,b_slice);
Im_comb(:,:,:,1)=id_short;

w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end

noise_num=1;
filename={['D:\imagereconstruction\4D\data_les_fbp\Im_fbp_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_les_5st0\Im_ncatles_nAS_s3t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_les_5st0\Im_ncatles_A_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_les_5st0\Im_ncatles_AS_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_les_5st2\Im_ncatles_nAS_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_les_5st2\Im_ncatles_A_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_les_5st2\Im_ncatles_AS_s1t2_n' num2str(noise_num) '.mat']};
for k=1:7
    load(filename{k})
    if k>1
        Im_part=Im_maps(:,:,:,:);
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp,cview);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short;
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=Im_st121(:,:,:,:);
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp,cview);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short;
    end
end
Im_comb(Im_comb>immax)=immax;

Im_combs=zeros((M+2)*2,(N+2)*4,G);
im_index=[1 2 3 6 4 7 5 8];
for g=1:G
    for n=1:2
        for k=1:4
            Im_combs((n-1)*(M+2)+1:n*M+(n-1)*2,(k-1)*(N+2)+1:k*N+(k-1)*2,g)=Im_comb(:,:,g,im_index((k-1)*2+n));
        end
    end
end
save(['data_4Dmovieles_' cview],'Im_combs');

%load data_4Dmovienor_hlong.mat
G=16
for g=1:16 
interp_im(:,:,g)=interp2(Im_combs(:,:,g),2);
end
interp_im=interp_im/max(interp_im(:));
for g=1:G
    figure('position',[100 100 800 500]),imagesc(interp_im(:,:,g),[0 1]),clinicalcolor,axis image%colormap(gray)
    axis off,im(g)=getframe;pause(.1)
    close
end
figure('position',[100 100 800 400]),movie(im,10)
movie2avi(im,['4Dmovieles_' cview])


%%%%%%%%%%%%%%%%%%%%%
%clinical
load myo_template
myo_ind=find(myo_temp>0);

p_slice=8;M=30;N=28;G=16;
b_slice=8;
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,:);
immax=max(Im_part(myo_ind));
Im_part(Im_part>immax)=immax;
Im_part=Im_part/immax;
for n=1:G
    temp=Im_part(:,:,:,n);
    temp=cardiac3drot(temp,'hlong');
    id_short(:,:,n)=temp(:,:,p_slice);
end
Im_comb=zeros(M,N,G,b_slice);%should normalize here to avoid intensity loss of thin heart wall!
Im_comb(:,:,:,1)=id_short/max(id_short(:));

w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end

noise_num=1;
filename={['D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_nAS_s3t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_A_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_AS_s4t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_nAS_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_A_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_AS_s1t2_n' num2str(noise_num) '.mat']};
for k=1:7
    load(filename{k})
    if k>1
        Im_part=Im_maps(:,:,:,:);
        immax=max(Im_part(myo_ind));
        Im_part(Im_part>immax)=immax;
        Im_part=Im_part/immax;
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp,'hlong');
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short/max(id_short(:));
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=Im_st121(:,:,:,:);
        immax=max(Im_part(myo_ind));
        Im_part(Im_part>immax)=immax;
        Im_part=Im_part/immax;
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp,'hlong');
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short/max(id_short(:));
    end
end

Im_combs=zeros((M+2)*2,(N+2)*4,G);
im_index=[1 2 3 6 4 7 5 8];
for g=1:G
    for n=1:2
        for k=1:4
            Im_combs((n-1)*(M+2)+1:n*M+(n-1)*2,(k-1)*(N+2)+1:k*N+(k-1)*2,g)=Im_comb(:,:,g,im_index((k-1)*2+n));
        end
    end
end
save data_4Dmovienor_clincal_hlong Im_combs

G=16
for g=1:16 
interp_im(:,:,g)=interp2(Im_combs(:,:,g),2);
end
interp_im=interp_im/max(interp_im(:));
for g=1:G
    figure('position',[100 100 800 500]),imagesc(interp_im(:,:,g),[0 1]),clinicalcolor,axis image%colormap(gray)
    axis off,im(g)=getframe;pause
    close
end
figure('position',[100 100 800 400]),movie(im,10)
movie2avi(im,'4Dmovienor_clinical_hlong')





%original
% clear org_short;p_slice=37;
% for g=1:16
%     filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
%     fid=fopen(filename,'rb');
%     temp=fread(fid,64^3,'single');    
%     fclose(fid);
%     temp=reshape(temp,[64 64 64]);    
% %     temp=temp(23:52,16:43,29:48);
%     temp=cardiac3drot(temp);
%     org_short(:,:,g)=temp(:,:,p_slice);
% end
% for g=1:G
%     figure('position',[100 100 400 400]),imagesc(interp2(org_short(:,:,g),1)),clinicalcolor,axis image%7:26
%     axis off,im(g)=getframe;
%     close
% end
% figure('position',[100 100 400 400]),movie(im,10)
% 
% movie2avi(im,'NCAT_org')
% movie2avi(im,'NCAT_org_trans')
% 
% clear org_short
% for g=1:16
%     filename=['D:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat128_act_' num2str(g) '.bin'];
%     fid=fopen(filename,'rb');
%     temp=fread(fid,128^2*64,'single');    
%     fclose(fid);
%     temp=reshape(temp,[128 128 64]);    
% %     temp=temp(23:52,16:43,29:48);
% % %     temp=cardiac3drot(temp);
%     org_short(:,:,g)=temp(:,:,63);
% end
% for g=1:G
%     figure('position',[100 100 400 400]),imagesc(interp2(org_short(:,:,g),1)),clinicalcolor,axis image%7:26
%     axis off,im(g)=getframe;
%     close
% end
% figure('position',[100 100 400 400]),movie(im,10)
% movie2avi(im,'NCAT_org_trans128')