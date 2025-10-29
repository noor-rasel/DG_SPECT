%ncat 4D results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROI
%(1) myocardium
for g=1:16
    %     filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    d_name='D:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\';
    filename=[d_name 'myoncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    myo_temp(:,:,:,g)=temp(23:52,16:43,29:48);    
end
for g=1:16
    dsp(myo_temp(:,:,:,g)>.6*myo_temp(:,:,:,g));pause;close
end
save myo_template myo_temp
% temp=reshape(temp,64,64,64);
% temp=temp(23:52,16:43,29:48);
% %get rid of liver
% temp_sift=zeros(30,28);temp_sift(13:20,6:9)=1;
% temp(:,:,1)=temp(:,:,1).*temp_sift;
% temp_sift=zeros(30,28);temp_sift(9:23,5:11)=1;temp_sift(18:22,12:13)=1;
% temp(:,:,2)=temp(:,:,2).*temp_sift;
% temp_sift=zeros(30,28);temp_sift(8:24,4:22)=1;temp_sift(8,16:18)=0;
% temp(:,:,3)=temp(:,:,3).*temp_sift;
% temp_sift=zeros(30,28);temp_sift(7:25,4:23)=1;
% temp(:,:,4)=temp(:,:,4).*temp_sift;
% temp_sift=zeros(30,28);temp_sift(4:25,4:23)=1;
% for s=5:20
%     temp(:,:,s)=temp(:,:,s).*temp_sift;
% end
% dsp(temp>max(temp(:))*.5)
% ind_myo=find(temp>max(temp(:))*.5);
% save ncat_g1myo_index ind_myo

for g=1
    filename=['C:\simind\ncat16gles_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    templ=fread(fid,64^3,'single');
    fclose(fid);
end
templ=reshape(templ,64,64,64);
templ=templ(23:52,16:43,29:48);
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);

clear Im_ideal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%snr
%lesion
cat_name={'Im_ncatles_nAS_s';'Im_ncatles_A_s';'Im_ncatles_AS_s'};
snr_name1={'snr_s';'snr_sac';'snr_sasc'};
snr_name23={'snr_t';'snr_tac';'snr_tasc'};
snr_les=zeros(3,5,3,30);
t_index=[0 2 3];
for m=1:3
    for t=1:3
        for k=1:30
            filename=[cat_name{m} num2str(5) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
            if t==1
                load(filename,snr_name1{m});
                snr_les(m,:,t,k)=eval(snr_name1{m});
            else
                load(filename,snr_name23{m});
                snr_les(m,:,t,k)=eval(snr_name23{m});
            end            
        end
    end
end
%normal
cat_name={'./data_16g_5st0/Im_ncat16g_nAS_s';'./data_16g_5st0/Im_ncat16g_A_s';'./data_16g_5st0/Im_ncat16g_AS_s'};
snr_name1={'snr_s';'snr_sac';'snr_sasc'};
snr_name23={'snr_t';'snr_tac';'snr_tasc'};
snr_16g=zeros(3,5,1,30);
t_index=[0 2 3];
for m=1:3
    for t=1:3
        for k=1:30
            filename=[cat_name{m} num2str(5) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
            if t==1
                load(filename,snr_name1{m});
                snr_16g(m,:,t,k)=eval(snr_name1{m});
            else
                load(filename,snr_name23{m});
                snr_16g(m,:,t,k)=eval(snr_name23{m});
            end            
        end
    end
end

%myocardium pixels
load ncat_g1myo_index ind_myo
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
dir_name='D:\imagereconstruction\4D\data_les_5st';
cat_name={'\Im_ncatles_nAS_s';...
    '\Im_ncatles_A_s';...
    '\Im_ncatles_AS_s'};
snr_les_myo=zeros(3,5,3,30);
t_index=[0 2 3];
for m=1:3
    for s=1:5
        for t=1:3
            for k=1:30
                filename=[dir_name num2str(t_index(t)) cat_name{m}...
                    num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
                load(filename,'Im_maps');
                Imroi=Im_maps(:,:,:,1);
                snr_les_myo(m,s,t,k)=10*log10(ncatg1roi(ind_myo)'*ncatg1roi(ind_myo)...
                    /sum((ncatg1roi(ind_myo)-Imroi(ind_myo)).^2));
            end
        end
    end
end
save snr_les_myo snr_les_myo
%fbp
load ncat_g1myo_index ind_myo
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
dir_name='D:\imagereconstruction\4D\data_les_fbp\Im_fbp_n';
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    Imroi=Im_fbp(:,:,:,1);
    snr_lm_fbp(n)=10*log10(ncatg1roi(ind_myo)'*ncatg1roi(ind_myo)...
        /sum((ncatg1roi(ind_myo)-Imroi(ind_myo)).^2));
end
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
    Im_st121=reshape(Im_st121,[30 28 20 16]);
    Imroi=Im_st121(:,:,:,1);
    snr_lm_st121(n)=10*log10(ncatg1roi(ind_myo)'*ncatg1roi(ind_myo)...
        /sum((ncatg1roi(ind_myo)-Imroi(ind_myo)).^2));
end

%normal gate 1
load ncat_g1myo_index ind_myo
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
dir_name='D:\imagereconstruction\4D\data_16g_5st';
cat_name={'\Im_ncat16g_nAS_s';...
    '\Im_ncat16g_A_s';...
    '\Im_ncat16g_AS_s'};
snr_16g_myo=zeros(3,5,3,30);
t_index=[0 2 3];
for m=1:3
    for s=1:5
        for t=1:3
            for k=1:30
                filename=[dir_name num2str(t_index(t)) cat_name{m}...
                    num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
                load(filename,'Im_maps');
                Imroi=Im_maps(:,:,:,1);
                snr_16g_myo(m,s,t,k)=10*log10(ncatg1roi(ind_myo)'*ncatg1roi(ind_myo)...
                    /sum((ncatg1roi(ind_myo)-Imroi(ind_myo)).^2));
            end
        end
    end
end
save snr_16g_myo snr_16g_myo
%fbp
load ncat_g1myo_index ind_myo
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
dir_name='D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n';
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    Imroi=Im_fbp(:,:,:,1);
    snr_nm_fbp(n)=10*log10(ncatg1roi(ind_myo)'*ncatg1roi(ind_myo)...
        /sum((ncatg1roi(ind_myo)-Imroi(ind_myo)).^2));
end
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
    Im_st121=reshape(Im_st121,[30 28 20 16]);
    Imroi=Im_st121(:,:,:,1);
    snr_nm_st121(n)=10*log10(ncatg1roi(ind_myo)'*ncatg1roi(ind_myo)...
        /sum((ncatg1roi(ind_myo)-Imroi(ind_myo)).^2));
end


%SC optimal CF
for n=2:10
    filename=['D:\imagereconstruction\4D\optSC\Im_optSC_CF10_n' num2str(n) '.mat'];
    load(filename,'snr_tasc');
    snr_mn(n-1,:)=snr_tasc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TAC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHO
nor_dir_name={'data_16g_5st0';'data_16g_5st2';'data_16g_5st3'};
les_dir_name={'data_les_5st0';'data_les_5st2';'data_les_5st3'};
wmap_les=zeros(64,64,30);wmap_nles=wmap_les;
for n=1:30
    filename=['I:\Work\R1_Data\im4dData\Data_Final\def\t0\Imrecovdef_s4e-005_t0_' num2str(n) '.mat'];
    load(filename);
    wmap_les(:,:,n)=Imrecon(:,:,37,1);
    filename=['I:\Work\R1_Data\im4dData\Data_Final\nodef\t0\Imrecov_s4e-005_t0_' num2str(n) '.mat'];
    load(filename);
    wmap_nles(:,:,n)=Imrecon(:,:,37,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%images
%normal
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_ideal=Im_ideal(23:52,16:43,29:48,:);

load('D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_nAS_s1t2_n1.mat')
Im_t2n=Im_maps;
load('D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_A_s1t2_n1.mat')
Im_t2a=Im_maps;
load('D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_AS_s1t2_n1.mat')
Im_t2sc=Im_maps;
filename=['D:\imagereconstruction\4D\optSC\Im_optSC_CF1_n' num2str(2) '.mat'];
load(filename,'Im_maps');
Im_t2sco=Im_maps;
load('D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n1.mat');
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end
Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
Im_st121=reshape(Im_st121,[30 28 20 16]);
for n=1:20
    dsp(Im_st121(:,:,n,:),1);title(['frame_' num2str(n)]);
    pause;close
end
%lesion
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_ideal=Im_ideal(23:52,16:43,29:48,:);
for n=1:20
    dsp(Im_ideal(:,:,n,:),1);title(['frame_' num2str(n)]);
    pause;close
end
load('D:\imagereconstruction\4D\data_les_5st2\Im_ncatles_A_s1t2_n1.mat')
Im_t2a=Im_maps;

load('D:\imagereconstruction\4D\data_les_fbp\Im_fbp_n1.mat');
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end
Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
Im_st121=reshape(Im_st121,[30 28 20 16]);
for n=1:20
    dsp(Im_st121(:,:,n,:),1);title(['frame_' num2str(n)]);
    pause;close
end

%short axis Images
load ncat_g1myo_index ind_myo
g=1;
filename=['C:\simind\ncatles2_act_' num2str(g) '.bin'];
fid=fopen(filename,'rb');
temp=fread(fid,'single');
fclose(fid);
temp=temp*5e5/sum(temp(:));
temp=reshape(temp,[64 64 64]);
temp=temp(23:52,16:43,29:48);
im_max=max(temp(ind_myo));
dsp(temp,1,2)
id_short_les=cardiac3drot(temp);
dsp(id_short_les(:,:,15:23),1,2)

filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
fid=fopen(filename,'rb');
temp=fread(fid,'single');
fclose(fid);
temp=temp*5e5/sum(temp(:));
temp=reshape(temp,[64 64 64]);
temp=temp(23:52,16:43,29:48);
im_max=max(temp(ind_myo));
dsp(temp,1,2)
id_short=cardiac3drot(temp);
dsp(id_short(:,:,15:23),1,2)

filename=['C:\simind\ncat16gles_act_' num2str(g) '.bin'];
fid=fopen(filename,'rb');
temp=fread(fid,'single');
fclose(fid);
temp=temp*5e5/sum(temp(:));
temp=reshape(temp,[64 64 64]);
temp=temp(23:52,16:43,29:48);
im_max=max(temp(ind_myo));
dsp(temp,1,2)
id_short_les=cardiac3drot(temp);
dsp(id_short_les(:,:,15:23),1,2)


%Transversal slice
s=1;t=2;n=1;frame=9;
filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_nAS_s'...
    num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
load(filename,'Im_maps');
temp=Im_maps(:,:,frame,1);
dsp(temp,1)

filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_AS_s'...
    num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
load(filename,'Im_maps');
temp=Im_maps(:,:,frame,1);
dsp(temp,1)

filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_A_s'...
    num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
load(filename,'Im_maps');
temp=Im_maps(:,:,frame,1);
dsp(temp,1)

filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_nAS_s'...
    num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
load(filename,'Im_maps');
temp=Im_maps(:,:,frame,1);
dsp(temp,1)

filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_AS_s'...
    num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
load(filename,'Im_maps');
temp=Im_maps(:,:,frame,1);
dsp(temp,1)

filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_A_s'...
    num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
load(filename,'Im_maps');
temp=Im_maps(:,:,frame,1);
dsp(temp,1)
