%data for Umass to view motion

%defect
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,:);
save D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_ideal Im_part -v6

% %binary
% for g=1:16
%     temp=single(Im_part(:,:,:,g));
%     filename=['d:\data4umass\bin\def\ncat_def_ideal' num2str(g) '.bin'];
%     fid=fopen(filename,'wb');
%     fwrite(fid,temp,'single');
%     fclose(fid);
% end


w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
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

load(filename{1});
Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
Im_st121=reshape(Im_st121,[30 28 20 16]);
Im_part=single(Im_st121);
save D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_st121 Im_part -v6

sname={['D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_mapsNC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_mapsAC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_mapsASC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_maptNC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_maptAC.mat'];...    
    ['D:\imagereconstruction\4D\data4umass\mat\defect\ncat_def_maptASC.mat']};

for n=1:6
    load(filename{n+1});
    Im_part=Im_maps;
    save(sname{n},'Im_part','-v6');
end

%normal
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,:);
save D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_ideal Im_part -v6

% %binary
% for g=1:16
%     temp=single(Im_part(:,:,:,g));
%     filename=['d:\data4umass\bin\def\ncat_def_ideal' num2str(g) '.bin'];
%     fid=fopen(filename,'wb');
%     fwrite(fid,temp,'single');
%     fclose(fid);
% end


w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
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

load(filename{1});
Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
Im_st121=reshape(Im_st121,[30 28 20 16]);
Im_part=single(Im_st121);
save D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_st121 Im_part -v6

sname={['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_mapsNC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_mapsAC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_mapsASC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_maptNC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_maptAC.mat'];...    
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_maptASC.mat']};

for n=1:6
    load(filename{n+1});
    Im_part=Im_maps;
    save(sname{n},'Im_part','-v6');
end

%view
lname={['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_ideal.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_st121.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_mapsNC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_mapsAC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_mapsASC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_maptNC.mat'];...
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_maptAC.mat'];...    
    ['D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_maptASC.mat']};

p_slice=18;M=20;G=16;b_slice=8;
Im_comb=zeros(M,M,G,b_slice);
cc
temp=Im_comb(:,:,:,1);immax=max(temp(:));
Im_comb(Im_comb>immax)=immax;

for k=1:8
    temp=Im_comb(:,:,:,k);
    Im_comb(:,:,:,k)=temp/max(temp(:));
end

Im_combs=zeros((M+2)*2,(M+2)*4,G);
im_index=[1 2 3 6 4 7 5 8];
for g=1:G
    for n=1:2
        for k=1:4
            Im_combs((n-1)*(M+2)+1:n*M+(n-1)*2,(k-1)*(M+2)+1:k*M+(k-1)*2,g)=Im_comb(:,:,g,im_index((k-1)*2+n));
        end
    end
end

for g=1:G
    figure('position',[100 100 800 400]),imagesc(interp2(Im_combs(:,:,g),2)),clinicalcolor,axis image
    axis off,im(g)=getframe;
    close
end
figure('position',[100 100 800 400]),movie(im,10)

%%%%%%%%%%%%%%%%%%%
%get full 64x64

dir_name='D:\imagereconstruction\4D\data4umass\mat\normal\ncat_nor_';
filename={'ideal';'mapsNC';'mapsAC';'mapsASC';...
    'maptNC';'maptAC';'maptASC';'st121'};

for n=1:8
    l_name=[dir_name filename{n}];
    load(l_name)
    Im_full=zeros(64,64,64,16);Im_full=single(Im_full);
    Im_full(23:52,16:43,29:48,:)=Im_part;
    save([l_name 'full'],'Im_full','-v6');
end