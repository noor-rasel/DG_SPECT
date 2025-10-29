%%%%%%%%%%%%%%%%%%%%%
%Transverse images: slice 36(8 in sub-images). gate 1 and 6
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,[1 6]);
Im_part=permute(squeeze(Im_part(3:30,:,8,:)),[2 1 3]);
immax=max(Im_part(:));
M=28;b_slice=8;
Im_comb=zeros(M,M,2,b_slice);
Im_comb(:,:,:,1)=Im_part;

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
for k=1:7
    load(filename{k})
    if k>1
        Im_part=permute(squeeze(Im_maps(3:30,:,b_slice,[1 6])),[2 1 3]);
        Im_comb(:,:,:,k+1)=Im_part;        
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=permute(squeeze(Im_st121(3:30,:,b_slice,[1 6])),[2 1 3]);
        Im_comb(:,:,:,k+1)=Im_part;        
    end
end
Im_comb(Im_comb>immax)=immax;

Im_combs=zeros(M*2,M*8);
for n=1:2
    for k=1:8
        Im_combs((n-1)*M+1:n*M,(k-1)*M+1:k*M)=Im_comb(:,:,n,k);
    end
end

figure('position',[100 100 800 300]),imagesc(interp2(Im_combs,2)),clinicalcolor,axis image
set(gca,'Xtick',[221+55:111:864],'Ytick',[55:110:216],'XTickLabel',...
    {'NC';'AC';'ASC';'NC';'AC';'ASC'},...
    'YTickLabel',{'ED';'ES'});
hold on, plot(110,1:221,'w.',221,1:221,'w.',553,1:221,'w.')

%%%%%%%%%%%%%%%%%%%%%%%
%short axis images
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,[1 6]);
p_slice=18;M=20;
for n=1:2
    temp=Im_part(:,:,:,n);
    temp=cardiac3drot(temp);
    id_short(:,:,n)=temp(:,:,p_slice);
end
immax=max(id_short(:));
Im_comb=zeros(M,M,2,b_slice);
Im_comb(:,:,:,1)=id_short(7:26,:,:);

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
for k=1:7
    load(filename{k})
    if k>1
        Im_part=Im_maps(:,:,:,[1 6]);
        for n=1:2
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short(7:26,:,:);
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=Im_st121(:,:,:,[1 6]);
        for n=1:2
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short(7:26,:,:);
    end
end
Im_comb(Im_comb>immax)=immax;
save reim_short_p18 Im_comb

Im_combs=zeros((M+2)*2,(M+2)*8);
for n=1:2
    for k=1:8
        Im_combs((n-1)*(M+2)+1:n*M+(n-1)*2,(k-1)*(M+2)+1:k*M+(k-1)*2)=Im_comb(:,:,n,k);
    end
end

figure('position',[100 100 800 300]),imagesc(interp2(Im_combs,2)),clinicalcolor,axis image
set(gca,'Xtick',[220:88:701],'Ytick',[44:88:173],'XTickLabel',...
    {'NC';'AC';'ASC';'NC';'AC';'ASC'},...
    'YTickLabel',{'ED';'ES'});
hold on, plot(83,1:221,'w.',168,1:221,'w.',438,1:221,'w.')



%%%%%%%%%%%%%%%%%%%%%
%short axis profile
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,[1 6]);

p_slice=18;p_line=17;%16-21 for slice (18)
id_short_p=zeros(20,length(p_slice),2);
for n=1:2
    temp=Im_part(:,:,:,n);
    id_short(:,:,:,n)=cardiac3drot(temp);
    id_short_p(:,:,n)=squeeze(id_short(p_line,:,p_slice,n));
end
id_short_p=squeeze(id_short_p);
% for p=1:6
%     figure,plot(squeeze(id_short_p(:,p,:))),legend('1','6')
% end

p_slice=18;p_line=17;%16-21 for slice (18)
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end

for noise_num=1:30;%2 is good
filename={['D:\imagereconstruction\4D\data_16g_5st0\Im_ncat16g_A_s3t0_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_5st2\Im_ncat16g_A_s1t2_n' num2str(noise_num) '.mat'];...
    ['D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n' num2str(noise_num) '.mat']};

re_short_p=zeros(20,length(p_slice),2,3);
for k=1:3
    load(filename{k})
    if k<3
        Im_part=Im_maps(:,:,:,[1 6]);
        for n=1:2
            temp=Im_part(:,:,:,n);
            re_short(:,:,:,n)=cardiac3drot(temp);
            re_short_p(:,:,n,k)=squeeze(re_short(p_line,:,p_slice,n));
        end
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=Im_st121(:,:,:,[1 6]);
        for n=1:2
            temp=Im_part(:,:,:,n);
            re_short(:,:,:,n)=cardiac3drot(temp);
            re_short_p(:,:,n,k)=squeeze(re_short(p_line,:,p_slice,n));
        end
    end
end
re_short_p_30noise(:,:,:,noise_num)=squeeze(re_short_p);
end
% save 4D_short_lineprof re_short_p_30noise id_short_p%AS
save 4D_shortA_lineprof re_short_p_30noise id_short_p
%plot
pt=1:20;
frame_name={'ED';'ES'};
re_short_p=mean(re_short_p_30noise,4);
for n=1:2
    figure,plot(pt,id_short_p(:,n),'-',pt,re_short_p(:,n,3),':',...
        pt,re_short_p(:,n,1),'-.',pt,re_short_p(:,n,2)),'--';
    title(frame_name{n}),ylabel('Activity'),legend('Ideal','ST-121','MAP-S','MAP-T')
end
for noise_num=1:30
for n=1:3
    cc1=corrcoef(id_short_p(:,1),re_short_p_30noise(:,1,n,noise_num));
    cc_ed(n,noise_num)=cc1(2);
    cc2=corrcoef(id_short_p(:,2),re_short_p_30noise(:,2,n,noise_num));
    cc_es(n,noise_num)=cc2(2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%all gates show
b_slice=8;
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_part=Im_ideal(23:52,16:43,29:48,:);

p_slice=18;M=20;G=16;
for n=1:G
    temp=Im_part(:,:,:,n);
    temp=cardiac3drot(temp);
    id_short(:,:,n)=temp(:,:,p_slice);
end
immax=max(id_short(:));
Im_comb=zeros(M,M,G,b_slice);
Im_comb(:,:,:,1)=id_short(7:26,:,:);

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
            temp=cardiac3drot(temp);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short(7:26,:,:);
    else
        Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
        Im_st121=reshape(Im_st121,[30 28 20 16]);
        Im_part=Im_st121(:,:,:,:);
        for n=1:G
            temp=Im_part(:,:,:,n);
            temp=cardiac3drot(temp);
            id_short(:,:,n)=temp(:,:,p_slice);
        end
        Im_comb(:,:,:,k+1)=id_short(7:26,:,:);
    end
end
Im_comb(Im_comb>immax)=immax;

Im_combs=zeros((M+2)*G,(M+2)*8);
for g=1:G
    for n=1:8
            Im_combs((g-1)*(M+2)+1:g*M+(g-1)*2,(n-1)*(M+2)+1:n*M+(n-1)*2)=Im_comb(:,:,g,n);        
    end
end
save 4D16gate_shortim Im_combs

im_show=zeros((M+2)*8,(M+2)*8);
for n=1:8
    im_show((n-1)*(M+2)+1:n*(M+2),:)=Im_combs((2*(n-1))*(M+2)+1:(2*(n-1)+1)*(M+2),:);
end
im_show=interp2(im_show,2);
figure('position',[100 100 800 800]),imagesc(im_show),clinicalcolor,axis image
set(gca,'Xtick',[220:88:701],'Ytick',[44:88:701],'XTickLabel',...
    {'NC';'AC';'ASC';'NC';'AC';'ASC'},...
    'YTickLabel',{'Gate #1';'Gate #3';'Gate #5';'Gate #7';'Gate #9';'Gate #11';'Gate #13';'Gate #15';});
hold on, plot(83,1:701,'w.',168,1:701,'w.',438,1:701,'w.')
