% load myo_template
% temp=myo_temp(:,:,8,1)';
% myo_ind=find(temp>0);%.7*max(temp(:)));

temp=zeros(size(Im_part));temp(11:12,24:25)=1;def_roi=find(temp>0);
temp=zeros(size(Im_part));temp(9:10,24:25)=1;nor_roi=find(temp>0);%apex
% temp=zeros(size(Im_part));temp(14-2:15-2,23:24)=1;def_roi=find(temp>0);
% temp=zeros(size(Im_part));temp(16-2:17-2,22:23)=1;nor_roi=find(temp>0);%base

load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
gate=1;
Im_part=Im_ideal(23:52,16:43,29:48,gate);
bv_slice=8;
Im_part=Im_part(:,:,bv_slice)';
figure('position',[100 100 500 500]),imagesc(Im_part),colormap(gray),axis image,axis off
hold on,plot(24,11:12,'kx',25,11:12,'kx',...
    24,9:10,'ko',25,9:10,'ko',...
    23,14:15,'k*',24,14:15,'k*',...
    22,16:17,'k^',23,16:17,'k^',...
    20,18:20,'k+',21,18:20,'k+',22,18:20,'k+')
temp=Im_part>0.7*max(Im_part(:));
temp(13,24)=0;
myo_ind=find(temp>0);

%st121
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end
dir_name='D:\imagereconstruction\4D\data_les_fbp\Im_fbp_n';
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
    Im_st121=reshape(Im_st121,[30 28 20 16]);
    Im_st121=Im_st121(:,:,:,gate);
%     immax=max(Im_st121(ind_myo));%immax=1;
    Im_part=Im_st121;
    Im_part=Im_part(:,:,bv_slice)';
    st121_def(n)=mean(Im_part(def_roi));
    st121_nor(n)=mean(Im_part(nor_roi));
    st121_std(n)=std(Im_part(myo_ind));
end

%MAP
dir_name='D:\imagereconstruction\4D\data_les_5st';
cat_name={'\Im_ncatles_nAS_s';...
    '\Im_ncatles_A_s';...
    '\Im_ncatles_AS_s'};
t_index=[0 2 3];
for m=1:3
    for s=1:5
        for t=1:3
            for k=1:30
                filename=[dir_name num2str(t_index(t)) cat_name{m}...
                    num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
                load(filename,'Im_maps');
                Im_maps=Im_maps(:,:,:,gate);
                Im_part=Im_maps;
                Im_part=Im_part(:,:,bv_slice)';
                map_def(k,m,s,t)=mean(Im_part(def_roi));
                map_nor(k,m,s,t)=mean(Im_part(nor_roi));
                map_std(k,m,s,t)=std(Im_part(myo_ind));
            end
        end
    end
end

% save CNR_g1_slice8 st121_def st121_nor st121_std map_def map_nor map_std%roi base
% save  CNR_g6_slice8 st121_def st121_nor st121_std map_def map_nor map_std%roi base
% save CNR_g1_slice8apex st121_def st121_nor st121_std map_def map_nor map_std

save CNR_g1_slice8apex_nm4std st121_def st121_nor st121_std map_def map_nor map_std

load CNR_g1_slice8apex
st121_cnr=abs(st121_nor-st121_def)./st121_std;
map_diff=abs(map_nor-map_def);
map_cnr=map_diff./map_std;

st121_cnr_ave=mean(st121_cnr);
map_cnr_ave=squeeze(mean(map_cnr));
st121_cnr_std=std(st121_cnr);
map_cnr_std=squeeze(std(map_cnr,[],1));

%plot
sbeta=[0 5e-4 1e-3 1e-2 3e-2];figure,
errorbar(sbeta,squeeze(map_cnr_ave(1,:,1)),squeeze(map_cnr_std(1,:,1)),'r*-'),hold on
errorbar(sbeta,squeeze(map_cnr_ave(1,:,2)),squeeze(map_cnr_std(1,:,2)),'b+-')
errorbar(sbeta,squeeze(map_cnr_ave(1,:,3)),squeeze(map_cnr_std(1,:,3)),'go-')
plot(sbeta,repmat(st121_cnr_ave,1,5),'k-')
xlabel('\beta_s'),ylabel('CNR'),title('NC')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','ST 121')

sbeta=[0 5e-5 1e-4 1e-3 3e-3];figure
errorbar(sbeta,squeeze(map_cnr_ave(2,:,1)),squeeze(map_cnr_std(2,:,1)),'r*-'),hold on
errorbar(sbeta,squeeze(map_cnr_ave(2,:,2)),squeeze(map_cnr_std(2,:,2)),'b+-')
errorbar(sbeta,squeeze(map_cnr_ave(2,:,3)),squeeze(map_cnr_std(2,:,3)),'go-')
plot(sbeta,repmat(st121_cnr_ave,1,5),'k-')
xlabel('\beta_s'),ylabel('CNR'),title('AC')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','ST 121')

sbeta=[0 5e-4 1e-3 2e-3 4e-3];figure
errorbar(sbeta,squeeze(map_cnr_ave(3,:,1)),squeeze(map_cnr_std(3,:,1)),'r*-'),hold on
errorbar(sbeta,squeeze(map_cnr_ave(3,:,2)),squeeze(map_cnr_std(3,:,2)),'b+-')
errorbar(sbeta,squeeze(map_cnr_ave(3,:,3)),squeeze(map_cnr_std(3,:,3)),'go-')
plot(sbeta,repmat(st121_cnr_ave,1,5),'k-')
xlabel('\beta_s'),ylabel('CNR'),title('ASC')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','ST 121')

%%%%%%%%%%%%%%%%%%
%
% nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
% les_dir_name='D:\imagereconstruction\4D\data_les_5st';
% t_ind=[0 2 3];
% S=5;T=3;N=30;
% ch_nAs_les=zeros(4,N,S,T);ch_nAs_nor=zeros(4,N,S,T);
% s=2;t=1;
% for n=1:N
%     filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_nAS_s'...
%         num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
%     load(filename,'Im_maps');
%     temp=Im_maps(:,:,frame,gate);
%     temp=postproc_cho(temp,140);
%     im_les(:,:,n)=temp;
%     filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_nAS_s'...
%         num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
%     load(filename,'Im_maps');
%     temp=Im_maps(:,:,frame,gate);
%     temp=postproc_cho(temp,140);    
%     im_nor(:,:,n)=temp;
% end
% em_nor_nAS=im_nor;
% em_les_nAS=im_les;