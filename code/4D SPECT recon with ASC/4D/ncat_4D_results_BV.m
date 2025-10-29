load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
gate=1;
Im_part=Im_ideal(23:52,16:43,29:48,gate);
% load ncat_g1myo_index ind_myo
% immax=max(Im_part(ind_myo));immax=1;
% Im_part=Im_part/immax;

% Im_part=Im_part(:,:,8,:);
% Im_part=permute(squeeze(Im_part),[2 1 3]);
% dsp(Im_part,1,0)
%slice 8: defect [11:12 24:25];normal [9:10 24:25]: good
%slice 9: defect [10:11 24:25];normal [8:9 24:25]: no good
bv_slice=8;
Im_part=Im_part(:,:,bv_slice)';tot_g1s8=sum(Im_part(:));
figure('position',[100 100 500 500]),imagesc(Im_part),colormap(gray)
% hold on,plot(24,11:12,'kx',25,11:12,'kx',...
%     24,9:10,'ko',25,9:10,'ko',20,18:20,'k+',21,18:20,'k+',22,18:20,'k+')
hold on,plot(23,14:15,'k*',24,14:15,'k*',...
     22,16:17,'k^',23,16:17,'k^',20,18:20,'k+',21,18:20,'k+',22,18:20,'k+')
% hold on,plot(23,14:15,'kx',24,14:15,'kx',...
%     22,16:17,'ko',23,16:17,'ko',18,20:22,'w+',19,20:22,'w+',20,20:22,'w+','linewidth',2)
axis image,axis off,title('x: Defect; o: Normal; +: TAC')

% temp=zeros(size(Im_part));temp(11:12,24:25)=1;def_roi=find(temp>0);
% temp=zeros(size(Im_part));temp(9:10,24:25)=1;nor_roi=find(temp>0);
temp=zeros(size(Im_part));temp(14:15,23:24)=1;def_roi=find(temp>0);
temp=zeros(size(Im_part));temp(16:17,22:23)=1;nor_roi=find(temp>0);%roi2

id_def=mean(Im_part(def_roi));
id_nor=mean(Im_part(nor_roi));

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
    Im_part=Im_part*tot_g1s8/sum(Im_part(:));
    st121_def(n)=mean(Im_part(def_roi));
    st121_nor(n)=mean(Im_part(nor_roi));
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
%                 immax=max(Im_maps(ind_myo));%immax=1;
%                 Im_part=Im_maps/immax;
                Im_part=Im_maps;
                Im_part=Im_part(:,:,bv_slice)';
                Im_part=Im_part*tot_g1s8/sum(Im_part(:));
                map_def(k,m,s,t)=mean(Im_part(def_roi));
                map_nor(k,m,s,t)=mean(Im_part(nor_roi));
            end
        end
    end
end
% save BV_slice8 id_def id_nor st121_def st121_nor map_def map_nor
% save BV_g1_slice8 id_def id_nor st121_def st121_nor map_def map_nor%normalized
% load BV_g1_slice8

% save BV_g1_slice8_org id_def id_nor st121_def st121_nor map_def map_nor%not normalized (apex)
% save BV_g1_slice8_rn id_def id_nor st121_def st121_nor map_def map_nor%normalization on current slice (apex)
% save BV_g1_slice8_roi2 id_def id_nor st121_def st121_nor map_def map_nor%not normalized (base)
save BV_g1_slice8_roi2N id_def id_nor st121_def st121_nor map_def map_nor%normalized on current slice (base)

% load BV_g1_slice8_rn
load BV_g1_slice8_roi2
bias_d_st121=abs(id_def-mean(st121_def))/id_def;
bias_n_st121=abs(id_nor-mean(st121_nor))/id_nor;
std_d_st121=std(st121_def)/id_def;
std_n_st121=std(st121_nor)/id_nor;

bias_d_mapnAS=abs(id_def-squeeze(mean(map_def(:,1,:,:),1)))/id_def;
bias_n_mapnAS=abs(id_nor-squeeze(mean(map_nor(:,1,:,:),1)))/id_nor;
std_d_mapnAS=squeeze(std(map_def(:,1,:,:),[],1))/id_def;
std_n_mapnAS=squeeze(std(map_nor(:,1,:,:),[],1))/id_nor;

bias_d_mapA=abs(id_def-squeeze(mean(map_def(:,2,:,:),1)))/id_def;
bias_n_mapA=abs(id_nor-squeeze(mean(map_nor(:,2,:,:),1)))/id_nor;
std_d_mapA=squeeze(std(map_def(:,2,:,:),[],1))/id_def;
std_n_mapA=squeeze(std(map_nor(:,2,:,:),[],1))/id_nor;

bias_d_mapAS=abs(id_def-squeeze(mean(map_def(:,3,:,:),1)))/id_def;
bias_n_mapAS=abs(id_nor-squeeze(mean(map_nor(:,3,:,:),1)))/id_nor;
std_d_mapAS=squeeze(std(map_def(:,3,:,:),[],1))/id_def;
std_n_mapAS=squeeze(std(map_nor(:,3,:,:),[],1))/id_nor;

figure,plot(bias_d_mapnAS(:,1)*100,std_d_mapnAS(:,1)*100,'+-',...
    bias_d_mapnAS(:,2)*100,std_d_mapnAS(:,2)*100,'*-',...
    bias_d_mapnAS(:,3)*100,std_d_mapnAS(:,3)*100,'.-')
hold on,plot(bias_d_st121*100,std_d_st121*100,'^')
xlabel('Bias(%)'),ylabel('Std(%)'),title('NC Defect')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','ST 121')
axis([10 45 0 40])

figure,plot(bias_n_mapnAS(:,1)*100,std_n_mapnAS(:,1)*100,'+-',...
    bias_n_mapnAS(:,2)*100,std_n_mapnAS(:,2)*100,'*-',...
    bias_n_mapnAS(:,3)*100,std_n_mapnAS(:,3)*100,'.-')
hold on,plot(bias_n_st121*100,std_n_st121*100,'^')
xlabel('Bias(%)'),ylabel('Std(%)'),title('NC Normal')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','ST 121')
axis([20 52 0 30])

figure,plot(bias_d_mapA(:,1)*100,std_d_mapA(:,1)*100,'+-',...
    bias_d_mapA(:,2)*100,std_d_mapA(:,2)*100,'*-',...
    bias_d_mapA(:,3)*100,std_d_mapA(:,3)*100,'.-')
hold on,plot(bias_d_st121*100,std_d_st121*100,'^')
xlabel('Bias(%)'),ylabel('Std(%)'),title('AC Defect')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','ST 121')
axis([10 45 0 40])

figure,plot(bias_n_mapA(:,1)*100,std_n_mapA(:,1)*100,'+-',...
    bias_n_mapA(:,2)*100,std_n_mapA(:,2)*100,'*-',...
    bias_n_mapA(:,3)*100,std_n_mapA(:,3)*100,'.-')
hold on,plot(bias_n_st121*100,std_n_st121*100,'^')
xlabel('Bias(%)'),ylabel('Std(%)'),title('AC Normal')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','ST 121')
axis([20 52 0 30])

figure,plot(bias_d_mapAS(:,1)*100,std_d_mapAS(:,1)*100,'+-',...
    bias_d_mapAS(:,2)*100,std_d_mapAS(:,2)*100,'*-',...
    bias_d_mapAS(:,3)*100,std_d_mapAS(:,3)*100,'.-')
hold on,plot(bias_d_st121*100,std_d_st121*100,'^')
xlabel('Bias(%)'),ylabel('Std(%)'),title('ASC Defect')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','ST 121')
axis([10 45 0 40])

figure,plot(bias_n_mapAS(:,1)*100,std_n_mapAS(:,1)*100,'+-',...
    bias_n_mapAS(:,2)*100,std_n_mapAS(:,2)*100,'*-',...
    bias_n_mapAS(:,3)*100,std_n_mapAS(:,3)*100,'.-')
hold on,plot(bias_n_st121*100,std_n_st121*100,'^')
xlabel('Bias(%)'),ylabel('Std(%)'),title('ASC Normal')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','ST 121')
axis([20 52 0 30])



% axis([0 0.5 0 0.5])
