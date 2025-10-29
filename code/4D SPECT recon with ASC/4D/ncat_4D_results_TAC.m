%TACs
%dsp(Im_ideal(:,:,1),1)%[11:12 22:23] fixed ROI
%tip%[17:18,19:20]
%defect ROI: good for paper
tac_slice=8:12;
clear lesion
for n=1:16
    fn=['D:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesion16g_act_' num2str(n) '.bin'];
    fid=fopen(fn,'rb');
    temp=fread(fid,'single');
    fclose(fid);
    temp=temp>0;%.8*max(temp(:));
    temp=reshape(temp,64,64,64);
    temp=temp(23:52,16:43,29:48);%ideal_c=sum(Im_ideal(:));
%     lesion(:,:,:,n)=temp;
    lesion(:,:,n)=temp(:,:,tac_slice);
end
lesion=zeros(30,28,16);
% lesion(21:22,18:19,:)=1;
lesion(20:22,18:20,:)=1;%15:16,10:11,:)bad;%15:16,10:11bad%20:22,18:20good
% lesion(18:20,20:22,:)=1;%not good for ideal recon.
%real phantom has no PVE

% lesion=squeeze(lesion(:,:,3,:));
% lesion=permute(lesion,[2 1 3]);
for n=1:16
    nz_def(n)=nnz(lesion(:,:,n));%max(temp)*.7
end

tac_slice=8;%:12;
%load ncatlestew_Im_idealUm
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg16=Im_ideal(23:52,16:43,29:48,:);
ncatg16=squeeze(ncatg16(:,:,tac_slice,:));
% tac_id=squeeze(sum(sum(sum(ncatg16.*lesion))))./nz_def';
tac_id=squeeze(sum(sum(ncatg16.*lesion)))./nz_def';
plot(tac_id/max(tac_id))

%st 121
w_st121=[0.5 0.25 zeros(1,13) 0.25];
for n=1:16
    wm_st121(:,n)=circshift(w_st121',n-1);
end
dir_name='D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n';
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    Im_st121=reshape(Im_fbp,[30*28*20 16])*wm_st121;
    Im_st121=reshape(Im_st121,[30 28 20 16]);    
    Im_st121=squeeze(Im_st121(:,:,tac_slice,:));
%     tac_st121(n,:)=squeeze(sum(sum(sum(Im_st121.*lesion))))./nz_def';
    tac_st121(n,:)=squeeze(sum(sum(Im_st121.*lesion)))./nz_def';
end

dir_name='D:\imagereconstruction\4D\data_16g_5st';
cat_name={'\Im_ncat16g_nAS_s';...
    '\Im_ncat16g_A_s';...
    '\Im_ncat16g_AS_s'};
t_index=[0 2 3];

%NC map-s
m=1;s=3;t=1;
for k=1:30
    filename=[dir_name num2str(t_index(t)) cat_name{m}...
        num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
    load(filename,'Im_maps');    
    Im_maps=squeeze(Im_maps(:,:,tac_slice,:));
%     tac_maps(k,:)=squeeze(sum(sum(sum(Im_maps.*lesion))))./nz_def';
    tac_Nmaps(k,:)=squeeze(sum(sum(Im_maps.*lesion)))./nz_def';
end
%NC map-t
m=1;s=1;t=2;
for k=1:30
    filename=[dir_name num2str(t_index(t)) cat_name{m}...
        num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
    load(filename,'Im_maps');
    Im_maps=squeeze(Im_maps(:,:,tac_slice,:));
%     tac_mapt(k,:)=squeeze(sum(sum(sum(Im_maps.*lesion))))./nz_def';
    tac_Nmapt(k,:)=squeeze(sum(sum(Im_maps.*lesion)))./nz_def';
end

%A map-s
m=2;s=4;t=1;
for k=1:30
    filename=[dir_name num2str(t_index(t)) cat_name{m}...
        num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
    load(filename,'Im_maps');    
    Im_maps=squeeze(Im_maps(:,:,tac_slice,:));
%     tac_maps(k,:)=squeeze(sum(sum(sum(Im_maps.*lesion))))./nz_def';
    tac_maps(k,:)=squeeze(sum(sum(Im_maps.*lesion)))./nz_def';
end
%A map-t
m=2;s=1;t=2;
for k=1:30
    filename=[dir_name num2str(t_index(t)) cat_name{m}...
        num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
    load(filename,'Im_maps');
    Im_maps=squeeze(Im_maps(:,:,tac_slice,:));
%     tac_mapt(k,:)=squeeze(sum(sum(sum(Im_maps.*lesion))))./nz_def';
    tac_mapt(k,:)=squeeze(sum(sum(Im_maps.*lesion)))./nz_def';
end

%AS map-s
m=3;s=4;t=1;
for k=1:30
    filename=[dir_name num2str(t_index(t)) cat_name{m}...
        num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
    load(filename,'Im_maps');    
    Im_maps=squeeze(Im_maps(:,:,tac_slice,:));
%     tac_ASmaps(k,:)=squeeze(sum(sum(sum(Im_maps.*lesion))))./nz_def';
    tac_ASmaps(k,:)=squeeze(sum(sum(Im_maps.*lesion)))./nz_def';
end
%AS map-t
m=3;s=1;t=2;
for k=1:30
    filename=[dir_name num2str(t_index(t)) cat_name{m}...
        num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
    load(filename,'Im_maps');    
    Im_maps=squeeze(Im_maps(:,:,tac_slice,:));
%     tac_ASmapt(k,:)=squeeze(sum(sum(sum(Im_maps.*lesion))))./nz_def';
    tac_ASmapt(k,:)=squeeze(sum(sum(Im_maps.*lesion)))./nz_def';
end

% save tac_base_f8_30noise tac_id tac_st121 tac_Nmaps tac_Nmapt tac_maps tac_mapt tac_ASmaps tac_ASmapt%best SNR
% save tac_base_f8_30noiseMaxReg tac_id tac_st121 tac_Nmaps tac_Nmapt tac_maps tac_mapt tac_ASmaps tac_ASmapt
save tac_base_f8_30noiseSROI tac_id tac_st121 tac_Nmaps tac_Nmapt tac_maps tac_mapt tac_ASmaps tac_ASmapt

load tac_base_f8_30noiseMaxReg
%noise 1
noise_num=10;
figure,plot(1:16,tac_id/max(tac_id),'-',...
    1:16,tac_Nmaps(noise_num,:)/max(tac_Nmaps(noise_num,:)),'--',...
    1:16,tac_Nmapt(noise_num,:)/max(tac_Nmapt(noise_num,:)),'-.',...
    1:16,tac_st121(noise_num,:)/max(tac_st121(noise_num,:)),':','linewidth',1.5)
title('NC'),xlabel('Gate number'),ylabel('Normalized activity')
axis([1 16 0.3 1]),
legend('Ideal','MAP-S','MAP-T','ST 121')

figure,plot(1:16,tac_id/max(tac_id),'-',...
    1:16,tac_maps(noise_num,:)/max(tac_maps(noise_num,:)),'--',...
    1:16,tac_mapt(noise_num,:)/max(tac_mapt(noise_num,:)),'-.',...
    1:16,tac_st121(noise_num,:)/max(tac_st121(noise_num,:)),':','linewidth',1.5)
title('AC'),xlabel('Gate number'),ylabel('Normalized activity')
axis([1 16 0.3 1]),
legend('Ideal','MAP-S','MAP-T','ST 121')

figure,plot(1:16,tac_id/max(tac_id),'-',...
    1:16,tac_ASmaps(noise_num,:)/max(tac_ASmaps(noise_num,:)),'--',...
    1:16,tac_ASmapt(noise_num,:)/max(tac_ASmapt(noise_num,:)),'-.',...
    1:16,tac_st121(noise_num,:)/max(tac_st121(noise_num,:)),':','linewidth',1.5)
% title([num2str(noise_num) 'ASC'])
title('ASC'),xlabel('Gate number'),ylabel('Normalized activity')
axis([1 16 0.3 1]),
legend('Ideal','MAP-S','MAP-T','ST 121')


%mean tac
figure,plot(1:16,tac_id/max(tac_id),'-',1:16,mean(tac_st121,1)/max(mean(tac_st121,1)),':',...
    1:16,mean(tac_Nmaps,1)/max(mean(tac_Nmaps,1)),'--',1:16,mean(tac_Nmapt,1)/max(mean(tac_Nmapt,1)),'-.','linewidth',1.5)
title('NC'),xlabel('Gate number'),ylabel('Normalized activity')
axis([1 16 0.5 1]),legend('Ideal','ST 121','MAP-S','MAP-T')

figure,plot(1:16,tac_id/max(tac_id),'-',1:16,mean(tac_st121,1)/max(mean(tac_st121,1)),':',...
    1:16,mean(tac_maps,1)/max(mean(tac_maps,1)),'--',1:16,mean(tac_mapt,1)/max(mean(tac_mapt,1)),'-.','linewidth',1.5)
title('A'),xlabel('Gate number'),ylabel('Normalized activity')
axis([1 16 0.5 1]),legend('Ideal','ST 121','MAP-S','MAP-T')

figure,plot(1:16,tac_id/max(tac_id),'-',1:16,mean(tac_st121,1)/max(mean(tac_st121,1)),':',...
    1:16,mean(tac_ASmaps,1)/max(mean(tac_ASmaps,1)),'--',1:16,mean(tac_ASmapt,1)/max(mean(tac_ASmapt,1)),'-.','linewidth',1.5)
title('AS'),xlabel('Gate number'),ylabel('Normalized activity')
axis([1 16 0.5 1]),legend('Ideal','ST 121','MAP-S','MAP-T')


tac_rec=cat(2,tac_id,mean(tac_st121)',mean(tac_maps)',mean(tac_mapt)');
corrcoef(tac_rec)
for n=1:30
    % figure,plot(1:16,tac_id/max(tac_id),'-',1:16,tac_st121(n,:)/max(tac_st121(n,:)),':',...
    %     1:16,tac_maps(n,:)/max(tac_maps(n,:)),'--',1:16,tac_mapt(n,:)/max(tac_mapt(n,:)),'-.')
    tac_recN=cat(2,tac_id,tac_st121(n,:)',tac_Nmaps(n,:)',tac_Nmapt(n,:)');
    tac_recA=cat(2,tac_id,tac_st121(n,:)',tac_maps(n,:)',tac_mapt(n,:)');
    tac_recAS=cat(2,tac_id,tac_st121(n,:)',tac_ASmaps(n,:)',tac_ASmapt(n,:)');
    temp=corrcoef(tac_recN);
    ccN(n,:)=temp(1,2:4);
    temp=corrcoef(tac_recA);
    ccA(n,:)=temp(1,2:4);
    temp=corrcoef(tac_recAS);
    ccAS(n,:)=temp(1,2:4);
end

save tac_base_f8 tac_recN tac_recA tac_recAS


