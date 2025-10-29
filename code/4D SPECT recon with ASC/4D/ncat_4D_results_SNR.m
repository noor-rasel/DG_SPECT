%data from ncat_4D_results.m
% load snr_les_myo%MAP lesion gate 1: 30 noise
% load snr_16g_myo%MAP normal gate 1: 30 noise
% load snr_myo_fbp121%FBP & ST-121 gate 1 (both lesion and normal): 30 noise
% 
% noise_num=1:30;
% figure,plot(noise_num,snr_lm_fbp,'*-',noise_num,snr_lm_st121,'+-',...
%     noise_num,squeeze(snr_les_myo(2,3,1,:)),'^-',noise_num,squeeze(snr_les_myo(2,1,2,:)),'o-')
% legend('FBP','ST121','MAP-S','MAP-T'),xlabel('Noise number'),ylabel('dB')
% 
% noise_num=1:30;
% figure,plot(noise_num,snr_nm_fbp,'*-',noise_num,snr_nm_st121,'+-',...
%     noise_num,squeeze(snr_16g_myo(2,3,1,:)),'^-',noise_num,squeeze(snr_16g_myo(2,1,2,:)),'o-')
% legend('FBP','ST121','MAP-S','MAP-T'),xlabel('Noise number'),ylabel('dB')

%16 gates
load myo_template
for g=1:16
    ind_myo{g}=find(myo_temp(:,:,:,g)>.6*myo_temp(:,:,:,g));
end
%myocardium pixels
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg16=Im_ideal(23:52,16:43,29:48,:);
dir_name='D:\imagereconstruction\4D\data_les_5st';
cat_name={'\Im_ncatles_nAS_s';...
    '\Im_ncatles_A_s';...
    '\Im_ncatles_AS_s'};
snr_les_myo=zeros(3,5,3,16,30);
t_index=[0 2 3];
for m=1:3
    for s=1:5
        for t=1:3
            for k=1:30
                filename=[dir_name num2str(t_index(t)) cat_name{m}...
                    num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
                load(filename,'Im_maps');
                for g=1:16
                    Imroi=Im_maps(:,:,:,g);
                    ncatg1roi=ncatg16(:,:,:,g);
                    snr_les_myo(m,s,t,g,k)=10*log10(ncatg1roi(ind_myo{g})'*ncatg1roi(ind_myo{g})...
                        /sum((ncatg1roi(ind_myo{g})-Imroi(ind_myo{g})).^2));
                end
            end
        end
    end
end
save snr_les_myo16g snr_les_myo
%fbp
dir_name='D:\imagereconstruction\4D\data_les_fbp\Im_fbp_n';
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    for g=1:16
        Imroi=Im_fbp(:,:,:,g);
        ncatg1roi=ncatg16(:,:,:,g);
        snr_lm_fbp(g,n)=10*log10(ncatg1roi(ind_myo{g})'*ncatg1roi(ind_myo{g})...
            /sum((ncatg1roi(ind_myo{g})-Imroi(ind_myo{g})).^2));
    end
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
    for g=1:16
        Imroi=Im_st121(:,:,:,g);
        ncatg1roi=ncatg16(:,:,:,g);
        snr_lm_st121(g,n)=10*log10(ncatg1roi(ind_myo{g})'*ncatg1roi(ind_myo{g})...
            /sum((ncatg1roi(ind_myo{g})-Imroi(ind_myo{g})).^2));
    end
end
save snr_les_fbp16g snr_lm_fbp snr_lm_st121
%normal gate 1
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg16=Im_ideal(23:52,16:43,29:48,:);
dir_name='D:\imagereconstruction\4D\data_16g_5st';
cat_name={'\Im_ncat16g_nAS_s';...
    '\Im_ncat16g_A_s';...
    '\Im_ncat16g_AS_s'};
snr_nor_myo=zeros(3,5,3,16,30);
t_index=[0 2 3];
for m=1:3
    for s=1:5
        for t=1:3
            for k=1:30
                filename=[dir_name num2str(t_index(t)) cat_name{m}...
                    num2str(s) 't' num2str(t_index(t)) '_n' num2str(k) '.mat'];
                load(filename,'Im_maps');
                for g=1:16
                    Imroi=Im_maps(:,:,:,g);
                    ncatg1roi=ncatg16(:,:,:,g);
                    snr_nor_myo(m,s,t,g,k)=10*log10(ncatg1roi(ind_myo{g})'*ncatg1roi(ind_myo{g})...
                        /sum((ncatg1roi(ind_myo{g})-Imroi(ind_myo{g})).^2));
                end
            end
        end
    end
end
save snr_nor_myo16g snr_nor_myo
%fbp
dir_name='D:\imagereconstruction\4D\data_16g_fbp\Im_fbp_n';
for n=1:30
    filename=[dir_name num2str(n) '.mat'];
    load(filename,'Im_fbp');
    for g=1:16
        Imroi=Im_fbp(:,:,:,g);
        ncatg1roi=ncatg16(:,:,:,g);
        snr_nm_fbp(g,n)=10*log10(ncatg1roi(ind_myo{g})'*ncatg1roi(ind_myo{g})...
            /sum((ncatg1roi(ind_myo{g})-Imroi(ind_myo{g})).^2));
    end
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
    for g=1:16
        Imroi=Im_st121(:,:,:,g);
        ncatg1roi=ncatg16(:,:,:,g);
        snr_nm_st121(g,n)=10*log10(ncatg1roi(ind_myo{g})'*ncatg1roi(ind_myo{g})...
            /sum((ncatg1roi(ind_myo{g})-Imroi(ind_myo{g})).^2));
    end
end
save snr_nor_fbp16g snr_nm_fbp snr_nm_st121

%plot
load snr_nor_myo16g
snr_nor_myo16=mean(mean(snr_nor_myo,5),4);
% snr_les_myo16=mean(mean(snr_les_myo16,5),4);
%NC
sbeta=[0 5e-4 1e-3 1e-2 3e-2];
figure,plot(sbeta,snr_nor_myo16(1,:,1),'*-',sbeta,snr_nor_myo16(1,:,2),'+-',...
    sbeta,snr_nor_myo16(1,:,3),'o-')
hold on,plot(sbeta,repmat(mean(mean(snr_nm_fbp)),1,5),'k-',sbeta,repmat(mean(mean(snr_nm_st121)),1,5),'k--')
axis([0 3e-2 6.5 13])
xlabel('\beta_s'),ylabel('SNR:dB')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','FBP','ST 121')
%A
sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,snr_nor_myo16(2,:,1),'*-',sbeta,snr_nor_myo16(2,:,2),'+-',...
    sbeta,snr_nor_myo16(2,:,3),'o-')
hold on,plot(sbeta,repmat(mean(mean(snr_nm_fbp)),1,5),'k-',sbeta,repmat(mean(mean(snr_nm_st121)),1,5),'k--')
axis([0 3e-3 6.5 13])
xlabel('\beta_s'),ylabel('SNR:dB')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','FBP','ST 121')
%A+S
sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,snr_nor_myo16(3,:,1),'*-',sbeta,snr_nor_myo16(3,:,2),'+-',...
    sbeta,snr_nor_myo16(3,:,3),'o-')
hold on,plot(sbeta,repmat(mean(mean(snr_nm_fbp)),1,5),'k-',sbeta,repmat(mean(mean(snr_nm_st121)),1,5),'k--')
axis([0 4e-3 6.5 13])
xlabel('\beta_s'),ylabel('SNR:dB')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','FBP','ST 121')

%%%%%%%%%%%%%%%
%no fbp
%NC
load snr_nor_myo16g
snr_nor_myo16=mean(mean(snr_nor_myo,5),4);
load snr_nor_fbp16g
sbeta=[0 5e-4 1e-3 1e-2 3e-2];
figure,plot(sbeta,snr_nor_myo16(1,:,1),'*-',sbeta,snr_nor_myo16(1,:,2),'+-',...
    sbeta,snr_nor_myo16(1,:,3),'o-')
hold on,plot(sbeta,repmat(mean(mean(snr_nm_st121)),1,5),'k-')
axis([0 3e-2 6.5 13])
xlabel('\beta_s'),ylabel('SNR:dB')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','ST 121')
%A
sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,snr_nor_myo16(2,:,1),'*-',sbeta,snr_nor_myo16(2,:,2),'+-',...
    sbeta,snr_nor_myo16(2,:,3),'o-')
hold on,plot(sbeta,repmat(mean(mean(snr_nm_st121)),1,5),'k-')
axis([0 3e-3 6.5 13])
xlabel('\beta_s'),ylabel('SNR:dB')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','ST 121')
%A+S
sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,snr_nor_myo16(3,:,1),'*-',sbeta,snr_nor_myo16(3,:,2),'+-',...
    sbeta,snr_nor_myo16(3,:,3),'o-')
hold on,plot(sbeta,repmat(mean(mean(snr_nm_st121)),1,5),'k-')
axis([0 4e-3 6.5 13])
xlabel('\beta_s'),ylabel('SNR:dB')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','ST 121')


%16 gate profile(AS): map-s:2e-3; map-t:2e-3
load snr_nor_myo16g
%AS
snr_nor_myo16=mean(snr_nor_myo,5);
figure,plot(1:16,squeeze(snr_nor_myo16(3,4,1,:)),'*-',1:16,squeeze(snr_nor_myo16(3,1,2,:)),'+-',...
    1:16,mean(snr_nm_st121,2),'o-')
axis([1 16 8.5 13.5])
xlabel('Gate number'),ylabel('SNR:dB')
legend('MAP-S','MAP-T','ST 121')
%A
snr_nor_myo16=mean(snr_nor_myo,5);
figure,plot(1:16,squeeze(snr_nor_myo16(2,4,1,:)),'*-',1:16,squeeze(snr_nor_myo16(2,1,2,:)),'+-',...
    1:16,mean(snr_nm_st121,2),'o-')
% axis([1 16 8.5 13.5])
xlabel('Gate number'),ylabel('SNR:dB')
legend('MAP-S','MAP-T','ST 121')
%nAs
snr_nor_myo16=mean(snr_nor_myo,5);
figure,plot(1:16,squeeze(snr_nor_myo16(1,4,1,:)),'*-',1:16,squeeze(snr_nor_myo16(1,1,2,:)),'+-',...
    1:16,mean(snr_nm_st121,2),'o-')
% axis([1 16 8.5 13.5])
xlabel('Gate number'),ylabel('SNR:dB')
legend('MAP-S','MAP-T','ST 121')