m_name={'ncatlestew_iMCnASC';'ncatlestew_iMC_AC';'ncatlestew_iMC_ASC';...
    'ncatlestew_mapsnASC';'ncatlestew_maps_AC';'ncatlestew_maps_ASC'};
snr_name={'snr_t';'snr_tac';'snr_tasc';'snr_s';'snr_sac';'snr_sasc'};
%first part
snr_comp1=zeros(6,7);
n=7;
for s=1:6
        filename=['.\isbi07\' m_name{s} num2str(n) '.mat'];        
        load(filename,snr_name{s});
        snr_comp1(s,:)=eval(snr_name{s});
end
%second part
snr_comp2=zeros(6,9);
%snr_sasc 7:16 recomputed
for n=1:9
    filename=['.\isbi07\' m_name{6} num2str(n+7) '.mat'];        
        load(filename,snr_name{6});
        temp=eval(snr_name{6});
        snr_comp2(6,n)=temp(n);
end
n=16;
for s=1:5
        filename=['.\isbi07\' m_name{s} num2str(n) '.mat'];        
        load(filename,snr_name{s});
        snr_comp2(s,:)=eval(snr_name{s});
end
snr_comp=cat(2,snr_comp2,snr_comp1);%7:16 1:6
%temporal and spatial only
figure,plot([0.001:0.001:0.009 0.01:0.01:0.07],snr_comp(1,:),'-',...
    [0.001:0.001:0.009 0.01:0.01:0.07],snr_comp(4,:),'--')
%temporal and spatial with AC
figure,plot([0.0001:0.0001:0.0009 0.001:0.001:0.007],snr_comp(2,:),'-',...
    [0.0001:0.0001:0.0009 0.001:0.001:0.007],snr_comp(5,:),'--')
%temporal and spatial with AC and SC
figure,plot([0.0004:0.0001:0.0009 0.001:0.001:0.003 0.004:0.001:0.01],snr_comp(3,:),'-',...
    [0.0002:0.0001:0.001 0.002:0.001:0.008],snr_comp(6,:),'--')
%TEW
m_name={'ncatlestewE_iMC_ASC';'ncatlestewE_maps_ASC'};
snr_name={'snr_tascE';'snr_sascE'};
n=9;
for s=1:2
    filename=['.\isbi07\' m_name{s} num2str(n) '.mat'];
    load(filename,snr_name{s});
    temp=eval(snr_name{s});
    snr_tew(s,:)=temp;
end


snr_comp=cat(2,snr_comp1,snr_comp2);
[max_snr,max_snr_i]=max(snr_comp,[],2)
%images
load ncatlestew_Im_idealA
Im_ideal=Im_ideal(25:49,16:40,34:38,:);
Im_ideal=squeeze(Im_ideal(:,:,3,:));
Im_ideal=permute(Im_ideal,[2 1 3]);
dsp(Im_ideal)

m_name={'ncatlestew_iMCnASC';'ncatlestew_iMC_AC';'ncatlestew_iMC_ASC';...
    'ncatlestew_mapsnASC';'ncatlestew_maps_AC';'ncatlestew_maps_ASC'};
for s=1:6
    filename=['.\isbi07\' m_name{s} num2str(max_snr_i(s)) '.mat'];
    load(filename,'Im_maps')
    temp=squeeze(Im_maps(:,:,3,:));temp=permute(temp,[2 1 3]);
    Im_comb(:,:,:,s)=temp;
end

m_name={'ncatlestewE_iMC_ASC';'ncatlestewE_maps_ASC'};
ind_tew=[8,9];
for s=1:2
    filename=['.\isbi07\' m_name{s} num2str(ind_tew(s)) '.mat'];
    load(filename,'Im_maps')
    temp=squeeze(Im_maps(:,:,3,:));temp=permute(temp,[2 1 3]);
    Im_comb_tew(:,:,:,s)=temp;
end

%optimal SNR
%load ncatlestew_Im_idealA
load ncatlestew_Im_ideal
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
m_name={'ncatlestew_eMCnASCopt';'ncatlestew_eMC_ACopt';'ncatlestew_eMC_ASCopt';...
    'ncatlestew_mapsnASCopt';'ncatlestew_maps_ACopt';'ncatlestew_maps_ASCopt'};
% m_name={'ncatlestew_eMCnASCopt20';'ncatlestew_eMC_ACopt20';'ncatlestew_eMC_ASCopt20';...
%     'ncatlestew_mapsnASCopt20';'ncatlestew_maps_ACopt20';'ncatlestew_maps_ASCopt20'};
Im_comb=zeros(25,25,16,6);
for s=1:6
    filename=['.\isbi07\' m_name{s} '.mat'];
    load(filename,'Im_maps')
    Im_maps=Im_maps*8e6/sum(Im_maps(:));
    Im_maps=Im_maps(25:49,16:40,34:38,:);
    temp=squeeze(Im_maps(:,:,3,:));temp=permute(temp,[2 1 3]);
    Im_comb(:,:,:,s)=temp;
    Imroi=Im_maps(:,:,:,1);
    snr_opt(s)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
end
%use 10 iterations!!!

%TACs
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
%roi_def=[45:47,26:28]%slice 36
%roi_def=[45:46,25:28]%slice 37
sli=37;
tac_ideal=squeeze(mean(mean(Im_ideal(45:46,25:28,sli,:))));
figure,plot(tac_ideal)
m_name={'ncatlestew_eMCnASCopt';'ncatlestew_eMC_ACopt';'ncatlestew_eMC_ASCopt';...
    'ncatlestew_mapsnASCopt';'ncatlestew_maps_ACopt';'ncatlestew_maps_ASCopt'};
for s=1:6
    filename=['.\isbi07\' m_name{s} '.mat'];
    load(filename,'Im_maps')
    Im_maps=Im_maps*8e6/sum(Im_maps(:));
    tac_rec(:,s)=squeeze(mean(mean(Im_maps(45:46,25:28,sli,:))));
end
figure,plot(1:16,tac_ideal,'s-',1:16,tac_rec(:,1),'^-',1:16,tac_rec(:,3),'o-',...
    1:16,tac_rec(:,4),'+-',1:16,tac_rec(:,6),'*-')
mse_tac=mean((repmat(tac_ideal,1,6)-tac_rec).^2);
figure,plot(1:16,tac_ideal/max(tac_ideal),'s-',1:16,tac_rec(:,3)/max(tac_ideal),'o-',...
    1:16,tac_rec(:,6)/max(tac_ideal),'*-')

% %phantom
% for n=1:16
% fn=['C:\simind\ncat16gles_act_' num2str(n) '.bin'];
% fid=fopen(fn,'rb');
% temp=fread(fid,'single');
% fclose(fid);
% ncat_phantom(:,:,:,n)=reshape(temp,64,64,64);
% end
% dsp(temp*5e5/sum(temp(:)))
% %analytic noiseless
% load ncatlestew_Im_idealA
% Im_idealA=Im_ideal;
% %air SIMIND noiseless: umass blur
% load ncatlestew_Im_idealUm
% %air SIMIND noiseless: 1.3cm FWHM blur
% Im_idealUm=Im_ideal;
% load ncatlestew_Im_ideal
% 
% for n=1:5
% Im_ideal_roiA=Im_idealA(25:49,16:40,34:38,:);
% Im_ideal_roiA=squeeze(Im_ideal_roiA(:,:,n,:));
% Im_ideal_roiA=permute(Im_ideal_roiA,[2 1 3]);
% dsp(Im_ideal_roiA,1)
% 
% Im_ideal_roi=Im_ideal(25:49,16:40,34:38,:);
% Im_ideal_roi=squeeze(Im_ideal_roi(:,:,n,:));
% Im_ideal_roi=permute(Im_ideal_roi,[2 1 3]);
% dsp(Im_ideal_roi,1)
% end
% %SIMIND has more blur than analytic reconstruction. (may need to increase
% %FWHM for Gaussian kenerl. Use Umass kenerl!!!)