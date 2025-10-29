%SNR recalculation for images without SC
%all images scale to the same counts as the ideal image. The SNRs in
%Im_*.mat are not right!!! Use SNR results in this routine! First set!
%Mingwu Jin, Oct 10, 2006
load ncat_simind_16g_nsino
sinoT=sinoP+sinoS;
%clear sinoP sionS
load Im_ideal
count_g1=sum(sum(sum(sinoT(:,:,:,1))));clear sinoP sinoS
ncatg1roi=Im_ideal(25:49,16:40,34:38)*count_g1/sum(Im_ideal(:));

for n=1:7%spatial, no AC
    filename=['Im_mapsnASC' num2str(n) '.mat'];
    load(filename,'Im_maps');
    Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1/sum(sum(sum(Im_maps(:,:,:,1))));
    snr_s(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
end
for n=1:6%spatial, AC
    filename=['Im_maps_AC' num2str(n) '.mat'];
    load(filename,'Im_maps');
    Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1/sum(sum(sum(Im_maps(:,:,:,1))));
    snr_sac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
end
for n=1:6%temporal, no AC    
    filename=['Im_idealMCnAC' num2str(n) '.mat'];
    load(filename,'Im_maps');
    Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1/sum(sum(sum(Im_maps(:,:,:,1))));
    snr_t(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
end
for n=1:6%temporal, AC    
    filename=['Im_idealMC_AC' num2str(n) '.mat'];
    load(filename,'Im_maps');
    Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1/sum(sum(sum(Im_maps(:,:,:,1))));
    snr_tac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
end

load Im_maps_ASC6 snr_sasc
load Im_idealMC_ASC6 snr_tasc

gbeta_t=[0.01:0.01:0.06];
gbeta_tac=[0.001:0.001:0.006];
gbeta_tasc=[0.001:0.001:0.006];
sbeta_s=[0.002:0.002:0.008 0.01:0.01:0.03];
sbeta_sac=[0.0005 0.001:0.001:0.005];
sbeta_sasc=[0.001:0.001:0.006];
save snr_results1 snr_s sbeta_s snr_sac sbeta_sac snr_sasc sbeta_sasc...
    snr_t gbeta_t snr_tac gbeta_tac snr_tasc gbeta_tasc