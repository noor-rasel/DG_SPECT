function ncat16g_4Dtp_optSC(k)


sub_num=16;it_num=10;OF_tag=0;blur=1;
pause(mod(k,10)*60);
%get noise sinogram
% tic,

%get sinoT and tew_scat
load ncat16gtew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
clear Im_ideal
load roi

%TEW method
gbeta=[1e-3 2e-3 3e-3];%[0.001:0.0002:0.0028];%0.002
sbeta=[0 5e-4 1e-3 2e-3 4e-3];%0.002
filename=['ncat16g_n' num2str(k+1) '.mat'];
load(filename)
sinoL=reshape(sinoL,[64^2 64 16])*4;
cf=0.1:0.1:1;
for com_numb=1:10
    tew_scat=buttlpf(sinoL,cf(com_numb)*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(1),gbeta(2),blur,1,tew_scat,0,k-1);
    Im_maps=Im_maps*8e6/sum(Im_maps(:));
    Im_maps=Im_maps(23:52,16:43,29:48,:);
    Imroi=Im_maps(:,:,:,1);
    snr_tasc(com_numb)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
    filename=['./optSC/Im_optSC_CF' num2str(com_numb) '_n' num2str(k+1) '.mat'];
    save(filename,'Im_maps','snr_tasc');
end
%conclusion: cf=0.1 is best in terms of SNR!