function ncatd64_4Dsp_asc_30noise(com_numb)

%gbeta=[0.005:0.001:0.01 0.02:0.01:0.05];
sub_num=16;it_num=10;OF_tag=0;blur=1;
[s,k]=ind2sub([3,30],com_numb+1);
filename=['ncatd64_n' num2str(k) '.mat'];
load(filename)
%get sinoT and tew_scat
load ncatd64_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
clear Im_ideal
load roi

if s==1%spatial only, no AC and SC:     
    gbeta=0;sbeta=[0 5e-4 1e-3 1e-2 3e-2];%[0.001:0.002:0.02];%0.007
    for n=1:5
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,0,0,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(23:52,16:43,29:48,:);
        Imroi=Im_maps(:,:,:,1);
        snr_s(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['./data_d64_5st0/Im_d64_nAS_s' num2str(n) 't0_n' num2str(k) '.mat'];
        save(filename,'Im_maps','snr_s');
    end
elseif s==2%spatial only, AC, no SC: 
    gbeta=0;sbeta=[0 1e-4 1e-3 5e-3 9e-3];%[0.0001:0.0002:0.002];%0.0008
    for n=1:5
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,0,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(23:52,16:43,29:48,:);
        Imroi=Im_maps(:,:,:,1);
        snr_sac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['./data_d64_5st0/Im_d64_A_s' num2str(n) 't0_n' num2str(k) '.mat'];
        save(filename,'Im_maps','snr_sac');
    end
elseif s==3%spatial only, AC and SC:
    %TEW method
    sinoL=reshape(sinoL,[64^2 64 16])*4;
    tew_scat=buttlpf(sinoL,0.2*.634,3);
%     tew_scat=buttlpf(sinoL,0.4*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    gbeta=0;sbeta=[0 5e-4 2e-3 6e-3 1e-2];%[0.0001:0.0002:0.002];%0.001
    for n=1:5
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,tew_scat,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(23:52,16:43,29:48,:);
        Imroi=Im_maps(:,:,:,1);
        snr_sasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['./data_d64_5st0/Im_d64_AS.2LPF_s' num2str(n) 't0_n' num2str(k) '.mat'];
%         filename=['./data_d64l_5st0/Im_d64l_AS_s' num2str(n) 't0_n' num2str(k) '.mat'];
        save(filename,'Im_maps','snr_sasc');
    end
end