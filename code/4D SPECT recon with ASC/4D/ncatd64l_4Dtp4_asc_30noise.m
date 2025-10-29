function ncatd64l_4Dtp4_asc_30noise(com_numb)

%gbeta=[0.005:0.001:0.01 0.02:0.01:0.05];
sub_num=16;it_num=10;OF_tag=0;blur=1;
[s,k]=ind2sub([2,30],com_numb+1);
s=s+1;
%get noise sinogram
filename=['ncatd64l_n' num2str(k) '.mat'];
load(filename)
%get sinoT and tew_scat
load ncatd64l_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
clear Im_ideal
load roi


if s==1%spatial only, no AC and SC:
    gbeta=[5e-3 1e-2 2e-2];%[0.001:0.002:0.02];%0.01
    sbeta=[0 5e-4 1e-3 1e-2 3e-2];
    for tt=3%1:3
        for n=1:5
            Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta(tt),blur,0,0,0,k-1);
            Im_maps=Im_maps*8e6/sum(Im_maps(:));
            Im_maps=Im_maps(23:52,16:43,29:48,:);
            Imroi=Im_maps(:,:,:,1);
            snr_t(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
            filename=['./data_d64l_5st3/Im_ncatd64l_nAS_s' num2str(n) 't' num2str(tt) '_n' num2str(k) '.mat'];
            save(filename,'Im_maps','snr_t');
        end
    end
elseif s==2%spatial only, AC, no SC:
    gbeta=[5e-4 1e-3 2e-3 3e-3];%[0.0001:0.0002:0.002];%0.001
    sbeta=[0 1e-4 1e-3 5e-3 9e-3];%[0 5e-5 1e-4 1e-3 3e-3];
    for tt=4%1:3
        for n=1:5
            Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta(tt),blur,1,0,0,k-1);
            Im_maps=Im_maps*8e6/sum(Im_maps(:));
            Im_maps=Im_maps(23:52,16:43,29:48,:);
            Imroi=Im_maps(:,:,:,1);
            snr_tac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
            filename=['./data_d64l_5st4/Im_ncatd64l_A_s' num2str(n) 't' num2str(tt) '_n' num2str(k) '.mat'];
            save(filename,'Im_maps','snr_tac');
        end
    end
elseif s==3%spatial only, AC and SC:
    %TEW method
    sinoL=reshape(sinoL,[64^2 64 16])*4;
    tew_scat=buttlpf(sinoL,0.2*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    gbeta=[1e-3 2e-3 3e-3 4e-3];%[0.001:0.0002:0.0028];%0.002
    sbeta=[0 5e-4 2e-3 6e-3 1e-2];%[0 5e-4 1e-3 2e-3 4e-3];%0.002
    for tt=4%1:3
        for n=1:5
            Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta(tt),blur,1,tew_scat,0,k-1);
            Im_maps=Im_maps*8e6/sum(Im_maps(:));
            Im_maps=Im_maps(23:52,16:43,29:48,:);
            Imroi=Im_maps(:,:,:,1);
            snr_tasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
            filename=['./data_d64l_5st4/Im_ncatd64l_AS_s' num2str(n) 't' num2str(tt) '_n' num2str(k) '.mat'];
            save(filename,'Im_maps','snr_tasc');
        end
    end
end