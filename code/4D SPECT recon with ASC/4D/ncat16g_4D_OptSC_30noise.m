function ncat16g_4D_OptSC_30noise(com_numb)

sub_num=16;it_num=10;OF_tag=0;blur=1;
[s,k]=ind2sub([4,30],com_numb+1);
pause(mod(com_numb,24)*10);
filename=['ncat16g_n' num2str(k) '.mat'];
load(filename)
%get sinoT and tew_scat
load ncat16gtew_Im_idealUm

Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
ncatg1roi=Im_ideal(23:52,16:43,29:48,1);
clear Im_ideal
load roi

if s==1%spatial only, AC and SC:
    %TEW method
    sinoL=reshape(sinoL,[64^2 64 16])*4;
    tew_scat=buttlpf(sinoL,0.2*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    gbeta=0;sbeta=[0 5e-4 1e-3 2e-3 4e-3];%[0.0001:0.0002:0.002];%0.001
    for n=1:5
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,tew_scat,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(23:52,16:43,29:48,:);
        Imroi=Im_maps(:,:,:,1);
        snr_sasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['./data_16g_optSC/Im_ncat16g_AS_s' num2str(n) 't0_n' num2str(k) '.mat'];
        save(filename,'Im_maps','snr_sasc');
    end
elseif s==2%temporal 1, AC and SC:
    %TEW method
    sinoL=reshape(sinoL,[64^2 64 16])*4;
    tew_scat=buttlpf(sinoL,0.2*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    gbeta=[1e-3 2e-3 3e-3];%[0.001:0.0002:0.0028];%0.002
    sbeta=[0 5e-4 1e-3 2e-3 4e-3];%0.002
    for tt=1%1:3
        for n=1:5
            Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta(tt),blur,1,tew_scat,0,k-1);
            Im_maps=Im_maps*8e6/sum(Im_maps(:));
            Im_maps=Im_maps(23:52,16:43,29:48,:);
            Imroi=Im_maps(:,:,:,1);
            snr_tasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
            filename=['./data_16g_optSC/Im_ncat16g_AS_s' num2str(n) 't' num2str(tt) '_n' num2str(k) '.mat'];
            save(filename,'Im_maps','snr_tasc');
        end
    end
elseif s==3%temporal 2, AC and SC:
    %TEW method
    sinoL=reshape(sinoL,[64^2 64 16])*4;
    tew_scat=buttlpf(sinoL,0.2*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    gbeta=[1e-3 2e-3 3e-3];%[0.001:0.0002:0.0028];%0.002
    sbeta=[0 5e-4 1e-3 2e-3 4e-3];%0.002
    for tt=2%1:3
        for n=1:5
            Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta(tt),blur,1,tew_scat,0,k-1);
            Im_maps=Im_maps*8e6/sum(Im_maps(:));
            Im_maps=Im_maps(23:52,16:43,29:48,:);
            Imroi=Im_maps(:,:,:,1);
            snr_tasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
            filename=['./data_16g_optSC/Im_ncat16g_AS_s' num2str(n) 't' num2str(tt) '_n' num2str(k) '.mat'];
            save(filename,'Im_maps','snr_tasc');
        end
    end
elseif s==4%temporal 3, AC and SC:
    %TEW method
    sinoL=reshape(sinoL,[64^2 64 16])*4;
    tew_scat=buttlpf(sinoL,0.2*.634,3);
    tew_scat=reshape(tew_scat,[64 64 64 16]);
    gbeta=[1e-3 2e-3 3e-3];%[0.001:0.0002:0.0028];%0.002
    sbeta=[0 5e-4 1e-3 2e-3 4e-3];%0.002
    for tt=3%1:3
        for n=1:5
            Im_maps=mbsrem4dv3(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta(tt),blur,1,tew_scat,0,k-1);
            Im_maps=Im_maps*8e6/sum(Im_maps(:));
            Im_maps=Im_maps(23:52,16:43,29:48,:);
            Imroi=Im_maps(:,:,:,1);
            snr_tasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
            filename=['./data_16g_optSC/Im_ncat16g_AS_s' num2str(n) 't' num2str(tt) '_n' num2str(k) '.mat'];
            save(filename,'Im_maps','snr_tasc');
        end
    end
end