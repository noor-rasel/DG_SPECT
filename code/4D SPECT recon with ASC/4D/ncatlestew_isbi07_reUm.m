function ncatlestew_isbi07_reUm(s)

%gbeta=[0.005:0.001:0.01 0.02:0.01:0.05];
sub_num=16;it_num=10;OF_tag=0;blur=1;
pause(s*10);
load ncatlestew_n1 sino_p sino_s
sinoP=sino_p;sinoS=sino_s;clear sino_p sino_s 
sinoT=sinoP+sinoS;clear sinoP
%clear sinoP sionS
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
ncatg1roi=Im_ideal(25:49,16:40,34:38,1);
clear Im_ideal
load roi
if s==0%MC temporal only, no AC and SC: 
    gbeta=[0.001:0.002:0.02];sbeta=0;%0.01
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,0,0,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_t(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestew_iMCnASC_Um' num2str(n) '.mat'];
        save(filename,'Im_maps','snr_t');
    end
elseif s==1%MC temporal only, AC, no SC:     
    gbeta=[0.0001:0.0002:0.002];sbeta=0;%0.001
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,1,0,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_tac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestew_iMC_AC_Um' num2str(n) '.mat'];
        save(filename,'Im_maps','snr_tac');
    end
elseif s==2%MC temporal only, AC and SC:
    gbeta=[0.001:0.0002:0.0028];sbeta=0;%0.002
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,1,sinoS,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_tasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestew_iMC_ASC_Um' num2str(n) '.mat'];
        save(filename,'Im_maps','snr_tasc');
    end
elseif s==3%spatial only, no AC and SC:     
    gbeta=0;sbeta=[0.001:0.002:0.02];%0.007
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,0,0,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_s(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestew_mapsnASC_Um' num2str(n) '.mat'];
        save(filename,'Im_maps','snr_s');
    end
elseif s==4%spatial only, AC, no SC: 
    gbeta=0;sbeta=[0.0001:0.0002:0.002];%0.0008
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,0,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_sac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestew_maps_AC_Um' num2str(n) '.mat'];
        save(filename,'Im_maps','snr_sac');
    end
elseif s==5%spatial only, AC and SC:    
    gbeta=0;sbeta=[0.0001:0.0002:0.002];%0.001
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,sinoS,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_sasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestew_maps_ASC_Um' num2str(n) '.mat'];
        save(filename,'Im_maps','snr_sasc');
    end
end