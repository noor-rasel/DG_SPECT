function idealMC(s)

%gbeta=[0.005:0.001:0.01 0.02:0.01:0.05];
sub_num=16;it_num=10;OF_tag=0;blur=1;
pause(s*10);
load ncat_simind_16g_nsino
sinoT=sinoP+sinoS;
%clear sinoP sionS
load Im_ideal
count_g1P=sum(sum(sum(sinoP(:,:,:,1))));clear sinoP
ncatg1roi=Im_ideal(25:49,16:40,34:38)*count_g1P/sum(Im_ideal(:));
count_g1=sum(sum(sum(sinoT(:,:,:,1))));
clear Im_ideal
load roi
if s==0%MC temporal only, no AC and SC: 0.04
%     gbeta=[0.01:0.01:0.06];sbeta=0;
%     for n=1:6
%         Im_maps=mbsrem4d(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,0,0);
%         Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1P/sum(sum(sum(Im_maps(:,:,:,1))));
%         snr_t(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
%         filename=['Im_idealMCnAC' num2str(n) '.mat'];
%         save(filename,'Im_maps','snr_t');
%     end
% elseif s==1%MC temporal only, AC, no SC: 0.005    
%     gbeta=[0.001:0.001:0.006];sbeta=0;
%     for n=1:6
%         Im_maps=mbsrem4d(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,1,0);
%         Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1P/sum(sum(sum(Im_maps(:,:,:,1))));
%         snr_tac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
%         filename=['Im_idealMC_AC' num2str(n) '.mat'];
%         save(filename,'Im_maps','snr_tac');
%     end
% elseif s==2%MC temporal only, AC and SC: ???    
    gbeta=[0.007:0.001:0.01 0.02 0.03];sbeta=0;
    for n=1:6
        Im_maps=mbsrem4d(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,1,sinoS);
        Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1P/sum(sum(sum(Im_maps(:,:,:,1))));
        snr_tasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['Im_idealMC_ASC' num2str(n+6) '.mat'];
        save(filename,'Im_maps','snr_tasc');
    end
elseif s==1%s==3%spatial only, no AC and SC: 0.01    
    gbeta=0;%sbeta=[0.002:0.002:0.008 0.01:0.01:0.03];
    sbeta=[0.04:0.01:0.1];
    for n=1:7
        Im_maps=mbsrem4d(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,0,0);
        Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1P/sum(sum(sum(Im_maps(:,:,:,1))));
        snr_s(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['Im_mapsnASC' num2str(n+7) '.mat'];
        save(filename,'Im_maps','snr_s');
    end
% elseif s==4%spatial only, AC, no SC: 0.002
%     gbeta=0;sbeta=[0.0005 0.001:0.001:0.005];
%     for n=1:6
%         Im_maps=mbsrem4d(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,0);
%         Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1P/sum(sum(sum(Im_maps(:,:,:,1))));
%         snr_sac(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
%         filename=['Im_maps_AC' num2str(n) '.mat'];
%         save(filename,'Im_maps','snr_sac');
%     end
% elseif s==5%spatial only, AC and SC: 0.003    
%     gbeta=0;sbeta=[0.001:0.001:0.006];
%     for n=1:6
%         Im_maps=mbsrem4d(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,sinoS);
%         Imroi=Im_maps(25:49,16:40,34:38,1)*count_g1P/sum(sum(sum(Im_maps(:,:,:,1))));
%         snr_sasc(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
%         filename=['Im_maps_ASC' num2str(n) '.mat'];
%         save(filename,'Im_maps','snr_sasc');
%     end
end
