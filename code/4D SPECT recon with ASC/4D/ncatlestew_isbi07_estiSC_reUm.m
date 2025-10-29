function ncatlestew_isbi07_estiSC_re(s)

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
load ncatlestew_esti_scat tew_scat

% if s==0%MC temporal only, AC and SC:
%     gbeta=[0.001:0.0002:0.0028];sbeta=0;
%     for n=1:10
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta(n),blur,1,tew_scat,0);
%         Im_maps=Im_maps*8e6/sum(Im_maps(:));
%         Im_maps=Im_maps(25:49,16:40,34:38,:);
%         Imroi=Im_maps(:,:,:,1);
%         snr_tascE(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
%         filename=['ncatlestewE_iMC_ASC_Um' num2str(n) '.mat'];
%         save(filename,'Im_maps','snr_tascE');
%     end
% elseif s==1%spatial only, AC and SC:
    gbeta=0;sbeta=[0.0021:0.0002:0.004];%sbeta=[0.0001:0.0002:0.002];
    for n=1:10
        Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta(n),gbeta,blur,1,tew_scat,0);
        Im_maps=Im_maps*8e6/sum(Im_maps(:));
        Im_maps=Im_maps(25:49,16:40,34:38,:);
        Imroi=Im_maps(:,:,:,1);
        snr_sascE(n)=10*log10(ncatg1roi(:)'*ncatg1roi(:)/sum((ncatg1roi(:)-Imroi(:)).^2));
        filename=['ncatlestewE_maps_ASC_Um' num2str(n+10) '.mat'];
        save(filename,'Im_maps','snr_sascE');
    end
% end