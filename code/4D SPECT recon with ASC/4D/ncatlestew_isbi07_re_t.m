function ncatlestew_isbi07_re_t(s)

%gbeta=[0.005:0.001:0.01 0.02:0.01:0.05];
sub_num=16;it_num=10;OF_tag=0;blur=1;
pause(s*10);
load ncatlestew_n1 sino_p sino_s
sinoP=sino_p;sinoS=sino_s;clear sino_p sino_s 
sinoT=sinoP+sinoS;clear sinoP
%clear sinoP sionS
load roi

if s==0%temporal only, no AC and SC:     
%     gbeta=0.01;sbeta=0;
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,0,0,0);
%         filename=['ncatlestew_eMCnASCopt.mat'];
%         save(filename,'Im_maps');
%         Im_maps=mbsrem4dv2(sinoT,Im_maps,sub_num,it_num,OF_tag,sbeta,gbeta,blur,0,0,0);
%         filename=['ncatlestew_eMCnASCopt20.mat'];
%         save(filename,'Im_maps');
% elseif s==1%temporal only, AC, no SC: 10 For motion estimation
    gbeta=0.001;sbeta=0;
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,0,0);
        filename=['ncatlestew_eMC_ACopt.mat'];
%         save(filename,'Im_maps');
        load(filename,'Im_maps');
        Im_maps=mbsrem4dv2(sinoT,Im_maps,sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,0,0);
        filename=['ncatlestew_eMC_ACopt20.mat'];
        save(filename,'Im_maps');
elseif s==1%temporal only, AC and SC:    
    gbeta=0.002;sbeta=0;
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,sinoS,0);
        filename=['ncatlestew_eMC_ASCopt.mat'];
%         save(filename,'Im_maps');
        load(filename,'Im_maps');
        Im_maps=mbsrem4dv2(sinoT,Im_maps,sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,sinoS,0);
        filename=['ncatlestew_eMC_ASCopt20.mat'];
        save(filename,'Im_maps');
end