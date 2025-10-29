function ncatlestew_isbi07_re_s(s)

%gbeta=[0.005:0.001:0.01 0.02:0.01:0.05];
sub_num=16;it_num=10;OF_tag=0;blur=1;
pause(s*10);
load ncatlestew_n1 sino_p sino_s
sinoP=sino_p;sinoS=sino_s;clear sino_p sino_s 
sinoT=sinoP+sinoS;clear sinoP
%clear sinoP sionS
load roi

if s==0%spatial only, no AC and SC:     
%     gbeta=0;sbeta=[0.007];
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,0,0,0);
%         filename=['ncatlestew_mapsnASCopt.mat'];
%         save(filename,'Im_maps');
%         Im_maps=mbsrem4dv2(sinoT,Im_maps,sub_num,it_num,OF_tag,sbeta,gbeta,blur,0,0,0);
%         filename=['ncatlestew_mapsnASCopt20.mat'];
%         save(filename,'Im_maps');
% elseif s==1%spatial only, AC, no SC: 10 For motion estimation
    gbeta=0;sbeta=0.0008;
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,0,0);
        filename=['ncatlestew_maps_ACopt.mat'];
%         save(filename,'Im_maps');
        load(filename,'Im_maps');
        Im_maps=mbsrem4dv2(sinoT,Im_maps,sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,0,0);
        filename=['ncatlestew_maps_ACopt20.mat'];
        save(filename,'Im_maps');
elseif s==1%spatial only, AC and SC:    
    gbeta=0;sbeta=0.001;
%         Im_maps=mbsrem4dv2(sinoT,repmat(roi,[1,1,64,16]),sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,sinoS,0);
        filename=['ncatlestew_maps_ASCopt.mat'];
%         save(filename,'Im_maps');
        load(filename,'Im_maps');
        Im_maps=mbsrem4dv2(sinoT,Im_maps,sub_num,it_num,OF_tag,sbeta,gbeta,blur,1,sinoS,0);
        filename=['ncatlestew_maps_ASCopt20.mat'];
        save(filename,'Im_maps');
end