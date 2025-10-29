%CHO
%read frame 8 
loc=[-8 10];frame=8;gate=6;
load im_max4cho
%gate 6
% for g=6
%     filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
%     fid=fopen(filename,'rb');
%     temp=fread(fid,64^3,'single');
%     fclose(fid);
% end
% temp=reshape(temp,64,64,64);
% normal_g6=temp(23:52,16:43,29:48);
% temp=normal_g6(:,:,8)>max(normal_g6(:))*.5;
% ind_myo=find(temp);

% load ncatlestew_Im_idealUm
% Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% lesion_g1=Im_ideal(23:52,16:43,29:48,gate);
% immax_l=max(max(lesion_g1(:,:,frame)));
% load ncat16gtew_Im_idealUm
% Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% normal_g1=Im_ideal(23:52,16:43,29:48,gate);
% immax_n=max(max(normal_g1(:,:,frame)));%uniformly normalization is not good

%ind_myo=find(normal_g6.*repmat(temp,[1 1 20]));
% loc=[-14 5];frame=7;
% loc=loc+[1 2];
nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
les_dir_name='D:\imagereconstruction\4D\data_les_5st';
t_ind=[0 2 3];
Fi=Qfreq(64);
S=5;T=3;N=30;
int_noise=(32.4/256)^2/12;
Az_nAs=zeros(S,T);
ch_nAs_les=zeros(4,N,S,T);ch_nAs_nor=zeros(4,N,S,T);
immax_mapnAs_n=mean(immax_mapnAs_n,3);
immax_mapnAs_l=mean(immax_mapnAs_l,3);
for s=1:S
    for t=1:T        
        for n=1:N
            filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_nAS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax_mapnAs_l(s,t,n)=max(temp(ind_myo));
            temp=temp*immax_mapnAs_n(s,t)/immax_mapnAs_l(s,t);
            ch_nAs_les(:,n,s,t)=cho_feature(Fi,temp,loc);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_nAS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax_mapnAs_n(s,t,n)=max(temp(ind_myo));
%             temp=temp/immax_n;
            ch_nAs_nor(:,n,s,t)=cho_feature(Fi,temp,loc);
        end
        Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(:,:,s,t))',squeeze(ch_nAs_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_nAS ch_nAs_les ch_nAs_nor Az_nAs
% save cho_result_nAS8_g6 ch_nAs_les ch_nAs_nor Az_nAs
%internal noise~=(30/256)^2/12=0.0011
%%%%%%%%%%%%%%%
%Attenuation+scatter correction
Az_As=zeros(S,T);
ch_As_les=zeros(4,N,S,T);ch_As_nor=zeros(4,N,S,T);
immax_mapAs_n=mean(immax_mapAs_n,3);
immax_mapAs_l=mean(immax_mapAs_l,3);
for s=1:S
    for t=1:T        
        for n=1:N
            filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_AS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax_mapAs_l(s,t,n)=max(temp(ind_myo));
            temp=temp*immax_mapAs_n(s,t)/immax_mapAs_l(s,t);
            ch_As_les(:,n,s,t)=cho_feature(Fi,temp,loc);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_AS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax_mapAs_n(s,t,n)=max(temp(ind_myo));
%             temp=temp/immax_n;
            ch_As_nor(:,n,s,t)=cho_feature(Fi,temp,loc);
        end
        Az_As(s,t)=Azfinder(squeeze(ch_As_nor(:,:,s,t))',squeeze(ch_As_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_AS ch_As_les ch_As_nor Az_As
% save cho_result_AS8_g6 ch_As_les ch_As_nor Az_As
%%%%%%%%%%%%%%%
%Attenuation correction
Az_A=zeros(S,T);
ch_A_les=zeros(4,N,S,T);ch_A_nor=zeros(4,N,S,T);
immax_mapA_n=mean(immax_mapA_n,3);
immax_mapA_l=mean(immax_mapA_l,3);
for s=1:S
    for t=1:T        
        for n=1:N
            filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_A_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax_mapA_l(s,t,n)=max(temp(ind_myo));
            temp=temp*immax_mapA_n(s,t)/immax_mapA_l(s,t);
            ch_A_les(:,n,s,t)=cho_feature(Fi,temp,loc);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_A_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax_mapA_n(s,t,n)=max(temp(ind_myo));
%             temp=temp/immax_n;
            ch_A_nor(:,n,s,t)=cho_feature(Fi,temp,loc);
        end
        Az_A(s,t)=Azfinder(squeeze(ch_A_nor(:,:,s,t))',squeeze(ch_A_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_AS ch_As_les ch_As_nor Az_As
% save cho_result_A8_g6 ch_A_les ch_A_nor Az_A
%%%%%%%%%%%%%%%%%%%%%%%
%ST-121

% loc=[-8 10];frame=8; % [19 8]
% cen_gate=6;l_gate=5;r_gate=7;
% Fi=Qfreq(64);
% int_noise=(32.4/256)^2/12;N=30;
% nor_dir_name='D:\imagereconstruction\4D\data_16g_fbp';
% les_dir_name='D:\imagereconstruction\4D\data_les_fbp';
% ch_st121_les=zeros(4,N);ch_st121_nor=zeros(4,N);
% for n=1:N
%     filename=[les_dir_name '\Im_fbp_n'...
%         num2str(n) '.mat'];
%     load(filename,'Im_fbp');
%     temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);
% %     immax_l(n)=max(temp(ind_myo));
%     temp=temp*mean(immax_n)/mean(immax_l);
%     ch_st121_les(:,n)=cho_feature(Fi,temp,loc);
%     filename=[nor_dir_name '\Im_fbp_n'...
%         num2str(n) '.mat'];
%     load(filename,'Im_fbp');
%     temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);
% %     immax_n(n)=max(temp(ind_myo));
% %     temp=temp/immax_n(n);
%     ch_st121_nor(:,n)=cho_feature(Fi,temp,loc);
% end
% Az_st121=Azfinder(ch_st121_nor(:,:)',ch_st121_les(:,:)',int_noise);
% save im_max4cho immax*