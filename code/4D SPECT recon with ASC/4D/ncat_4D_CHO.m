%CHO
%read frame 8 (less obvious) [46,27] or 9 (obvious) [46,27]
loc=[-14 5];frame=8;%[47,28]
loc=[-14 6];
% loc=[-11 4];%loc=[-11 4]
% loc=[-13 7];frame=9;%[45,26]
load ncat_g1myo_index ind_myo
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
%ind_myo=find(normal_g6.*repmat(temp,[1 1 20]));
% loc=[-14 5];frame=7;
% loc=loc+[1 2];
nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
les_dir_name='D:\imagereconstruction\4D\data_les_5st';
t_ind=[0 2 3];
Fi=Qfreq(64);
S=5;T=3;N=30;gate=1;
int_noise=0.0011;
Az_nAs=zeros(S,T);
ch_nAs_les=zeros(4,N,S,T);ch_nAs_nor=zeros(4,N,S,T);
for s=1:S
    for t=1:T        
        for n=1:N
            filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_nAS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax=max(temp(ind_myo));
%             temp=temp(:,:,frame)/immax;
            ch_nAs_les(:,n,s,t)=cho_feature(Fi,temp,loc);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_nAS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax=max(temp(ind_myo));
%             temp=temp(:,:,frame)/immax;
            ch_nAs_nor(:,n,s,t)=cho_feature(Fi,temp,loc);
        end
        Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(:,:,s,t))',squeeze(ch_nAs_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_nAS ch_nAs_les ch_nAs_nor Az_nAs
%save cho_nAS_64 ch_nAs_les ch_nAs_nor Az_nAs
save cho_nAS_64_center ch_nAs_les ch_nAs_nor Az_nAs
%internal noise~=(30/256)^2/12=0.0011
%%%%%%%%%%%%%%%
%Attenuation+scatter correction
Az_As=zeros(S,T);
ch_As_les=zeros(4,N,S,T);ch_As_nor=zeros(4,N,S,T);
for s=1:S
    for t=1:T        
        for n=1:N
            filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_AS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax=max(temp(ind_myo));
%             temp=temp(:,:,frame)/immax;
            ch_As_les(:,n,s,t)=cho_feature(Fi,temp,loc);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_AS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax=max(temp(ind_myo));
%             temp=temp(:,:,frame)/immax;
            ch_As_nor(:,n,s,t)=cho_feature(Fi,temp,loc);
        end
        Az_As(s,t)=Azfinder(squeeze(ch_As_nor(:,:,s,t))',squeeze(ch_As_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_AS ch_As_les ch_As_nor Az_As
%save cho_AS_64 ch_As_les ch_As_nor Az_As
save cho_AS_64_center ch_As_les ch_As_nor Az_As
%%%%%%%%%%%%%%%
%Attenuation correction
Az_A=zeros(S,T);
ch_A_les=zeros(4,N,S,T);ch_A_nor=zeros(4,N,S,T);
for s=1:S
    for t=1:T        
        for n=1:N
            filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_A_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax=max(temp(ind_myo));
%             temp=temp(:,:,frame)/immax;
            ch_A_les(:,n,s,t)=cho_feature(Fi,temp,loc);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_A_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
%             immax=max(temp(ind_myo));
%             temp=temp(:,:,frame)/immax;
            ch_A_nor(:,n,s,t)=cho_feature(Fi,temp,loc);
        end
        Az_A(s,t)=Azfinder(squeeze(ch_A_nor(:,:,s,t))',squeeze(ch_A_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_AS ch_As_les ch_As_nor Az_As
%save cho_A_64 ch_A_les ch_A_nor Az_A
save cho_A_64_center ch_A_les ch_A_nor Az_A
%%%%%%%%%%%%%%%%%%%%%%%
%ST-121
% loc=[-14 5];frame=8;
cen_gate=1;l_gate=16;r_gate=2;
Fi=Qfreq(64);
int_noise=0.0011;;N=30;
nor_dir_name='D:\imagereconstruction\4D\data_16g_fbp';
les_dir_name='D:\imagereconstruction\4D\data_les_fbp';
ch_st121_les=zeros(4,N);ch_st121_nor=zeros(4,N);
for n=1:N
    filename=[les_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,:,cen_gate)+.25*Im_fbp(:,:,:,l_gate)+.25*Im_fbp(:,:,:,r_gate);
    temp=temp(:,:,frame);
    ch_st121_les(:,n)=cho_feature(Fi,temp,loc);
    temp=Im_fbp(:,:,:,cen_gate);
    temp=temp(:,:,frame);
    ch_fbp_les(:,n)=cho_feature(Fi,temp,loc);
    filename=[nor_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,:,cen_gate)+.25*Im_fbp(:,:,:,l_gate)+.25*Im_fbp(:,:,:,r_gate);
    temp=temp(:,:,frame)/1.1189;
    ch_st121_nor(:,n)=cho_feature(Fi,temp,loc);
    temp=Im_fbp(:,:,:,cen_gate);
    temp=temp(:,:,frame)/1.1189;
    ch_fbp_nor(:,n)=cho_feature(Fi,temp,loc);
end
Az_st121=Azfinder(ch_st121_nor(:,:)',ch_st121_les(:,:)',int_noise);
Az_fbp=Azfinder(ch_fbp_nor(:,:)',ch_fbp_les(:,:)',int_noise);
% save cho_st121_64 ch_st121_les ch_st121_nor Az_st121
% save cho_fbp_64 ch_fbp_les ch_fbp_nor Az_fbp
save cho_st121_64_center ch_st121_les ch_st121_nor Az_st121
save cho_fbp_64_center ch_fbp_les ch_fbp_nor Az_fbp
% %0.9950 for ST121; 0.9903 for FBP : gate 1 frame 7
% %gate 1 frame 8
% %loc=[-14 5]: Az=0.8996 (cho_result_A8 better)
% %loc=[-13 7]: Az=0.9499
% %loc=[-12 4]; Az=0.9336
% %loc=[-11 4]; Az=0.9044 (cho_result_A8_3 worse)
% %loc=[-13 8]; Az=0.8865##########(cho_result_A8_2 best)
% %loc=[-12 8]; Az=0.9398


%plot
% load cho_fbp_64
% load cho_st121_64
% load cho_nAS_64
% load cho_A_64
% load cho_AS_64
% 
int_noise=0;
for s=1:S
    for t=1:T
        Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(:,:,s,t))',squeeze(ch_nAs_les(:,:,s,t))',int_noise);
    end
end
for s=1:S
    for t=1:T 
        Az_As(s,t)=Azfinder(squeeze(ch_As_nor(:,:,s,t))',squeeze(ch_As_les(:,:,s,t))',int_noise);
    end
end
for s=1:S
    for t=1:T        
        Az_A(s,t)=Azfinder(squeeze(ch_A_nor(:,:,s,t))',squeeze(ch_A_les(:,:,s,t))',int_noise);
    end
end
Az_st121=Azfinder(ch_st121_nor(:,:)',ch_st121_les(:,:)',int_noise);
Az_fbp=Azfinder(ch_fbp_nor(:,:)',ch_fbp_les(:,:)',int_noise);

az_fbp=Az_fbp;
az_st121=Az_st121;%0.5*(1+erf(snr_st121/2));
az_mapnAS=Az_nAs';
az_mapA=Az_A';%0.5*(1+erf(snr_mapA/2));
az_mapAS=Az_As';%0.5*(1+erf(snr_mapAS/2));

sbeta=[0 5e-4 1e-3 1e-2 3e-2];
figure,plot(sbeta,az_mapnAS(1,:),'*-',sbeta,az_mapnAS(2,:),'+-',sbeta,az_mapnAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','FBP','ST 121')
title('NC'),xlabel('\beta_s'),ylabel('SNR')

sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,az_mapA(1,:),'*-',sbeta,az_mapA(2,:),'+-',sbeta,az_mapA(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','FBP','ST 121')
title('AC'),xlabel('\beta_s'),ylabel('SNR')

sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,az_mapAS(1,:),'*-',sbeta,az_mapAS(2,:),'+-',sbeta,az_mapAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','FBP','ST 121')
title('ASC'),xlabel('\beta_s'),ylabel('SNR')
%detail
sbeta=[0 5e-4 1e-3 1e-2 3e-2];
figure,plot(sbeta,az_mapnAS(2,:),'+-',sbeta,az_mapnAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','ST 121')
title('NC'),xlabel('\beta_s'),ylabel('SNR')

sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,az_mapA(2,:),'+-',sbeta,az_mapA(3,:),'o-')
hold on,plot(sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','ST 121')
title('AC'),xlabel('\beta_s'),ylabel('SNR')

sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,az_mapAS(2,:),'+-',sbeta,az_mapAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','ST 121')
title('ASC'),xlabel('\beta_s'),ylabel('SNR')


%isbi07
% loc=[-14 5];frame=7;
% loc=loc+[1 2];%[1 2] (24(y) 11(x)) [1 1] (24(y) 12(x))
% nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
% les_dir_name='D:\imagereconstruction\4D\data_les_5st';
% st_comb=[4 0;1 2];
% Fi=Qfreq(64);
% S=2;N=30;
% int_noise=0.0011;
% Az_nAs=zeros(S,1);
% ch_nAs_les=zeros(4,N,S);ch_nAs_nor=zeros(4,N,S);
% for s=1:S
%     for n=1:N
%         filename=[les_dir_name num2str(st_comb(s,2)) '\Im_ncatles_nAS_s'...
%             num2str(st_comb(s,1)) 't' num2str(st_comb(s,2)) '_n' num2str(n) '.mat'];
%         load(filename,'Im_maps');
%         temp=Im_maps(:,:,frame,1);
%         ch_nAs_les(:,n,s)=cho_feature(Fi,temp,loc);
%         filename=[nor_dir_name num2str(st_comb(s,2)) '\Im_ncat16g_nAS_s'...
%             num2str(st_comb(s,1)) 't' num2str(st_comb(s,2)) '_n' num2str(n) '.mat'];
%         load(filename,'Im_maps');
%         temp=Im_maps(:,:,frame,1);
%         ch_nAs_nor(:,n,s)=cho_feature(Fi,temp,loc);
%     end
%     Az_nAs(s)=Azfinder(squeeze(ch_nAs_nor(:,:,s))',squeeze(ch_nAs_les(:,:,s))',int_noise);
% end
% 
% Az_AS=zeros(S,1);
% ch_AS_les=zeros(4,N,S);ch_AS_nor=zeros(4,N,S);
% for s=1:S
%     for n=1:N
%         filename=[les_dir_name num2str(st_comb(s,2)) '\Im_ncatles_A_s'...
%             num2str(st_comb(s,1)) 't' num2str(st_comb(s,2)) '_n' num2str(n) '.mat'];
%         load(filename,'Im_maps');
%         temp=Im_maps(:,:,frame,1);
%         ch_AS_les(:,n,s)=cho_feature(Fi,temp,loc);
%         filename=[nor_dir_name num2str(st_comb(s,2)) '\Im_ncat16g_A_s'...
%             num2str(st_comb(s,1)) 't' num2str(st_comb(s,2)) '_n' num2str(n) '.mat'];
%         load(filename,'Im_maps');
%         temp=Im_maps(:,:,frame,1);
%         ch_AS_nor(:,n,s)=cho_feature(Fi,temp,loc);
%     end
%     Az_AS(s)=Azfinder(squeeze(ch_AS_nor(:,:,s))',squeeze(ch_AS_les(:,:,s))',int_noise);
% end