frame=8;
% loc=[-33 18];%<=>[104 53](center [71 71])%no postfix
%loc=[-43 8];%[114 63] _loc64
% loc=[-41 14];%[112 57] _center
% loc=[-38 14];%[109 57] _centerlu
loc=[-39 14];%[110 57] _centerld
% loc=[-38 13];%[109 58] _centerru
% loc=[-39 13];%[110 58] _centerrd

gate=1;
iFi=Qfreq_s(140,loc);iFi=reshape(iFi,[140*140 4])';

nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
les_dir_name='D:\imagereconstruction\4D\data_les_5st';
t_ind=[0 2 3];
S=5;T=3;N=30;
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
            temp=postproc_cho(temp,140);
            ch_nAs_les(:,n,s,t)=iFi*temp(:);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_nAS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
            temp=postproc_cho(temp,140);
            ch_nAs_nor(:,n,s,t)=iFi*temp(:);
        end
        Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(:,:,s,t))',squeeze(ch_nAs_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_nAS140_g1f8 ch_nAs_les ch_nAs_nor Az_nAs
% save cho_result_nAS140_g1f8_loc64 ch_nAs_les ch_nAs_nor Az_nAs
save cho_result_nAS140_g1f8_centerld ch_nAs_les ch_nAs_nor Az_nAs
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
            temp=postproc_cho(temp,140);
%             avg_les(:,:,n)=temp;
            ch_As_les(:,n,s,t)=iFi*temp(:);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_AS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
            temp=postproc_cho(temp,140);
%             avg_nor(:,:,n)=temp;
            ch_As_nor(:,n,s,t)=iFi*temp(:);
        end
        Az_As(s,t)=Azfinder(squeeze(ch_As_nor(:,:,s,t))',squeeze(ch_As_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_AS140_g1f8 ch_As_les ch_As_nor Az_As
% save cho_result_AS140_g1f8_loc64 ch_As_les ch_As_nor Az_As
save cho_result_AS140_g1f8_centerld ch_As_les ch_As_nor Az_As
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
            temp=postproc_cho(temp,140);
            ch_A_les(:,n,s,t)=iFi*temp(:);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_A_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
            temp=postproc_cho(temp,140);
            ch_A_nor(:,n,s,t)=iFi*temp(:);
        end
        Az_A(s,t)=Azfinder(squeeze(ch_A_nor(:,:,s,t))',squeeze(ch_A_les(:,:,s,t))',int_noise);
    end
end
% save cho_result_A140_g1f8 ch_A_les ch_A_nor Az_A
% save cho_result_A140_g1f8_loc64 ch_A_les ch_A_nor Az_A
save cho_result_A140_g1f8_centerld ch_A_les ch_A_nor Az_A
%%%%%%%%%%%%%%%%%%%%%%%
%ST-121
cen_gate=1;l_gate=16;r_gate=2;
nor_dir_name='D:\imagereconstruction\4D\data_16g_fbp';
les_dir_name='D:\imagereconstruction\4D\data_les_fbp';
ch_st121_les=zeros(4,N);ch_st121_nor=zeros(4,N);
for n=1:N
    filename=[les_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);%Im_fbp(:,:,frame,cen_gate);%
    temp=postproc_cho(temp,140);
%     avg_les(:,:,n)=temp;
    ch_st121_les(:,n)=iFi*temp(:);
    filename=[nor_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);%Im_fbp(:,:,frame,cen_gate);%
    temp=temp/1.1189;
    temp=postproc_cho(temp,140);
% 	avg_nor(:,:,n)=temp;
    ch_st121_nor(:,n)=iFi*temp(:);
end
Az_st121=Azfinder(ch_st121_nor(:,:)',ch_st121_les(:,:)',int_noise);
% save cho_result_st121_140_g1f8 ch_st121_les ch_st121_nor Az_st121
% save cho_result_st121_140_g1f8_loc64 ch_st121_les ch_st121_nor Az_st121
save cho_result_st121_140_g1f8_centerld ch_st121_les ch_st121_nor Az_st121
cen_gate=1;
nor_dir_name='D:\imagereconstruction\4D\data_16g_fbp';
les_dir_name='D:\imagereconstruction\4D\data_les_fbp';
ch_fbp_les=zeros(4,N);ch_fbp_nor=zeros(4,N);
for n=1:N
    filename=[les_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=Im_fbp(:,:,frame,cen_gate);
    temp=postproc_cho(temp,140);
%     avg_les(:,:,n)=temp;
    ch_fbp_les(:,n)=iFi*temp(:);
    filename=[nor_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=Im_fbp(:,:,frame,cen_gate);
    temp=temp/1.1189;
    temp=postproc_cho(temp,140);
%     avg_nor(:,:,n)=temp;
    ch_fbp_nor(:,n)=iFi*temp(:);
end
Az_fbp=Azfinder(ch_fbp_nor(:,:)',ch_fbp_les(:,:)',int_noise);
% save cho_result_fbp_140_g1f8 ch_fbp_les ch_fbp_nor Az_fbp
% save cho_result_fbp_140_g1f8_loc64 ch_fbp_les ch_fbp_nor Az_fbp
save cho_result_fbp_140_g1f8_centerld ch_fbp_les ch_fbp_nor Az_fbp
%plot
% load cho_result_fbp_140_g1f8
% load cho_result_st121_140_g1f8
% load cho_result_nAS140_g1f8
% load cho_result_A140_g1f8
% load cho_result_AS140_g1f8
% 

% load cho_result_fbp_140_g1f8_center
% load cho_result_st121_140_g1f8_center
% load cho_result_nAS140_g1f8_center
% load cho_result_A140_g1f8_center
% load cho_result_AS140_g1f8_center
% 
% int_noise=0.00003;%0.0011
% for s=1:S
%     for t=1:T
%         Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(1:2,:,s,t))',squeeze(ch_nAs_les(1:2,:,s,t))',int_noise);
%     end
% end
% for s=1:S
%     for t=1:T 
%         Az_As(s,t)=Azfinder(squeeze(ch_As_nor(1:2,:,s,t))',squeeze(ch_As_les(1:2,:,s,t))',int_noise);
%     end
% end
% for s=1:S
%     for t=1:T        
%         Az_A(s,t)=Azfinder(squeeze(ch_A_nor(1:2,:,s,t))',squeeze(ch_A_les(1:2,:,s,t))',int_noise);
%     end
% end
% Az_st121=Azfinder(ch_st121_nor(1:2,:)',ch_st121_les(1:2,:)',int_noise);
% Az_fbp=Azfinder(ch_fbp_nor(1:2,:)',ch_fbp_les(1:2,:)',int_noise);

az_fbp=Az_fbp;
az_st121=Az_st121;%0.5*(1+erf(snr_st121/2));
az_mapnAS=Az_nAs';
az_mapA=Az_A';%0.5*(1+erf(snr_mapA/2));
az_mapAS=Az_As';%0.5*(1+erf(snr_mapAS/2));

sbeta=[0 5e-4 1e-3 1e-2 3e-2];
figure,plot(sbeta,az_mapnAS(1,:),'*-',sbeta,az_mapnAS(2,:),'+-',sbeta,az_mapnAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=1\times 10^{-2}','\beta_t=2\times 10^{-2}','FBP','ST 121')
title('NC'),xlabel('\beta_s'),ylabel('Az')

sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,az_mapA(1,:),'*-',sbeta,az_mapA(2,:),'+-',sbeta,az_mapA(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','FBP','ST 121')
title('AC'),xlabel('\beta_s'),ylabel('Az')

sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,az_mapAS(1,:),'*-',sbeta,az_mapAS(2,:),'+-',sbeta,az_mapAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_fbp,1,5),'k-',sbeta,repmat(az_st121,1,5),'k--')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','FBP','ST 121')
title('ASC'),xlabel('\beta_s'),ylabel('Az')

%channel output
% figure, plot(1:30,ch_nAs_nor(:,:,1,2))
% hold on, plot(1:30,ch_nAs_les(:,:,1,2),'--')
% 
% figure, plot(1:30,ch_As_nor(:,:,1,2))
% hold on, plot(1:30,ch_As_les(:,:,1,2),'--')
% 
% figure, plot(1:30,ch_A_nor(:,:,1,2))
% hold on, plot(1:30,ch_A_les(:,:,1,2),'--')
% 
% figure, plot(1:30,ch_st121_nor)
% hold on, plot(1:30,ch_st121_les,'--')
% 
% mean(ch_nAs_nor(:,:,1,2),2)-mean(ch_nAs_les(:,:,1,2),2)
% mean(ch_A_nor(:,:,1,2),2)-mean(ch_A_les(:,:,1,2),2)
% mean(ch_As_nor(:,:,1,2),2)-mean(ch_As_les(:,:,1,2),2)
% mean(ch_st121_nor,2)-mean(ch_st121_les,2)
% 
% covnoles=cov(ch_nAs_nor(:,:,1,2)');covles=cov(ch_nAs_les(:,:,1,2)');
% 0.5*(covnoles+covles)
% covnoles=cov(ch_A_nor(:,:,1,2)');covles=cov(ch_A_les(:,:,1,2)');
% 0.5*(covnoles+covles)
% covnoles=cov(ch_As_nor(:,:,1,2)');covles=cov(ch_As_les(:,:,1,2)');
% 0.5*(covnoles+covles)
% covnoles=cov(ch_st121_nor');covles=cov(ch_st121_les');
% 0.5*(covnoles+covles)