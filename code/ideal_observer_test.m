%ideal observer (pre-whitening template mathching)
% load ncatlestew_Im_idealUm
% Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% lesion_g1=Im_ideal(23:52,16:43,29:48,1);
% load ncat16gtew_Im_idealUm
% Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
% normal_g1=Im_ideal(23:52,16:43,29:48,1);

load ncat_g1myo_index ind_myo
for g=6
    filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
end
temp=reshape(temp,64,64,64);
normal_g1=temp(23:52,16:43,29:48);
immax=max(normal_g1(ind_myo));
normal_g1=normal_g1/immax;
for g=6
    filename=['C:\simind\ncat16gles_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
end
temp=reshape(temp,64,64,64);
lesion_g1=temp(23:52,16:43,29:48);
immax=max(lesion_g1(ind_myo));
lesion_g1=lesion_g1/immax;

frame=8;
sig_temp=normal_g1(:,:,frame)-lesion_g1(:,:,frame);
% sig_temp=sig_temp(23:25,11:13);
%ST-121
N=30;
c_gate=1;l_gate=16;r_gate=2;
nor_dir_name='D:\imagereconstruction\4D\data_16g_fbp';
les_dir_name='D:\imagereconstruction\4D\data_les_fbp';
ch_st121_les=zeros(N,1);ch_st121_nor=zeros(N,1);
temp_les=zeros(30,28,30);temp_nor=zeros(30,28,30);
for n=1:N
    filename=[les_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,:,c_gate)+.25*Im_fbp(:,:,:,l_gate)+.25*Im_fbp(:,:,:,r_gate);
    immax=max(temp(ind_myo));
    temp=temp/immax;
    temp_les(:,:,n)=temp(:,:,frame);
    
    filename=[nor_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,:,c_gate)+.25*Im_fbp(:,:,:,l_gate)+.25*Im_fbp(:,:,:,r_gate);
    immax=max(temp(ind_myo));
    temp=temp/immax;
    temp_nor(:,:,n)=temp(:,:,frame);
    
end

temp_les=reshape(temp_les,[30*28 30]);
temp_nor=reshape(temp_nor,[30*28 30]);
% % temp_nor=temp_nor(23:25,11:13,:);
% % temp_les=temp_les(23:25,11:13,:);
% % temp_les=reshape(temp_les,[9 30]);
% % temp_nor=reshape(temp_nor,[9 30]);
% % c_temp_les=cov(temp_les');
% % inv_c_les=pinv(c_temp_les);
% % c_temp_nor=cov(temp_nor');
% % inv_c_nor=pinv(c_temp_nor);
% temp=cat(2,temp_les,temp_nor);
% c_temp=cov(temp');
% inv_c=pinv(c_temp);

ch_st121_nor=sig_temp(:)'*temp_nor;%*inv_c
ch_st121_les=sig_temp(:)'*temp_les;%*inv_c
snr_st121=sqrt((mean(ch_st121_les)-mean(ch_st121_nor))^2/(.5*var(ch_st121_les)+.5*var(ch_st121_nor)));

% snr_st121=sig_temp(:)'*inv_c*sig_temp(:);
% Az_st121=0.5*(1+erf(snr_id/2));

%map
N=30;
%c_gate=1;l_gate=16;r_gate=2;
nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
les_dir_name='D:\imagereconstruction\4D\data_les_5st';
t_index=[0 2 3];

for t=1:3
    for s=1:5
        temp_les=zeros(30,28,30);temp_nor=zeros(30,28,30);
        for n=1:N
            filename=[les_dir_name num2str(t_index(t)) '\Im_ncatles_A_s' num2str(s) 't' num2str(t_index(t)) '_n'...
                num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,:,1);
            immax=max(temp(ind_myo));
            temp=temp/immax;
            temp_les(:,:,n)=temp(:,:,frame);

            filename=[nor_dir_name num2str(t_index(t)) '\Im_ncat16g_A_s' num2str(s) 't' num2str(t_index(t)) '_n'...
                num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,:,1);
            immax=max(temp(ind_myo));            
            temp=temp/immax;
            temp_nor(:,:,n)=temp(:,:,frame);

        end
        
        temp_les=reshape(temp_les,[30*28 30]);
        temp_nor=reshape(temp_nor,[30*28 30]);
% %         temp_nor=temp_nor(23:25,11:13,:);
% %         temp_les=temp_les(23:25,11:13,:);
% %         temp_les=reshape(temp_les,[9 30]);
% %         c_temp_les=cov(temp_les');
% %         inv_c_les=pinv(c_temp_les);        
% %         temp_nor=reshape(temp_nor,[9 30]);
% %         c_temp_nor=cov(temp_nor');
% %         inv_c_nor=pinv(c_temp_nor);
%         temp=cat(2,temp_les,temp_nor);
%         c_mapt=cov(temp');
%         inv_mapt=pinv(c_mapt);

        ch_mapt_nor=sig_temp(:)'*temp_nor;%*inv_mapt
        ch_mapt_les=sig_temp(:)'*temp_les;%*inv_mapt
        snr_mapA(t,s)=sqrt((mean(ch_mapt_les)-mean(ch_mapt_nor))^2/(.5*var(ch_mapt_les)+.5*var(ch_mapt_nor)));
%         snr_mapA(t,s)=sig_temp(:)'*inv_mapt*sig_temp(:);
    end
end
%Az_mapt=0.5*(1+erf(snr_id/2));

for t=1:3
    for s=1:5
        temp_les=zeros(30,28,30);temp_nor=zeros(30,28,30);
        for n=1:N
            filename=[les_dir_name num2str(t_index(t)) '\Im_ncatles_AS_s' num2str(s) 't' num2str(t_index(t)) '_n'...
                num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,:,1);
            immax=max(temp(ind_myo));
            temp=temp/immax;
            temp_les(:,:,n)=temp(:,:,frame);

            filename=[nor_dir_name num2str(t_index(t)) '\Im_ncat16g_AS_s' num2str(s) 't' num2str(t_index(t)) '_n'...
                num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,:,1);
            immax=max(temp(ind_myo));            
            temp=temp/immax;
            temp_nor(:,:,n)=temp(:,:,frame);

        end
        
        temp_les=reshape(temp_les,[30*28 30]);
        temp_nor=reshape(temp_nor,[30*28 30]);
% %         temp_nor=temp_nor(23:25,11:13,:);
% %         temp_les=temp_les(23:25,11:13,:);
% %         temp_les=reshape(temp_les,[9 30]);
% %         c_temp_les=cov(temp_les');
% %         inv_c_les=pinv(c_temp_les);        
% %         temp_nor=reshape(temp_nor,[9 30]);
% %         c_temp_nor=cov(temp_nor');
% %         inv_c_nor=pinv(c_temp_nor);
%         temp=cat(2,temp_les,temp_nor);
%         c_mapt=cov(temp');
%         inv_mapt=pinv(c_mapt);
        
        ch_mapt_nor=sig_temp(:)'*temp_nor;%*inv_mapt
        ch_mapt_les=sig_temp(:)'*temp_les;%*inv_mapt
        
        snr_mapAS(t,s)=sqrt((mean(ch_mapt_les)-mean(ch_mapt_nor))^2/(.5*var(ch_mapt_les)+.5*var(ch_mapt_nor)));
        %sig_temp(:)'*inv_mapt*sig_temp(:);%
    end
end
% save IO_snr snr_st121 snr_mapA snr_mapAS %good (pre-white)
save IO_snr_n snr_st121 snr_mapA snr_mapAS %no-whitening, not good for MAP-S

load IO_snr_n snr_st121 snr_mapA snr_mapAS
az_st121=snr_st121;%0.5*(1+erf(snr_st121/2));
az_mapA=snr_mapA;%0.5*(1+erf(snr_mapA/2));
az_mapAS=snr_mapAS;%0.5*(1+erf(snr_mapAS/2));

sbeta=[0 5e-5 1e-4 1e-3 3e-3];
figure,plot(sbeta,az_mapA(1,:),'*-',sbeta,az_mapA(2,:),'+-',sbeta,az_mapA(3,:),'o-')
hold on,plot(sbeta,repmat(az_st121,1,5),'k-')
legend('\beta_t=0','\beta_t=1\times 10^{-3}','\beta_t=2\times 10^{-3}','ST 121')
title('AC'),xlabel('\beta_s'),ylabel('Az')

sbeta=[0 5e-4 1e-3 2e-3 4e-3];
figure,plot(sbeta,az_mapAS(1,:),'*-',sbeta,az_mapAS(2,:),'+-',sbeta,az_mapAS(3,:),'o-')
hold on,plot(sbeta,repmat(az_st121,1,5),'k-')
legend('\beta_t=0','\beta_t=2\times 10^{-3}','\beta_t=3\times 10^{-3}','ST 121')
title('ASC'),xlabel('\beta_s'),ylabel('Az')