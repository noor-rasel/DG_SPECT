loc=[-71 25];frame=8;%<=>[198 104](center [129 129])
loc=[-72 33];%<=>[197 96] _loc2
gate=1;
iFi=Qfreq_s(256,loc);iFi=reshape(iFi,[256*256 4])';

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
            temp=postproc_cho(temp,256);
            ch_nAs_les(:,n,s,t)=iFi*temp(:);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_nAS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
            temp=postproc_cho(temp,256);
            ch_nAs_nor(:,n,s,t)=iFi*temp(:);
        end
        Az_nAs(s,t)=Azfinder(squeeze(ch_nAs_nor(:,:,s,t))',squeeze(ch_nAs_les(:,:,s,t))',int_noise);
    end
end
save cho_result_nAS256_g1f8_2 ch_nAs_les ch_nAs_nor Az_nAs
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
            temp=postproc_cho(temp,256);
            ch_As_les(:,n,s,t)=iFi*temp(:);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_AS_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
            temp=postproc_cho(temp,256);
            ch_As_nor(:,n,s,t)=iFi*temp(:);
        end
        Az_As(s,t)=Azfinder(squeeze(ch_As_nor(:,:,s,t))',squeeze(ch_As_les(:,:,s,t))',int_noise);
    end
end
save cho_result_AS256_g1f8_2 ch_As_les ch_As_nor Az_As
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
            temp=postproc_cho(temp,256);
            ch_A_les(:,n,s,t)=iFi*temp(:);
            filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_A_s'...
                num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
            load(filename,'Im_maps');
            temp=Im_maps(:,:,frame,gate);
            temp=postproc_cho(temp,256);
            ch_A_nor(:,n,s,t)=iFi*temp(:);
        end
        Az_A(s,t)=Azfinder(squeeze(ch_A_nor(:,:,s,t))',squeeze(ch_A_les(:,:,s,t))',int_noise);
    end
end
save cho_result_A256_g1f8_2 ch_A_les ch_A_nor Az_A
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
    temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);
    temp=postproc_cho(temp,256);
    ch_st121_les(:,n)=iFi*temp(:);
    filename=[nor_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);
    temp=temp/1.1189;
    temp=postproc_cho(temp,256);
    ch_st121_nor(:,n)=iFi*temp(:);
end
Az_st121=Azfinder(ch_st121_nor(:,:)',ch_st121_les(:,:)',int_noise);
save cho_result_st121_256_g1f8 ch_st121_les ch_st121_nor Az_st121
%save cho_result_st121_256_g1f8_2 ch_st121_les ch_st121_nor Az_st121