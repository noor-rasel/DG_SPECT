
% loc=[-14 5];%[25 13]%real [i,j]=[11-loc(1) 18-loc(2)]
frame=8;
% loc=[-12 7];

nor_dir_name='D:\imagereconstruction\4D\data_16g_5st';
les_dir_name='D:\imagereconstruction\4D\data_les_5st';
t_ind=[0 2 3];
S=5;T=3;N=30;gate=1;
int_noise=0.0011;
s=1;t=2;
Az_A=zeros(S,T);
ch_A_les=zeros(4,N);ch_A_nor=zeros(4,N);
for n=1:N
    filename=[les_dir_name num2str(t_ind(t)) '\Im_ncatles_A_s'...
        num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
    load(filename,'Im_maps');
    temp=Im_maps(:,:,frame,gate);
%     immax=max(temp(ind_myo));
    im_mapnAl(:,:,n)=temp;%/immax;
%     ch_A_les(:,n)=cho_feature(Fi,temp,loc);
    filename=[nor_dir_name num2str(t_ind(t)) '\Im_ncat16g_A_s'...
        num2str(s) 't' num2str(t_ind(t)) '_n' num2str(n) '.mat'];
    load(filename,'Im_maps');
    temp=Im_maps(:,:,frame,gate);
%     immax=max(temp(ind_myo));
    im_mapnAn(:,:,n)=temp;%/immax;
%     ch_A_nor(:,n)=cho_feature(Fi,temp,loc);
end
Az_A=Azfinder(ch_A_nor(2:4,:)',ch_A_les(2:4,:)',int_noise);

% loc=[-8 10];frame=8; % [19 8]
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
    temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);
    im_st121_l(:,:,n)=temp;
%     ch_st121_les(:,n)=cho_feature(Fi,temp,loc);
    filename=[nor_dir_name '\Im_fbp_n'...
        num2str(n) '.mat'];
    load(filename,'Im_fbp');
    temp=.5*Im_fbp(:,:,frame,cen_gate)+.25*Im_fbp(:,:,frame,l_gate)+.25*Im_fbp(:,:,frame,r_gate);    
    im_st121_n(:,:,n)=temp/1.1189;%multiply is better than division.; no!!
%     ch_st121_nor(:,n)=cho_feature(Fi,temp,loc);
end
Az_st121=Azfinder(ch_st121_nor(2:4,:)',ch_st121_les(2:4,:)',int_noise);

%256
loc=[-70 25];%<=>[199 104]
loc=[-71 25];%<=>[198 104]
iFi=Qfreq_s(256,loc);iFi=reshape(iFi,[256*256 4])';
tic
for n=1:30
%     temp=postproc_cho(im_mapnAl(:,:,n),256);
%     ch_nA_les(:,n)=iFi*temp(:);
%     temp=postproc_cho(im_mapnAn(:,:,n),256);
%     ch_nA_nor(:,n)=iFi*temp(:);
%     temp=postproc_cho(im_mapAl(:,:,n),256);
%     ch_A_les(:,n)=iFi*temp(:);
%     temp=postproc_cho(im_mapAn(:,:,n),256);
%     ch_A_nor(:,n)=iFi*temp(:);
    temp=postproc_cho(im_st121_l(:,:,n),256);
    ch_st121_les(:,n)=iFi*temp(:);
    temp=postproc_cho(im_st121_n(:,:,n),256);
    ch_st121_nor(:,n)=iFi*temp(:);
end
toc%13.5sec
figure,plot(1:30,ch_A_nor(1,:),1:30,ch_A_les(1,:))

temp_nl=postproc_cho(mean(im_mapnAl,3),256);
temp_nn=postproc_cho(mean(im_mapnAn,3),256);
temp_l=postproc_cho(mean(im_mapAl,3),256);
temp_n=postproc_cho(mean(im_mapAn,3),256);
temp_sl=postproc_cho(mean(im_st121_l,3),256);
temp_sn=postproc_cho(mean(im_st121_n,3),256);

dsp(temp_n(198-10:198+10,104-10:104+10)-temp_l(198-10:198+10,104-10:104+10),1,0)

dsp((temp_n(198-10:198+10,104-10:104+10)-temp_l(198-10:198+10,104-10:104+10))-...
	(temp_nn(198-10:198+10,104-10:104+10)-temp_nl(198-10:198+10,104-10:104+10)),1,0)
dsp(temp_n(198-10:198+10,104-10:104+10),1,0)

dsp((temp_n(198-10:198+10,104-10:104+10)-temp_l(198-10:198+10,104-10:104+10))-...
	(temp_sn(198-10:198+10,104-10:104+10)-temp_sl(198-10:198+10,104-10:104+10)),1,0)

dsp(temp_sn(198-10:198+10,104-10:104+10)-temp_sl(198-10:198+10,104-10:104+10),1,0)