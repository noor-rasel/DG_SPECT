%phantom
load ncatlestew_NL sino_p
dsp(flipud(sino_p(:,:,63,1)))

fid=fopen('C:\simind\ncatlestew1a.a02','rb');
temp=fread(fid,'single');
fclose(fid);
temp=reshape(temp,64,64,64);
dsp(flipud(rot90(temp(:,:,1))))
fid=fopen('C:\simind\ncatlestew1s.a02','rb');
temps=fread(fid,'single');
fclose(fid);
temps=reshape(temps,64,64,64);
dsp(flipud(rot90(temps(:,:,1))))
%Anterior projection of NCAT
pri=flipud(rot90(temp(:,:,1)-temps(:,:,1)));
sca=flipud(rot90(temps(:,:,1)));
antp=cat(2,pri,sca*3);
tot=pri+sca;
tot=tot/max(tot(:));
tot=uint8(tot*255);
imwrite(tot,'ncat_simind_ap.tif','tiff');
%roi
n=1;
fn=['C:\simind\ncat16gles_act_' num2str(n) '.bin'];
fid=fopen(fn,'rb');
temp=fread(fid,'single');
fclose(fid);
ncat_phantom=reshape(temp,64,64,64);
ncat_phantom=ncat_phantom(25:49,16:40,34:38);
ncat_phantom=ncat_phantom(:,:,3);
dsp(ncat_phantom')
ncat_phantom=ncat_phantom';ncat_phantom(11,23)=ncat_phantom(11,23)*.5;
ncat_phantom=ncat_phantom/max(ncat_phantom(:));
ncat_phantom=uint8(ncat_phantom*255);
imwrite(ncat_phantom,'ncat_isbi07_4droi.tif','tiff');
%SNR
m_name={'ncatlestew_iMCnASC_Um';'ncatlestew_iMC_AC_Um';'ncatlestew_iMC_ASC_Um';...
    'ncatlestew_mapsnASC_Um';'ncatlestew_maps_AC_Um';'ncatlestew_maps_ASC_Um'};
snr_name={'snr_t';'snr_tac';'snr_tasc';'snr_s';'snr_sac';'snr_sasc'};
%first part
snr_comp=zeros(6,10);
n=10;
for s=1:6
        filename=['.\isbi07\umass\' m_name{s} num2str(n) '.mat'];        
        load(filename,snr_name{s});
        snr_comp(s,:)=eval(snr_name{s});
end

[max_snr,max_snr_i]=max(snr_comp,[],2)
%     9.2229 6     
%    10.1118 6     
%    10.6176 6     
%     8.3807 3     
%     9.0206 4     
%     9.5776 6

%TEW
m_name={'ncatlestewE_iMC_ASC_Um';...
    'ncatlestewE_maps_ASC_Um'};
snr_name={'snr_tascE';'snr_sascE'};
%first part
snr_compE=zeros(2,10);
n=10;
for s=1:2
        filename=['.\isbi07\umass\' m_name{s} num2str(n) '.mat'];        
        load(filename,snr_name{s});
        snr_compE(s,:)=eval(snr_name{s});
end

[max_snr,max_snr_i]=max(snr_compE,[],2)
% 9.9970  8
% 7.9012  10

%images & SNR
load ncatlestew_Im_idealUm
Im_ideal=Im_ideal*8e6/sum(Im_ideal(:));
Im_ideal=Im_ideal(25:49,16:40,34:38,:);%ideal_c=sum(Im_ideal(:));
Im_ideal=squeeze(Im_ideal(:,:,3,:));
Im_ideal=permute(Im_ideal,[2 1 3]);
dsp(Im_ideal)
m_name={'ncatlestew_mapsnASC_Um';'ncatlestew_maps_ASC_Um';'ncatlestew_maps_ASC_Um';...    
    'ncatlestew_iMCnASC_Um';'ncatlestewE_iMC_ASC_Um';'ncatlestew_iMC_ASC_Um';...
    };%'ncatlestew_maps_AC_Um';'ncatlestew_iMC_AC_Um';
max_snr_i=[3 10 6 6 8 6];
Im_comb=zeros(25,25,16,6);
for s=1:6
    filename=['.\isbi07\umass\' m_name{s} num2str(max_snr_i(s)) '.mat'];
    load(filename,'Im_maps')
%     Im_maps=Im_maps*ideal_c/sum(Im_maps(:));
% if s==2 | s==3 | s==5 | s==6
%     Im_maps=Im_maps*.8;
% end
    temp=squeeze(Im_maps(:,:,3,:));temp=permute(temp,[2 1 3]);
    snr_slice(s)=10*log10(Im_ideal(:)'*Im_ideal(:)/sum((Im_ideal(:)-temp(:)).^2));
    Im_comb(:,:,:,s)=temp;
end
%SNR
% 9.6910   10.9943   12.0794   10.8557    8.4121    9.2308   10.2149   10.1173
Im_comb=cat(4,Im_ideal,Im_comb);
Im_show=zeros(25*6,25*7);
for n=1:7
    temp=Im_comb(:,:,1:3:16,n);
    temp=permute(temp,[1 3 2]);
    temp=reshape(temp,[25*6 25]);
    Im_show(:,(n-1)*25+1:n*25)=temp;
end
figure('position',[100 100 700 600])
imagesc(interp2(Im_show,2)),clinicalcolor,axis equal,axis tight
hold on, plot(102,1:597,'w.',402,1:597,'w.')
set(gca,'Xtick',[150:100:697],'Ytick',[50:100:597],'XTickLabel',...
    {'NC';'A+S';'A+PS';'NC';'A+S';'A+PS'},...
    'YTickLabel',{'1';'4';'7';'10';'13';'16'});
title('Ideal                                  MAP-S                                                         MAP-T                 ')
% title('Ideal                                  MAP-T                                                         MAP-S                 ')

%TACs
%dsp(Im_ideal(:,:,1),1)%[11:12 22:23] fixed ROI
%tip%[17:18,19:20]
%defect ROI: good for paper
for n=1:16
    fn=['D:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesion16g_act_' num2str(n) '.bin'];
    fid=fopen(fn,'rb');
    temp=fread(fid,'single');
    fclose(fid);
    temp=temp>0;    
    temp=reshape(temp,64,64,64);
    lesion(:,:,:,n)=temp(25:49,16:40,34:38);%ideal_c=sum(Im_ideal(:));
end
lesion=squeeze(lesion(:,:,3,:));
lesion=permute(lesion,[2 1 3]);
for n=1:16
    nz_def(n)=nnz(lesion(:,:,n));%max(temp)*.7
end
for n=1:7
%     tac_rec(:,n)=squeeze(mean(mean(Im_comb(17:18,19:20,:,n))));
    tac_rec(:,n)=squeeze(sum(sum(Im_comb(:,:,:,n).*lesion))/nz_def(n));
end
figure,plot(1:16,tac_rec(:,1)/max(tac_rec(:,1)),'s-',1:16,tac_rec(:,5)/max(tac_rec(:,5)),'^-',...
    1:16,tac_rec(:,9)/max(tac_rec(:,9)),'+-')
%1:16,tac_rec(:,2)/max(tac_rec(:,2)),'o-',...
%1:16,tac_rec(:,6)/max(tac_rec(:,6)),'*-',
axis([1 16 0.2 1.1]),xlabel('Gate frame number'),ylabel('Normalized counts')
legend('Ideal','MAP-T ASTC','MAP-S ASTC')
mean((tac_rec(:,2:end)-repmat(tac_rec(:,1),1,8)).^2)
figure,plot(1:16,tac_rec(:,1)/max(tac_rec(:,1)),'s-',1:16,tac_rec(:,2)/max(tac_rec(:,1)),'o-',...
    1:16,tac_rec(:,5)/max(tac_rec(:,1)),'+-')
    %     1:16,tac_rec(:,3)/max(tac_rec(:,3)),'*-',1:16,tac_rec(:,4)/max(ta
    %     c_rec(:,4)),'^-',...
figure,plot(1:16,tac_rec(:,1)/max(tac_rec(:,1)),'s-',1:16,tac_rec(:,6)/max(tac_rec(:,6)),'o-',...
%     1:16,tac_rec(:,7)/max(tac_rec(:,7)),'*-',1:16,tac_rec(:,8)/max(tac_rec(:,8)),'^-',...
    1:16,tac_rec(:,9)/max(tac_rec(:,9)),'+-')
corrcoef(tac_rec)
%mse
%12.9013   10.9543    3.8581    4.2482   19.9498   16.8162    7.6161    7.8561
%on defect
%0.4761   15.8737    2.4111    2.7614   30.3187   15.7751    4.9284   18.6393