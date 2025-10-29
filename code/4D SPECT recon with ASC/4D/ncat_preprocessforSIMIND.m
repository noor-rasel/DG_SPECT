%get rid of the activity and intensity outside of the ROI
load roi
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
%     filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
%     fid=fopen(filename,'rb');
%     temp1=fread(fid,64^3,'single');
%     fclose(fid);
%     temp1=reshape(temp1,64,64,64);
%     dsp(temp-temp1);pause;close
end
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
%     filename=['C:\simind\ncat16g_atn_' num2str(g) '.bin'];
%     fid=fopen(filename,'rb');
%     temp1=fread(fid,64^3,'single');
%     fclose(fid);
%     temp1=reshape(temp1,64,64,64);
%     dsp(temp-temp1);pause;close
end
%CT motion
for g=1:16    
    filename=['C:\simind\ncat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
	slice35(:,:,g)=temp(:,:,35);
    %slice35(:,:,g)=fliplr(rot90(temp(:,:,35),-1));
    
%     filename=['C:\simind\ncat16g_atn_' num2str(g) '.bin'];
%     fid=fopen(filename,'rb');
%     temp1=fread(fid,64^3,'single');
%     fclose(fid);
%     temp1=reshape(temp1,64,64,64);
%     dsp(temp-temp1);pause;close
end
%phantom with lesion
intensity_reduce_ratio=0.7;
load roi
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\lesion16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp1=fread(fid,64^3,'single');
    fclose(fid);
    temp1=reshape(temp1,64,64,64);
    temp=temp1-intensity_reduce_ratio*temp;
    ni(g)=any(find(temp<0));%negative value monitor
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncat16gles_act_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncat16gles_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end

%phantom with subtle lesion (anterior-lateral wall 3 pixels, 30 degree)
intensity_reduce_ratio=0.14;
load roi
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\lesion2_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);%dsp(temp(23:52,16:43,29:48),1);pause;close
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp1=fread(fid,64^3,'single');
    fclose(fid);
    temp1=reshape(temp1,64,64,64);
    temp=temp1-intensity_reduce_ratio*temp;
    ni(g)=any(find(temp<0));%negative value monitor
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncatles2_act_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncatles2_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
%show lesion
for g=1:16
    filename=['C:\simind\ncatles2_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
     temp=reshape(temp,64,64,64);
    dsp(temp(23:52,16:43,29:48),1);pause;close
end


%phantom with King's LAD lesion (anterior wall 4 pixels, 45 degree) 2002 TNS
intensity_reduce_ratio=0.35;
load roi
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\lesion_king_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);%dsp(temp(23:52,16:43,29:48),1);pause;close
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp1=fread(fid,64^3,'single');
    fclose(fid);
    temp1=reshape(temp1,64,64,64);
    temp=temp1-intensity_reduce_ratio*temp;
    ni(g)=any(find(temp<0));%negative value monitor
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncatlesK_act_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncatlesK_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
%show lesion
for g=1:16
    filename=['C:\simind\ncatlesK_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
     temp=reshape(temp,64,64,64);
%     dsp(temp(23:52,16:43,29:48),1);pause;close
dsp(temp(23:52,16:43,39),1);im(g)=getframe;pause(.1);close
end
movie(im,10)%transverse frame 39; short axis frame 40!!!(all data)

%phantom with Inferior Base lesion (inferior wall 4 pixels, 45 degree)
intensity_reduce_ratio=0.35;
load roi
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\lesionIB_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);%dsp(temp(23:52,16:43,29:48),1);pause;close
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp1=fread(fid,64^3,'single');
    fclose(fid);
    temp1=reshape(temp1,64,64,64);
    temp=temp1-intensity_reduce_ratio*temp;
    ni(g)=any(find(temp<0));%negative value monitor
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncatlesIB_act_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
for g=1:16
    filename=['D:\imagereconstruction\4D\ncat\wholencat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    temp=temp.*repmat(roi,[1 1 64]);
    filename=['C:\simind\ncatlesIB_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
end
%show lesion
for g=1:16
    filename=['C:\simind\ncatlesIB_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
     temp=reshape(temp,64,64,64);
    dsp(temp(23:52,16:43,29:48),1);pause;close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%128
filename='D:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat128_act_1.bin';
fid=fopen(filename,'rb');
temp=fread(fid,128^2*64,'single');
fclose(fid);
temp=reshape(temp,[128 128 64]);

% filename='D:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\wholencat16g_act_1.bin';
filename='C:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\wholencat16g_act_1.bin';
fid=fopen(filename,'rb');
temp64=fread(fid,64^3,'single');
fclose(fid);
temp64=reshape(temp64,[64 64 64]);
% temp64_s=cardiac3drot(temp64);
load ncat256_down64 ncat256_down64

filename='c:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat256_act_1.bin';
filename='C:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesion_S256_act_1.bin';
fid=fopen(filename,'rb');
temp=fread(fid,256^3,'single');
fclose(fid);
temp=reshape(temp,[256 256 256]);
temp=temp(1:4:256,1:4:256,1:4:256);dsp(temp(:,:,33:41+7))
lesion256=temp;
dsp(temp(:,:,33:41+7)-lesion256(:,:,33:41+7))
%tic,temp_s=cardiac3drot(temp);toc

filename='C:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesion_king_act_1.bin';
filename='C:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesionIB_act_1.bin';
fid=fopen(filename,'rb');
lesion64=fread(fid,64^3,'single');
fclose(fid);
lesion64=reshape(lesion64,[64 64 64]);

dsp(ncat256_down64(:,:,33:41+7))
dsp(ncat256_down64(:,:,33:41+7)-temp(:,:,33:41+7))
dsp(temp64(:,:,33-4:41+7-4))
dsp(temp64(:,:,33-4:41+7-4)-lesion64(:,:,33-4:41+7-4))

% %downsample
order=21;
[f1,f2] = freqspace(order,'meshgrid');
Hd = ones(order);
r = sqrt(f1.^2 + f2.^2);
Hd((r>0.25)) = 0;
h = fwind1(Hd,hamming(order));
ham21_cf25=h(11,:)/sum(h(11,:));

A=zeros(21,1,1);
A(:,1,1)=ham21_cf25;
B=convn(temp,A,'same');
A=zeros(1,21,1);
A(1,:,1)=ham21_cf25;
B=convn(B,A,'same');
A=zeros(1,1,21);
A(1,1,:)=ham21_cf25;
B=convn(B,A,'same');
lesion256d64=B(1:4:256,1:4:256,1:4:256);
ncat256d64=B(1:4:256,1:4:256,1:4:256);
dsp(ncat256d64(:,:,28:28+15)-lesion256d64(:,:,28:28+15))
% 
% save ncat256_lpf B
% ncat256_down64=B(1:4:256,1:4:256,1:4:256);
% save ncat256_down64 ncat256_down64

prefix={'c:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat256_act_';
    'C:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesion_S256_act_';
    'c:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat256_atn_'}
postfix='.bin';

intensity_reduce_ratio=0.20;
load roi
scal=10;
% for m=1:3
tic
    for g=1:16
        filename=['c:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat256_act_' num2str(g) '.bin'];
        fid=fopen(filename,'rb');
        temp=fread(fid,256^3,'single');
        fclose(fid);
        temp=reshape(temp,[256 256 256]);
        temp=ideal_downsample(temp);
        temp=temp.*repmat(roi,[1 1 64])*scal;
        filename=['C:\simind\ncatd64_act_' num2str(g) '.bin'];
        fid=fopen(filename,'wb');
        fwrite(fid,temp,'single');
        fclose(fid);
        txt=['Normal activity gate # ' num2str(g) ' finished...'];
        disp(txt);
            
        filename=['c:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\lesion_S256_act_' num2str(g) '.bin'];
        fid=fopen(filename,'rb');
        temp1=fread(fid,256^3,'single');
        fclose(fid);
        temp1=reshape(temp1,[256 256 256]);
        temp1=ideal_downsample(temp1);
        temp1=temp1.*repmat(roi,[1 1 64])*scal;
        temp=temp-intensity_reduce_ratio*temp1;
        ni(g)=any(find(temp<0));%negative value monitor
        filename=['C:\simind\ncatd64l_act_' num2str(g) '.bin'];
        fid=fopen(filename,'wb');
        fwrite(fid,temp,'single');
        fclose(fid);
        txt=['Defect activity gate # ' num2str(g) ' finished...'];
        disp(txt);
        
        filename=['c:\imagereconstruction\ET_PHANTOM\NCAT\ncat_routines\ncat256_atn_' num2str(g) '.bin'];
        fid=fopen(filename,'rb');
        temp=fread(fid,256^3,'single');
        fclose(fid);
        temp=reshape(temp,[256 256 256]);
        temp=ideal_downsample(temp);
        temp=temp.*repmat(roi,[1 1 64]);
        filename=['C:\simind\ncatd64l_atn_' num2str(g) '.bin'];
        fid=fopen(filename,'wb');
        fwrite(fid,temp,'single');
        fclose(fid);
        txt=['Attenuation gate # ' num2str(g) ' finished...'];
        disp(txt);
    end
    toc
% end
% get rid of negative values
for g=1:16
    filename=['C:\simind\ncatd64_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=abs(temp);
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
    txt=['Normal activity gate # ' num2str(g) ' finished...'];
    disp(txt);

    filename=['C:\simind\ncatd64l_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=abs(temp);
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
    txt=['Defect activity gate # ' num2str(g) ' finished...'];
    disp(txt);

    filename=['C:\simind\ncatd64l_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=abs(temp);
    fid=fopen(filename,'wb');
    fwrite(fid,temp,'single');
    fclose(fid);
    txt=['Attenuation gate # ' num2str(g) ' finished...'];
    disp(txt);
end

for g=1:16
    filename=['C:\simind\ncatd64l_atn_' num2str(g) '.bin'];
    outname=['C:\simind\ncatd64_atn_' num2str(g) '.bin'];
    copyfile(filename,outname);
    filename=['C:\simind\ncatD64_act_' num2str(g) '.bin'];
    outname=['C:\simind\ncatd64_act_' num2str(g) '.bin'];
    copyfile(filename,outname);
end