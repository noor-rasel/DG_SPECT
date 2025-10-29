%load
fid=fopen('D:\imagereconstruction\4D\ncat\wholencat16g_act_1.bin','rb');
temp=fread(fid,'single');
fclose(fid);
tempa=reshape(temp,[64 64 64]);
% dsp(temp)
fid=fopen('D:\imagereconstruction\4D\ncat\wholencat16g_atn_1.bin','rb');
temp=fread(fid,'single');
fclose(fid);
temp=reshape(temp,[64 64 64]);
% dsp(temp)

dsp(temp(:,:,34)')
dsp(tempa(:,:,34)')
%transverse atenuation map
[xi yi]=meshgrid(1:63/255:64,1:63/255:64);
temp=temp(:,:,34)';
temp=interp2(temp,xi,yi);
temp=temp/max(temp(:));
% temp=temp.^0.6;
dsp(temp)
tran=temp;tran=uint8(tran*255);
imwrite(tran,'ncat_tran34atn.tif','tiff');
%%activity
temp=tempa(:,:,34)';
temp=interp2(temp,xi,yi);
temp=temp/max(temp(:));
% temp=temp.^0.6;
dsp(temp)
tran=temp;tran=uint8(tran*255);
imwrite(tran,'ncat_tran34act.tif','tiff');


%coronal
temp=squeeze(x_phantom(:,36,:,1))';
temp=interp2(temp,2);
temp=temp/max(temp(:));
temp=temp.^0.6;
dsp(temp)
coro=temp;coro=uint8(coro*255);
imwrite(coro,'ncat_coro.tif','tiff');
%sagittal
temp=squeeze(x_phantom(36,:,:,1))';
temp=interp2(temp,2);
temp=temp/max(temp(:));
temp=temp.^0.6;
dsp(temp)
sagi=temp;sagi=uint8(sagi*255);
imwrite(sagi,'ncat_sagi.tif','tiff');

%%%%%%%%%%%%%%%%%%%%%
%raw data
ncat_act=zeros(64,64,64,16);
for g=1:16    
    filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    ncat_act(:,:,:,g)=temp;
end
ncat_atn=zeros(64,64,64,16);
for g=1:16
    filename=['C:\simind\ncat16g_atn_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
    temp=reshape(temp,64,64,64);
    ncat_atn(:,:,:,g)=temp;
end