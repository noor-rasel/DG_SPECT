%spatial and frequency filtering

for g=1
    filename=['C:\simind\ncat16g_act_' num2str(g) '.bin'];
    fid=fopen(filename,'rb');
    temp=fread(fid,64^3,'single');
    fclose(fid);
end
temp=reshape(temp,64,64,64);
temp=temp(:,:,37)';
dsp(temp)

Fi=Qfreq(64);Fi=Fi(:,:,1);

%frequency filtering
tempF=fft2(temp);
tempFF=tempF.*fftshift(Fi(:,:,1));
tempIF=ifft2(tempFF);
dsp(real(tempIF))

%spatial filtering
iFi=ifft2(fftshift(Fi));
tempC=conv2(temp,fftshift(iFi),'same');
dsp(real(tempC))