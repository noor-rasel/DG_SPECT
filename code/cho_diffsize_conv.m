temp=zeros(64);
temp(23:52,16:43)=ideal_f8;
iFi=Qfreq_s(64,0);
for n=1:4
    temp_f(:,:,n)=conv2(iFi(:,:,n),temp,'same');
end
for n=1:4
    dsp(temp_f(25:52,16:43,n))
end


temp=postproc_cho(ideal_f8,140);
iFi=Qfreq_s(140,0);
for n=1:4
    temp_f2(:,:,n)=conv2(iFi(:,:,n),temp,'same');
end
for n=1:4
    dsp(temp_f2(:,:,n))
end