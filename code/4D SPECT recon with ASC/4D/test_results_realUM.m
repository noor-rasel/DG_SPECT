load real_gate1_sino temp
sino=temp(23:86,:,:);
tc=sum(sino(:));

load fbp_realdata42 rim
% imir=rim(:,:,23:86);
imir=clinicfilt3d_yyjin(rim,2.4,.4);
imir=imir(:,:,23:86);
imir=imir*tc/sum(imir(:));
dsp(imir(:,:,20:28),1);

load osem5_attn
x_ac=x;
x_ac=x_ac*tc/sum(x_ac(:));
dsp(x_ac(:,:,20:28),1);
load osem5
x=x*tc/sum(x(:));
dsp(x(:,:,20:28),1);


load ./real_dataUM/realg1_map_s2 %1e-5
dsp(x(:,:,23:26),1);
load ./real_dataUM/realg1_map_s3 %5e-5
dsp(x(:,:,23:26),1);
load ./real_dataUM/realg1_map_s4 %1e-4
dsp(x(:,:,23:26),1);
load ./real_dataUM/realg1_map_s5 %5e-4
dsp(x(:,:,23:26),1);
%%%%%%%%%%%%%%%%%%%%
%FBP
y=zeros(128,128,64,8);
for g=1:8
    filename=['./real_dataUM/realg' num2str(g) '_fbp'];
    load(filename);
    for s=1:64
        y(:,:,s,g)=rot90(flipud(x(:,:,s)),1);
    end
end    
save real_data_fbp_oim y
%ST121
imir=reshape(y,[128^2*64 8]);
w_st121=[0.5 0.25 zeros(1,5) 0.25];G=8;
for n=1:G
    wm_st121(:,n)=circshift(w_st121',n-1);
end
imir121=imir*wm_st121;
imir121=reshape(imir121,[128 128 64 8]);
save real_data_st121 imir121
%MAP-S NC
y=zeros(128,128,64,8);
for g=1:8
    filename=['./real_dataUM/realg' num2str(g) '_mapNC_s5'];
    load(filename);
    for s=1:64
        y(:,:,s,g)=rot90(flipud(x(:,:,s)),1);
    end
end    
save real_data_mapsNC_oim y
%MAP-S AC
y=zeros(128,128,64,8);
for g=1:8
    filename=['./real_dataUM/realg' num2str(g) '_map_s4'];
    load(filename);
    for s=1:64
        y(:,:,s,g)=rot90(flipud(x(:,:,s)),1);
    end
end    
save real_data_mapsAC_oim y %1e-4

load real_data_mapsAC_oim y
yac=y;
load real_data_mapsNC_oim y
clear temp
slice=26;
for g=1:8
    temp(:,:,g)=y(52:81,59:88,slice,g);
    tempac(:,:,g)=yac(52:81,59:88,slice,g);
end
temp=temp/max(temp(:));
tempac=tempac/max(tempac(:));
temp=cat(2,temp,tempac);
for g=1:8
    figure('position',[100 100 600 350]), imagesc(interp2(temp(:,:,g),2),[0 1]), clinicalcolor, axis image
    mov(g)=getframe;pause(.1);close
end
figure('position',[100 100 600 350]),movie(mov,100)
%MAP-S ASC
ysas=zeros(128,128,64,8);
for g=1:8
    filename=['./real_dataUM/realg' num2str(g) '_mapASC_s4'];
    load(filename);
    for s=1:64
        ysas(:,:,s,g)=rot90(flipud(x(:,:,s)),1);
    end
end    
save real_data_mapsASC_oim ysas

%AC + MC 
load ./real_dataUM/real8gate_map_ST x
load ./real_dataUM/real8gate_map_s2t3 x
load ./real_dataUM/real8gate_map_s2t2 x
load ./real_dataUM/real8gate_map_s2t1 x

yst=zeros(128,128,64,8);
for g=1:8    
    for s=1:64
        yst(:,:,s,g)=rot90(flipud(x(:,:,s,g)),1);
    end
end    
% save real_data_maptAC_oim yt %from t3 (1e-4)
% save real_data_mapt2AC_oim yt2 %from t2 (5e-5)
save real_data_mapstAC_oim yst %from t2 (5e-5,5e-5)

%%%%%%%%%%%%%
%ASC+MC
load ./real_dataUM/real8gate_map_s2t2 x
yst=zeros(128,128,64,8);
for g=1:8    
    for s=1:64
        yst(:,:,s,g)=rot90(flipud(x(:,:,s,g)),1);
    end
end   
slice=26;
for g=1:8
    tempscst(:,:,g)=yst(52:81,59:88,slice,g);
end
% 
%composed st121 maps mapt2 mapt(3)
load real_data_st121 imir121
% load real_data_mapsAC_oim y
% load real_data_mapt2AC_oim yt2
% load real_data_maptAC_oim yt
load real_data_mapsASC_oim ysas
load real_data_mapsASC_s1t3 x
for g=1:8    
    for sli=1:64
        yt(:,:,sli,g)=rot90(flipud(x(:,:,sli,g)),1);
    end
end  
load real_data_mapsASC_s3t3 x
for g=1:8    
    for sli=1:64
        yst(:,:,sli,g)=rot90(flipud(x(:,:,sli,g)),1);
    end
end  

slice=26;
for g=1:8
    temp121(:,:,g)=imir121(52:81,59:88,slice,g);
    temps(:,:,g)=y(52:81,59:88,slice,g);
    tempt2(:,:,g)=yt2(52:81,59:88,slice,g);
    tempt(:,:,g)=yt(52:81,59:88,slice,g);
    tempsas(:,:,g)=ysas(52:81,59:88,slice,g);
end
tempsas=tempsas/max(tempsas(:));
temp121=temp121/max(temp121(:));
temps=temps/max(temps(:));
tempt2=tempt2/max(tempt2(:));
tempst=tempst/max(tempst(:));
tempt=tempt/max(tempt(:));
temp1=cat(2,temp121,temps);
temp2=cat(2,tempt2,tempst);
temp=cat(1,temp1,temp2);
for g=1:8
    figure('position',[100 100 600 650]), imagesc(interp2(temp(:,:,g),2),[0 1]), clinicalcolor, axis image
    axis off
    mov(g)=getframe;pause(.1);close
end
figure('position',[100 100 600 650]),movie(mov,100)

%images
for n=1:8
    im_com((n-1)*30+1:n*30,1:30)=temp121(:,:,n);
    im_com((n-1)*30+1:n*30,31:60)=temps(:,:,n);
    im_com((n-1)*30+1:n*30,61:90)=tempsas(:,:,n);
    im_com((n-1)*30+1:n*30,91:120)=tempt2(:,:,n);
    im_com((n-1)*30+1:n*30,121:150)=tempst(:,:,n);
end
dsp(interp2(im_com,2),1)
save images_data_may07 temp*

%%%%%%%%%%%%%%%%%%%%%%%%%
%movie and images for ASC
load real_data_st121 imir121
load real_data_mapsASC_oim ysas
load ./real_dataUM/real8g_mapASC_s1t5 x
for g=1:8    
    for sli=1:64
        yt(:,:,sli,g)=rot90(flipud(x(:,:,sli,g)),1);
    end
end  
load ./real_dataUM/real8g_mapASC_s3t3 x
for g=1:8    
    for sli=1:64
        yst(:,:,sli,g)=rot90(flipud(x(:,:,sli,g)),1);
    end
end  
slice=26;
for g=1:8
    temp121(:,:,g)=imir121(52:81,59:88,slice,g);
    temps(:,:,g)=ysas(52:81,59:88,slice,g);
    tempt(:,:,g)=yt(52:81,59:88,slice,g);
    tempst(:,:,g)=yst(52:81,59:88,slice,g);
end
temp121=temp121/max(temp121(:));
temps=temps/max(temps(:));
tempst=tempst/max(tempst(:));
tempt=tempt/max(tempt(:));
temp1=cat(2,temp121,temps);
temp2=cat(2,tempt,tempst);
temp=cat(1,temp1,temp2);
for g=1:8
    figure('position',[100 100 600 650]), imagesc(interp2(temp(:,:,g),2),[0 1]), clinicalcolor, axis image
    axis off
    mov(g)=getframe;pause(.1);close
end
figure('position',[100 100 600 650]),movie(mov,100)
movie2avi(mov,'real_st121_4DASC')

for n=1:8
    im_com((n-1)*30+1:n*30,1:30)=temp121(:,:,n);
    im_com((n-1)*30+1:n*30,31:60)=temps(:,:,n);
    im_com((n-1)*30+1:n*30,61:90)=tempt(:,:,n);
    im_com((n-1)*30+1:n*30,91:120)=tempst(:,:,n);
end
figure('position',[100 100 400 850]), imagesc(interp2(im_com,2)),clinicalcolor,axis image,axis off

%%%%%%%%%%%%%%%%
%image show
im_121=zeros(30*2,30*4);
for n=1:2
    for m=1:4
        im_121((n-1)*30+1:n*30,(m-1)*30+1:m*30)=temp121(:,:,(n-1)*4+m);
    end
end
figure('position',[100 100 800 450]),imagesc(interp2(im_121,2)),clinicalcolor,axis image,axis off


load ./real_dataUM/real8g_mapASC_s2t1 x
load ./real_dataUM/real8g_mapASC_s2t2 x
for s=1:5
    for t=1:5
        l_name=['./real_dataUM/real8g_mapASC_s' num2str(s) 't' num2str(t)];
    load(l_name)
yst=zeros(128,128,64,8);
for g=1:8    
    for sli=1:64
        yst(:,:,sli,g)=rot90(flipud(x(:,:,sli,g)),1);
    end
end   
slice=26;
for g=1:8
    tempscst(:,:,g)=yst(52:81,59:88,slice,g);
end

im_scst=zeros(30*2,30*4);
for n=1:2
    for m=1:4
        im_scst((n-1)*30+1:n*30,(m-1)*30+1:m*30)=tempscst(:,:,(n-1)*4+m);
    end
end
figure('position',[100 100 800 450]),imagesc(interp2(im_scst,2)),clinicalcolor,axis image,axis off

    end
end

tempscst=tempscst/max(tempscst(:));
tempscst=cat(2,tempscst,tempscst2);
for n=1:8
    imagesc(interp2(tempscst(:,:,n),2),[0 1]),clinicalcolor,axis image,axis off
    mov2(n)=getframe;pause(0.1);
end
figure,movie(mov,100)
figure,movie(mov2,100)
